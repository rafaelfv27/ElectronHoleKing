"""Command line: ``holeking run geom.xyz [more.xyz ...]`` or ``holeking analyze workdir``."""
import argparse
import csv
import multiprocessing
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

from . import analysis, fragments, io, run

HC_EV_NM = 1239.84198
F_TOL_ABS, F_TOL_REL = 2e-3, 1e-2      # stda writes single-precision vectors
BRIGHT_FRACTION = 0.05                  # first bright state: f >= 5 % of the largest f
NEAR_DEGENERATE_EV = 0.1


def _fragments(args, xyz, natm):
    if args.fragments:
        frags = fragments.read_fragment_file(args.fragments)
    elif args.core:
        frags = fragments.core_fragments(xyz, args.core, charge=args.charge)
    else:
        return None
    fragments.check_partition(frags, natm)
    return frags


def analyze_dir(work, args):
    work = Path(work)
    vec = work / ("sing_a" if (work / "sing_a").is_file() else "ciss_a")
    gs = io.read_xtb4stda_molden(work / "gs.molden")
    frags = _fragments(args, str(work / "geom.xyz"), gs.mol.natm)
    excs = io.read_stda_vectors(vec, gs.c_occ.shape[1], gs.c_vir.shape[1])
    tda = io.read_tda_dat(work / "tda.dat")[:len(excs)]
    results = analysis.Analyzer(gs, grid_level=args.grid_level, fragments=frags).analyze(excs)
    for k, ((ev, f), r) in enumerate(zip(tda, results)):
        if abs(r.f_check - f) > F_TOL_ABS + F_TOL_REL * f:
            raise RuntimeError(f"state {k + 1}: f from vectors {r.f_check:.4f} != stda {f:.4f}; "
                               "the molden and the vectors do not belong to the same orbitals")
    fmax = max(f for _, f in tda)
    bright = [k for k, (_, f) in enumerate(tda) if f >= BRIGHT_FRACTION * fmax]
    rows = []
    for k, ((ev, f), r) in enumerate(zip(tda, results)):
        row = {"state": k + 1, "energy_ev": ev, "wavelength_nm": HC_EV_NM / ev, "f": f,
               "first_bright": int(bool(bright) and k == bright[0]),
               "near_degenerate_bright": int(k in bright and any(
                   j != k and abs(tda[j][0] - ev) < NEAR_DEGENERATE_EV for j in bright)),
               "sr": r.sr, "d_index": r.d_index, "h_ct": r.h_ct, "h_index": r.h_index,
               "t_index": r.t_index, "d_x": r.d_vector[0], "d_y": r.d_vector[1], "d_z": r.d_vector[2],
               "sigma_hole": float(np.linalg.norm(r.sigma_hole)),
               "sigma_ele": float(np.linalg.norm(r.sigma_ele)),
               "hole_negative": r.hole_negative, "grid_discarded": r.grid_discarded}
        for part in r.frag_hole:
            q, dpop, ct = analysis.ifct(r.frag_hole[part], r.frag_ele[part])
            for i in range(len(frags)):
                row[f"hole_frag{i + 1}_{part}"] = r.frag_hole[part][i]
                row[f"ele_frag{i + 1}_{part}"] = r.frag_ele[part][i]
                row[f"dpop_frag{i + 1}_{part}"] = dpop[i]
            row[f"ct_intrinsic_{part}"] = ct
        rows.append(row)
    return rows


def limit_threads(n):
    """Cap OpenBLAS/MKL (threadpoolctl) and pyscf's OpenMP at n threads in this process."""
    from pyscf import lib
    from threadpoolctl import threadpool_limits
    lib.num_threads(n)
    threadpool_limits(n)


def run_one(xyz, work, args):
    run.run(xyz, work, solvent=args.solvent, rpa=args.rpa, nstates=args.nstates, ethr=args.ethr,
            charge=args.charge, nthreads=args.threads)
    return analyze_dir(work, args)


def write_csv(rows, out):
    fh = open(out, "w", newline="") if out else sys.stdout
    w = csv.DictWriter(fh, fieldnames=list(rows[0]))
    w.writeheader()
    for row in rows:
        w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in row.items()})
    if out:
        fh.close()


def main(argv=None):
    ap = argparse.ArgumentParser(prog="holeking", description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p_run = sub.add_parser("run", help="run xtb4stda + stda, then analyze")
    p_run.add_argument("xyz", nargs="+")
    p_run.add_argument("--workdir", default=None, help="default: <xyz stem>_holeking (single xyz only)")
    p_run.add_argument("--solvent", choices=run.SOLVENTS)
    p_run.add_argument("--rpa", action="store_true", help="sTD-DFT-xTB instead of sTDA-xTB")
    p_run.add_argument("--nstates", type=int, default=10)
    p_run.add_argument("--ethr", type=float, default=7.0, help="stda -e energy threshold (eV)")
    p_run.add_argument("--jobs", type=int, default=1, help="geometries processed in parallel (default 1)")
    p_ana = sub.add_parser("analyze", help="analyze an existing work directory")
    p_ana.add_argument("workdir")
    p_ana.set_defaults(jobs=1)
    for p in (p_run, p_ana):
        p.add_argument("--threads", type=int, default=None,
                       help="threads per job for xtb4stda, stda and the analysis (default: CPUs // jobs)")
        p.add_argument("--charge", type=int, default=0)
        p.add_argument("--fragments", help="Multiwfn-style fragment file (1-based atoms)")
        p.add_argument("--core", help="core SMARTS or 'btd'; the rest is split into connected pieces")
        p.add_argument("--grid-level", type=int, default=3,
                       help="pyscf Becke grid level for Sr and Becke populations (default 3)")
        p.add_argument("--csv", help="output file (default: stdout)")
    args = ap.parse_args(argv)
    if args.jobs < 1 or (args.threads is not None and args.threads < 1):
        ap.error("--jobs and --threads must be >= 1")
    if args.threads is None:
        args.threads = max(1, (os.cpu_count() or 1) // args.jobs)
    if args.cmd == "analyze":
        limit_threads(args.threads)
        write_csv(analyze_dir(args.workdir, args), args.csv)
        return
    if len(args.xyz) > 1 and args.workdir:
        ap.error("--workdir takes a single xyz; with several, each goes to <xyz stem>_holeking")
    works = [args.workdir or f"{Path(x).stem}_holeking" for x in args.xyz]
    if len(set(works)) < len(works):
        ap.error("two xyz files share a stem, so they would share a work directory")
    if args.jobs == 1:
        limit_threads(args.threads)
        tables = [run_one(x, w, args) for x, w in zip(args.xyz, works)]
    else:
        with ProcessPoolExecutor(args.jobs, mp_context=multiprocessing.get_context("spawn"),
                                 initializer=limit_threads, initargs=(args.threads,)) as pool:
            tables = list(pool.map(run_one, args.xyz, works, [args] * len(works)))
    if len(args.xyz) == 1:
        rows = tables[0]
    else:
        rows = [{"xyz": Path(x).stem, **row} for x, t in zip(args.xyz, tables) for row in t]
    write_csv(rows, args.csv)


if __name__ == "__main__":
    main()
