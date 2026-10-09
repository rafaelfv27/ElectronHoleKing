"""Run xtb4stda and stda for one geometry.

xtb4stda is run twice on the same input: once with -molden (ground-state orbitals in a format
we can read) and once without (the binary wfn.xtb that stda needs). The orbitals of both runs
must agree including phases; Analyzer checks this per state by recomputing stda's oscillator
strength from the vectors.

Executables: $XTB4STDA_BIN and $STDA_BIN (default: on PATH). Parameter files (.param_stda*.xtb,
.param_gbsa_<solvent>): $XTB4STDAHOME. xtb4stda truncates that path at 80 characters, so the
needed files are copied into the work directory and XTB4STDAHOME is set to ".".
"""
import os
import shutil
import subprocess
from pathlib import Path

SOLVENTS = ("acetone", "acetonitrile", "benzene", "ch2cl2", "chcl3", "cs2", "dmso", "ether",
            "h2o", "methanol", "thf", "toluene")


def _param_files(solvent):
    home = Path(os.environ.get("XTB4STDAHOME", Path.home()))
    names = [".param_stda1.xtb", ".param_stda2.xtb"]
    if solvent:
        names.append(f".param_gbsa_{solvent}")
    files = [home / n for n in names]
    missing = [str(f) for f in files if not f.is_file()]
    if missing:
        raise FileNotFoundError(f"xtb4stda parameter files not found (set XTB4STDAHOME): {missing}")
    return files


def run(xyz, workdir, solvent=None, rpa=False, nstates=10, ethr=7.0, charge=0, nthreads=1):
    """Return the paths needed by the analysis: gs.molden, vector file and tda.dat."""
    if solvent and solvent not in SOLVENTS:
        raise ValueError(f"solvent {solvent!r} not in xtb4stda GBSA set {SOLVENTS}")
    xtb4stda = os.environ.get("XTB4STDA_BIN", "xtb4stda")
    stda = os.environ.get("STDA_BIN", "stda")
    work = Path(workdir)
    work.mkdir(parents=True, exist_ok=True)
    for f in _param_files(solvent):
        shutil.copy(f, work / f.name)
    shutil.copy(xyz, work / "geom.xyz")
    env = dict(os.environ, XTB4STDAHOME=".", OMP_NUM_THREADS=str(nthreads), MKL_NUM_THREADS=str(nthreads))
    base = [xtb4stda, "geom.xyz", "-chrg", str(charge)] + (["-gbsa", solvent] if solvent else [])

    def call(cmd, log):
        with open(work / log, "w") as fh:
            subprocess.run(cmd, cwd=work, env=env, stdout=fh, stderr=subprocess.STDOUT, check=True)

    call(base + ["-molden"], "xtb4stda_molden.out")
    os.replace(work / "molden.input", work / "gs.molden")
    call(base, "xtb4stda.out")
    call([stda, "-xtb", "-e", str(ethr), "-vectm", str(nstates)] + (["-rpa"] if rpa else []), "stda.out")
    vectors = work / ("sing_a" if rpa else "ciss_a")
    for p in (work / "gs.molden", vectors, work / "tda.dat"):
        if not p.is_file():
            raise RuntimeError(f"{p} was not written; see the logs in {work}")
    return work / "gs.molden", vectors, work / "tda.dat"
