"""Tests on 2,1,3-benzothiadiazole (sTDA-xTB and sTD-DFT-xTB, stda 1.6.3, gas phase).

Multiwfn references were obtained by feeding the same X/Y (written as a Gaussian TD output)
and the same gs.molden to Multiwfn 2026.7.15, main function 18 -> 1, medium grid.
"""
import os
from pathlib import Path

import numpy as np
import pytest
from pyscf.tools import molden

from holeking import analysis, fragments, io

DATA = Path(__file__).parent / "data"
BTD = DATA / "btd"
# (vector file, tda.dat, [(Sr, D) per state from Multiwfn])
CASES = [("ciss_a", "tda.dat", [(0.89917, 0.056), (0.89884, 0.871)]),
         ("rpa/sing_a", "rpa/tda.dat", [(0.85894, 0.043), (0.87884, 0.871)])]


@pytest.fixture(scope="module")
def gs():
    return io.read_xtb4stda_molden(BTD / "gs.molden")


@pytest.fixture(scope="module")
def analyzer(gs):
    return analysis.Analyzer(gs, fragments=[[3, 4, 5, 6, 7], [0, 1, 2, 8, 9, 10, 11, 12]])


def test_raw_pyscf_reading_is_wrong_and_ours_is_orthonormal(gs):
    mol, _, c, _, _, _ = molden.load(str(BTD / "gs.molden"))
    raw = np.abs(c.T @ mol.intor("int1e_ovlp") @ c - np.eye(c.shape[1])).max()
    assert raw > 1.0
    ours = np.abs(np.hstack([gs.c_occ, gs.c_vir]).T @ gs.mol.intor("int1e_ovlp")
                  @ np.hstack([gs.c_occ, gs.c_vir]) - np.eye(gs.c_occ.shape[1] + gs.c_vir.shape[1])).max()
    assert ours < 1e-6


@pytest.mark.parametrize("vec, dat, ref", CASES)
def test_against_stda_and_multiwfn(gs, analyzer, vec, dat, ref):
    excs = io.read_stda_vectors(BTD / vec, gs.c_occ.shape[1], gs.c_vir.shape[1])
    tda = io.read_tda_dat(BTD / dat)
    results = analyzer.analyze(excs)
    for (ev, f), r in zip(tda, results):
        assert r.f_check == pytest.approx(f, abs=2e-4)          # layout, X+Y and phases
        # stda writes real*4 vectors: sum(X^2 - Y^2) = 1 only to ~1e-5
        assert r.hole_integral == pytest.approx(1, abs=1e-4)
        assert r.ele_integral == pytest.approx(1, abs=1e-4)
    for (sr, d), r in zip(ref, results):
        assert r.sr == pytest.approx(sr, abs=1e-3)
        assert r.d_index == pytest.approx(d, abs=0.03)          # Multiwfn uses a uniform grid


def test_analytic_centroids_match_grid(gs, analyzer):
    exc = io.read_stda_vectors(BTD / "ciss_a", gs.c_occ.shape[1], gs.c_vir.shape[1])[1]
    d_hole, d_ele = analysis.density_matrices(gs, exc)
    g = analyzer.grids
    ao = gs.mol.eval_gto("GTOval", g.coords)
    cents = []
    for dm in (d_hole, d_ele):
        rho = np.einsum("gp,pq,gq->g", ao, dm, ao) * g.weights
        cents.append(rho @ g.coords / rho.sum())
    d_grid = np.linalg.norm(cents[1] - cents[0]) * analysis.BOHR_ANG
    assert analyzer.analyze([exc])[0].d_index == pytest.approx(d_grid, abs=1e-3)


def test_fragment_populations_and_ifct(gs, analyzer):
    exc = io.read_stda_vectors(BTD / "ciss_a", gs.c_occ.shape[1], gs.c_vir.shape[1])[2]
    r = analyzer.analyze([exc])[0]
    for part in ("mulliken", "becke"):
        assert r.frag_hole[part].sum() == pytest.approx(1)
        q, dpop, ct = analysis.ifct(r.frag_hole[part], r.frag_ele[part])
        assert q.sum() == pytest.approx(1)
        assert dpop.sum() == pytest.approx(0, abs=1e-12)
    assert r.frag_hole["mulliken"][0] > 0.95                   # n -> pi* on the thiadiazole


def test_fragment_file(tmp_path):
    f = tmp_path / "frags.txt"
    f.write_text("1,3-5\n2,6\n")
    assert fragments.read_fragment_file(f) == [[0, 2, 3, 4], [1, 5]]


def test_btd_core_split():
    frags = fragments.core_fragments(str(DATA / "dad_geom.xyz"), "btd")
    assert [len(f) for f in frags] == [11, 19, 19]             # core with 2 H, two arms
    single = fragments.core_fragments(str(BTD / "geom.xyz"), "btd")
    assert len(single) == 1 and len(single[0]) == 13


def test_wrong_orbital_count_is_rejected(gs):
    with pytest.raises(ValueError, match="does not match"):
        io.read_stda_vectors(BTD / "ciss_a", gs.c_occ.shape[1] - 1, gs.c_vir.shape[1])


def _csv_rows(path):
    import csv
    with open(path) as fh:
        return list(csv.DictReader(fh))


def test_analyze_threads_option(tmp_path):
    import os
    from threadpoolctl import threadpool_info
    from holeking import cli
    one, two = tmp_path / "t1.csv", tmp_path / "t2.csv"
    cli.main(["analyze", str(BTD), "--threads", "2", "--csv", str(two)])
    cli.main(["analyze", str(BTD), "--threads", "1", "--csv", str(one)])
    assert all(i["num_threads"] == 1 for i in threadpool_info() if i["user_api"] == "blas")
    cli.limit_threads(os.cpu_count() or 1)
    a, b = _csv_rows(one), _csv_rows(two)
    assert len(a) == len(b) == 3
    for ra, rb in zip(a, b):
        for k in ("sr", "d_index", "t_index"):
            assert float(ra[k]) == pytest.approx(float(rb[k]), abs=1e-6)


def test_workdir_with_several_xyz_is_rejected():
    from holeking import cli
    with pytest.raises(SystemExit):
        cli.main(["run", "a.xyz", "b.xyz", "--workdir", "w"])


@pytest.mark.skipif(not (os.environ.get("XTB4STDA_BIN") and os.environ.get("STDA_BIN")),
                    reason="needs xtb4stda and stda ($XTB4STDA_BIN, $STDA_BIN, $XTB4STDAHOME)")
def test_run_several_xyz_in_parallel(tmp_path, monkeypatch):
    import shutil
    from holeking import cli
    for name in ("a.xyz", "b.xyz"):
        shutil.copy(BTD / "geom.xyz", tmp_path / name)
    monkeypatch.chdir(tmp_path)
    cli.main(["run", "a.xyz", "b.xyz", "--jobs", "2", "--threads", "1", "--nstates", "3", "--csv", "out.csv"])
    rows = _csv_rows(tmp_path / "out.csv")
    assert [r["xyz"] for r in rows] == ["a"] * 3 + ["b"] * 3
    for ra, rb in zip(rows[:3], rows[3:]):                       # same geometry in both jobs
        assert float(ra["energy_ev"]) == pytest.approx(float(rb["energy_ev"]), abs=1e-3)
        assert float(ra["sr"]) == pytest.approx(float(rb["sr"]), abs=1e-4)


def test_interactive_mode_builds_run_command(tmp_path, monkeypatch, capsys):
    from holeking import cli
    for name in ("a.xyz", "b.xyz"):
        (tmp_path / name).write_text("")
    monkeypatch.chdir(tmp_path)
    answers = iter(["", "salt water", "water", "rpa", "zero", "", "-1", "9", "2", "", "[#7]~[#16]~[#7]", "x", "2",
                    ""])
    monkeypatch.setattr("builtins.input", lambda prompt: next(answers))
    assert cli.ask() == ["run", "a.xyz", "b.xyz", "--solvent", "h2o", "--rpa", "--nstates", "10",
                         "--charge", "-1", "--core", "[#7]~[#16]~[#7]", "--jobs", "2", "--csv", "results.csv"]
    out = capsys.readouterr().out
    assert "unknown solvent" in out and "type 1, 2, 3 or 4" in out and "no default" in out
