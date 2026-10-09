"""Hole-electron analysis of one excitation (definitions of Multiwfn, Lu & Chen).

hole     D_h = Co (X X^T - Y Y^T) Co^T
electron D_e = Cv (X^T X - Y^T Y) Cv^T
Centroids and RMSDs are analytic (dipole and second-moment integrals); only Sr and the
Becke fragment populations need a grid, which is processed in blocks to bound memory.
"""
from dataclasses import dataclass, field

import numpy as np
from pyscf import dft

BOHR_ANG = 0.52917721092
GRID_BLOCK = 20000
NATURAL_TOL = 1e-6   # dropped hole/electron weight allowed in the grid quantities


@dataclass
class Result:
    f_check: float              # length-gauge f recomputed from X+Y and the orbitals
    sr: float
    d_index: float              # Angstrom
    d_vector: np.ndarray        # electron centroid - hole centroid, Angstrom
    h_ct: float
    h_index: float
    t_index: float
    sigma_hole: np.ndarray      # Angstrom, x/y/z
    sigma_ele: np.ndarray
    hole_integral: float
    ele_integral: float
    hole_negative: float        # integral of the negative part of the hole (sTD-DFT only)
    frag_hole: dict = field(default_factory=dict)   # partition -> fragment fractions
    frag_ele: dict = field(default_factory=dict)
    grid_discarded: float = 0.0  # hole/electron weight left out of Sr and Becke populations


def density_matrices(gs, exc):
    x, y = exc.x, exc.y
    d_hole = gs.c_occ @ (x @ x.T - y @ y.T) @ gs.c_occ.T
    d_ele = gs.c_vir @ (x.T @ x - y.T @ y) @ gs.c_vir.T
    return d_hole, d_ele


def oscillator_strength(gs, exc, r_int):
    d_trans = np.sqrt(2) * gs.c_occ @ (exc.x + exc.y) @ gs.c_vir.T
    mu = np.einsum("xpq,pq->x", r_int, d_trans)
    return 2 / 3 * exc.energy_au * float(mu @ mu)


def _moments(dm, r_int, rr_int):
    n = float(np.einsum("pq,pq", dm, rr_int[3]))  # overlap integral stored as index 3
    c = np.einsum("xpq,pq->x", r_int, dm) / n
    second = np.array([np.einsum("pq,pq", rr_int[k], dm) for k in range(3)]) / n
    return n, c, np.sqrt(np.clip(second - c ** 2, 0, None))


class Analyzer:
    def __init__(self, gs, grid_level=3, fragments=None):
        """fragments: list of lists of 0-based atom indices (optional)."""
        self.gs = gs
        mol = gs.mol
        self.r_int = mol.intor("int1e_r")
        rr = mol.intor("int1e_rr").reshape(3, 3, mol.nao, mol.nao)
        self.rr_int = np.array([rr[0, 0], rr[1, 1], rr[2, 2], mol.intor("int1e_ovlp")])
        self.grids = dft.gen_grid.Grids(mol)
        self.grids.level = grid_level
        self.grids.build()
        keep = self.grids.atm_idx >= 0   # pyscf pads the grid with zero-weight points (atm_idx -1)
        self.grids.coords = self.grids.coords[keep]
        self.grids.weights = self.grids.weights[keep]
        self.grids.atm_idx = self.grids.atm_idx[keep]
        self.fragments = fragments
        ao_atom = np.zeros(mol.nao, dtype=int)
        for a, (_, _, p0, p1) in enumerate(mol.aoslice_by_atom()):
            ao_atom[p0:p1] = a
        self.ao_atom = ao_atom

    def _mulliken_atoms(self, dm):
        pop = np.diag(dm @ self.rr_int[3])
        return np.bincount(self.ao_atom, weights=pop, minlength=self.gs.mol.natm)

    def analyze(self, excitations):
        """Analyze a list of excitations; the grid is traversed once for all of them."""
        mol, g = self.gs.mol, self.grids
        results, natural = [], []
        for exc in excitations:
            d_hole, d_ele = density_matrices(self.gs, exc)
            nh, ch, sh = _moments(d_hole, self.r_int, self.rr_int)
            ne, ce, se = _moments(d_ele, self.r_int, self.rr_int)
            dvec = (ce - ch) * BOHR_ANG
            d = float(np.linalg.norm(dvec))
            u = dvec / d if d > 1e-10 else np.zeros(3)
            hvec = (sh + se) / 2 * BOHR_ANG
            h_ct = float(np.linalg.norm(hvec * u))   # Multiwfn: norm of the component-wise product
            res = Result(oscillator_strength(self.gs, exc, self.r_int), 0.0, d, dvec, h_ct,
                         float(np.linalg.norm(hvec)), d - h_ct, sh * BOHR_ANG, se * BOHR_ANG,
                         nh, ne, 0.0)
            if self.fragments:
                res.frag_hole["mulliken"] = self._frag(self._mulliken_atoms(d_hole))
                res.frag_ele["mulliken"] = self._frag(self._mulliken_atoms(d_ele))
            x, y = exc.x, exc.y
            hole_c, hole_w, cut_h = _natural(x @ x.T - y @ y.T, self.gs.c_occ)
            ele_c, ele_w, cut_e = _natural(x.T @ x - y.T @ y, self.gs.c_vir)
            res.grid_discarded = max(cut_h, cut_e)
            natural.append((hole_c, hole_w, ele_c, ele_w))
            results.append(res)
        atom_h = np.zeros((len(results), mol.natm))
        atom_e = np.zeros((len(results), mol.natm))
        for p0 in range(0, len(g.weights), GRID_BLOCK):
            sl = slice(p0, p0 + GRID_BLOCK)
            ao = mol.eval_gto("GTOval", g.coords[sl])
            w, idx = g.weights[sl], g.atm_idx[sl]
            for k, (hole_c, hole_w, ele_c, ele_w) in enumerate(natural):
                rho_h = (ao @ hole_c) ** 2 @ hole_w
                rho_e = (ao @ ele_c) ** 2 @ ele_w
                results[k].sr += float(w @ np.sqrt(np.clip(rho_h, 0, None) * np.clip(rho_e, 0, None)))
                results[k].hole_negative += float(w @ np.minimum(rho_h, 0))
                atom_h[k] += np.bincount(idx, weights=w * rho_h, minlength=mol.natm)
                atom_e[k] += np.bincount(idx, weights=w * rho_e, minlength=mol.natm)
        if self.fragments:
            for k, res in enumerate(results):
                res.frag_hole["becke"] = self._frag(atom_h[k])
                res.frag_ele["becke"] = self._frag(atom_e[k])
        return results

    def _frag(self, atom_pop):
        return np.array([atom_pop[f].sum() for f in self.fragments]) / atom_pop.sum()


def _natural(p_mo, c, tol=NATURAL_TOL):
    """Eigenvectors of a hole/electron MO density matrix, dropping the smallest weights whose
    absolute sum stays below tol; returns AO coefficients, weights and the dropped weight."""
    val, vec = np.linalg.eigh(p_mo)
    order = np.argsort(-np.abs(val))
    val, vec = val[order], vec[:, order]
    tail = np.cumsum(np.abs(val)[::-1])[::-1]       # tail[k] = sum of |val[k:]|
    keep = int(np.searchsorted(-tail, -tol))           # first k with tail[k] < tol
    keep = max(keep, 1)
    return c @ vec[:, :keep], val[:keep], float(tail[keep]) if keep < len(val) else 0.0


def ifct(frag_hole, frag_ele):
    """Interfragment charge transfer (Multiwfn IFCT): Q[R,S] = hole_R * ele_S.

    Returns Q, net population change per fragment (ele - hole) and intrinsic CT fraction."""
    q = np.outer(frag_hole, frag_ele)
    return q, frag_ele - frag_hole, float(q.sum() - np.trace(q))
