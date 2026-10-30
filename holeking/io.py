"""Readers for xtb4stda/stda files.

- ``xtb4stda -molden`` ground-state molden: its contraction coefficients already include the
  primitive normalization, while pyscf expects coefficients of normalized primitives.
- ``stda -vectm`` eigenvectors (TURBOMOLE format): ``ciss_a`` (sTDA, X) or ``sing_a``
  (sTD-DFT, X+Y followed by X-Y), occupied index outer, normalized to sum(X^2 - Y^2) = 1.
- ``tda.dat``: excitation energies (eV) and length-gauge oscillator strengths.
"""
import math
import os
import tempfile
from dataclasses import dataclass

import numpy as np
from pyscf.tools import molden

_L = {"s": 0, "p": 1, "d": 2, "f": 3, "g": 4}
ORTHO_TOL = 1e-4


def _primitive_norm(alpha, l):
    double_fact = math.prod(range(2 * l - 1, 0, -2))
    return (2 * alpha / math.pi) ** 0.75 * (4 * alpha) ** (l / 2) / math.sqrt(double_fact)


@dataclass
class GroundState:
    mol: object          # pyscf Mole
    c_occ: np.ndarray    # AO x nocc
    c_vir: np.ndarray    # AO x nvir


def read_xtb4stda_molden(path):
    """Load the ground state written by ``xtb4stda -molden`` and check orthonormality."""
    lines = open(path).read().splitlines()
    out, i, in_gto = [], 0, False
    while i < len(lines):
        ln = lines[i]
        if ln.strip().startswith("["):
            in_gto = ln.strip().upper().startswith("[GTO]")
            if ln.strip().upper().startswith("[TITLE]"):    # pyscf prints "Unknown section TITLE"
                i += 1
                continue
        tok = ln.split()
        if in_gto and len(tok) == 3 and tok[0].lower() in _L:
            l, nprim = _L[tok[0].lower()], int(tok[1])
            out.append(ln)
            for prim in lines[i + 1:i + 1 + nprim]:
                alpha, coef = (float(x.replace("D", "E")) for x in prim.split())
                out.append(f"  {alpha:.14E}  {coef / _primitive_norm(alpha, l):.14E}")
            i += 1 + nprim
            continue
        out.append(ln)
        i += 1
    fd, tmp = tempfile.mkstemp(suffix=".molden")
    try:
        with os.fdopen(fd, "w") as fh:
            fh.write("\n".join(out) + "\n")
        mol, _, c, occ, _, _ = molden.load(tmp)
    finally:
        os.remove(tmp)
    if c.ndim == 3:
        raise ValueError(f"{path}: open-shell molden is not supported")
    err = np.abs(c.T @ mol.intor("int1e_ovlp") @ c - np.eye(c.shape[1])).max()
    if err > ORTHO_TOL:
        raise ValueError(f"{path}: orbitals not orthonormal (max|C^T S C - I| = {err:.1e}); "
                         "is this a molden written by xtb4stda -molden?")
    return GroundState(mol, c[:, occ > 0], c[:, occ == 0])


@dataclass
class Excitation:
    energy_au: float
    x: np.ndarray        # nocc x nvir
    y: np.ndarray        # nocc x nvir (zeros for sTDA)


def read_stda_vectors(path, nocc, nvir):
    """Read ``ciss_a`` or ``sing_a`` written by ``stda -vectm``."""
    states, cur, rpa = [], None, None
    for ln in open(path).read().splitlines():
        if ln.startswith("$scfinstab"):
            kind = ln.split()[1]
            if kind not in ("ciss", "rpas"):
                raise ValueError(f"{path}: unsupported vector type {kind!r} (singlets only)")
            rpa = kind == "rpas"
        if "eigenvalue" in ln:
            cur = []
            states.append((float(ln.split("=")[1].replace("D", "E")), cur))
        elif ln.startswith("$"):
            cur = None
        elif cur is not None:
            s = ln.rstrip()
            cur += [float(s[k:k + 20].replace("D", "E")) for k in range(0, len(s), 20)]
    if rpa is None:
        raise ValueError(f"{path}: no $scfinstab line, not a stda vector file")
    dim = nocc * nvir
    result = []
    for energy, values in states:
        v = np.array(values)
        if len(v) != (2 * dim if rpa else dim):
            raise ValueError(f"{path}: vector length {len(v)} does not match nocc*nvir = {dim} "
                             "(molden and vectors from different runs?)")
        if rpa:
            plus, minus = v[:dim].reshape(nocc, nvir), v[dim:].reshape(nocc, nvir)
            x, y = (plus + minus) / 2, (plus - minus) / 2
        else:
            x, y = v.reshape(nocc, nvir), np.zeros((nocc, nvir))
        result.append(Excitation(energy, x, y))
    return result


def read_tda_dat(path):
    """Return a list of (energy_eV, f_length) from stda's tda.dat."""
    block = open(path).read().split("DATXY")[1]
    rows = [ln.split() for ln in block.strip().splitlines() if ln.strip()]
    return [(float(r[1]), float(r[2])) for r in rows]
