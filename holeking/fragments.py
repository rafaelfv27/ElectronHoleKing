"""Fragment definitions for interfragment charge-transfer analysis.

Either a Multiwfn-style text file (one fragment per line, 1-based, ranges allowed: ``1,3,6-10``)
or a core SMARTS matched on the xyz geometry: the core (with its hydrogens) is one fragment and
each connected piece of the remainder is another.
"""
import numpy as np

# 2,1,3-benzothiadiazole ring atoms, written without bond orders so that perception is not needed.
CORE_SMARTS = {"btd": "[#6]1~[#6]~[#6]~[#6]2~[#7]~[#16]~[#7]~[#6]~2~[#6]~1"}


def read_fragment_file(path):
    frags = []
    for ln in open(path):
        ln = ln.split("#")[0].strip()
        if not ln:
            continue
        atoms = []
        for part in ln.replace(" ", "").split(","):
            if "-" in part:
                a, b = part.split("-")
                atoms += range(int(a) - 1, int(b))
            else:
                atoms.append(int(part) - 1)
        frags.append(sorted(atoms))
    return frags


def core_fragments(xyz_path, smarts, charge=0):
    """Return [core, piece1, piece2, ...] as lists of 0-based atom indices."""
    from rdkit import Chem
    from rdkit.Chem import rdDetermineBonds

    mol = Chem.MolFromXYZFile(xyz_path)
    rdDetermineBonds.DetermineConnectivity(mol, charge=charge)
    patt = Chem.MolFromSmarts(CORE_SMARTS.get(smarts, smarts))
    matches = mol.GetSubstructMatches(patt)
    if len(matches) != 1:
        raise ValueError(f"core pattern matched {len(matches)} times in {xyz_path}; expected 1")
    core = set(matches[0])
    core |= {n.GetIdx() for a in matches[0] for n in mol.GetAtomWithIdx(a).GetNeighbors()
             if n.GetAtomicNum() == 1}
    rest = [a.GetIdx() for a in mol.GetAtoms() if a.GetIdx() not in core]
    pieces, seen = [], set()
    for start in rest:
        if start in seen:
            continue
        stack, piece = [start], []
        seen.add(start)
        while stack:
            a = stack.pop()
            piece.append(a)
            for n in mol.GetAtomWithIdx(a).GetNeighbors():
                if n.GetIdx() not in core and n.GetIdx() not in seen:
                    seen.add(n.GetIdx())
                    stack.append(n.GetIdx())
        pieces.append(sorted(piece))
    return [sorted(core)] + pieces


def check_partition(frags, natm):
    flat = np.concatenate([np.asarray(f, dtype=int) for f in frags])
    if len(flat) != natm or len(set(flat.tolist())) != natm:
        raise ValueError(f"fragments must cover each of the {natm} atoms exactly once")
