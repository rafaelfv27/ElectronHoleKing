# HoleKing

Hole–electron and interfragment charge-transfer (IFCT) analysis of excitations computed with
sTDA-xTB or sTD-DFT-xTB (`xtb4stda` + `stda`/`std2`). Multiwfn cannot read `stda` output for this
analysis, and ORCA's `XTB` keyword does not write a usable wavefunction; this package works from
the full `stda` eigenvectors instead.

For each excited state it reports:

- Sr, D, H, H_CT and t indices (definitions of Multiwfn's hole–electron analysis);
- hole and electron populations on user fragments (Mulliken and Becke partitions) and the IFCT
  quantities (net population change per fragment, intrinsic CT);
- the length-gauge oscillator strength recomputed from the vectors, which must match `stda`'s.

## How it works

1. `xtb4stda geom.xyz -molden` gives the ground-state orbitals; `xtb4stda geom.xyz` gives the
   binary `wfn.xtb` that `stda` reads. Both runs give the same orbitals, phases included (checked
   for every state through the oscillator strength).
2. `stda -xtb -vectm N [-rpa]` writes the eigenvectors (`ciss_a` = X; `sing_a` = X+Y then X−Y).
3. Hole = Co(XXᵀ − YYᵀ)Coᵀ, electron = Cv(XᵀX − YᵀY)Cvᵀ. Centroids and spreads (D, H, t) come
   from dipole and second-moment integrals, so they are exact; Sr and the Becke populations are
   integrated on a pyscf Becke grid, using the eigenvectors of the hole/electron matrices
   (weights below 1e-6 dropped, reported as `grid_discarded`).

## Requirements

- Python ≥ 3.9, numpy, pyscf ≥ 2.3; rdkit only for `--core`.
- `xtb4stda` and `stda` (1.6.1, 1.6.3) or `std2` (2.0.1) — all three give identical results.
  `$XTB4STDA_BIN`, `$STDA_BIN` point to the executables, `$XTB4STDAHOME` to the folder with
  `.param_stda1.xtb`, `.param_stda2.xtb` and, for implicit solvent, `.param_gbsa_<solvent>`.

## Usage

```bash
holeking run geom.xyz --nstates 10 --core btd --csv out.csv            # sTDA-xTB, gas phase
holeking run geom.xyz --rpa --solvent h2o --fragments frags.txt        # sTD-DFT-xTB, GBSA water
holeking analyze geom_holeking/ --core "[#6]1~[#6]~[#6]~[#6]2~[#7]~[#16]~[#7]~[#6]~2~[#6]~1"
holeking run mol1.xyz mol2.xyz mol3.xyz --jobs 8 --threads 1 --core btd --csv all.csv   # batch
```

With several xyz files each geometry gets its own `<xyz stem>_holeking` work directory and the CSV gains an
`xyz` column.

### Choosing jobs, threads and grid level

- `--jobs N` (`run` only): geometries processed in parallel, one process each. Default 1.
- `--threads T`: threads per job for `xtb4stda`, `stda` and the analysis (OpenBLAS/MKL and pyscf's OpenMP).
  Default: number of CPUs // jobs. Keep `jobs × threads` at or below the number of logical CPUs.
- `--grid-level L`: pyscf Becke grid used for Sr and the Becke populations (default 3). D, H, t, centroids and the
  Mulliken populations are analytic and do not depend on it.

Measured on an i5-10300H (4 cores / 8 threads), `holeking analyze`, 20 to 60 sTDA-xTB work dirs of 25–72 atoms with
10 states each, against level 3:

| setting | time | largest change |
|---|---|---|
| level 3, 4 jobs × 2 threads (60 dirs) | 375 s | — |
| level 2, 4 jobs × 2 threads (60 dirs) | 276 s | Sr 4.8e-4, Becke fraction 4.8e-5 |
| level 1, 4 jobs × 2 threads (60 dirs) | 142 s | Sr 2.4e-3 (20 of 60 dirs above 1e-3) |
| level 3, 8 × 1 / 4 × 2 / 4 × 1 / 2 × 4 (20 dirs) | 117 / 129–133 / 132 / 170 s | identical results |

Level 1 is not recommended: its Sr error exceeds the 1e-3 agreement with Multiwfn used for validation. On this CPU
one job per logical CPU with one thread each was fastest; hyperthreads did not speed up the BLAS part.

Fragments: a Multiwfn-style file (one fragment per line, 1-based atoms, ranges allowed:
`1,3,6-10`), or `--core SMARTS|btd`, which takes the matched core plus its hydrogens as fragment 1
and every connected piece of the remainder as a further fragment.

Output columns include `first_bright` (lowest state with f ≥ 5 % of the largest f among the
computed states) and `near_degenerate_bright` (another bright state within 0.1 eV).

## Validation (`tests/`, plus the checks below)

| check | result |
|---|---|
| f recomputed from vectors vs `stda` | equal for every state tested (TDA and RPA, gas and GBSA, up to 135 atoms) |
| same X/Y analyzed by Multiwfn 2026.7.15 | Sr within 1e-3; D within 0.03 Å (Multiwfn's uniform grid; our D is analytic) |
| same index definitions on CAM-B3LYP TDA NTOs (ORCA) vs Multiwfn's own analysis | D, Sr, t within 0.01 |
| IFCT vs Multiwfn (Mulliken-like) | fragment fractions within 0.4 percentage points |
| `stda` 1.6.1 / 1.6.3 / `std2` 2.0.1 | identical to 1e-5 |
| 135-atom molecule, 10 states | 2 min 44 s, 0.75 GB |

## Known limitations (measured)

- **Symmetric D–A–D molecules:** D and t do not see the charge transfer (hole and electron
  centroids nearly coincide). Use the fragment populations / IFCT, e.g. `--core btd`.
- **Oscillator strengths:** sTDA-xTB gave about twice the CAM-B3LYP TDA value for CT states
  (0.80 vs 0.38; 1.56 vs 0.70) while sTDA on CAM-B3LYP orbitals did not; the inflation comes from
  the xTB part. Calibrate or use f only for ranking.
- **GBSA:** only the ground-state orbitals feel the solvent. Toluene already gives most of the
  water effect (Born-type 1 − 1/ε); linear-response CPCM with CAM-B3LYP also underestimates CT
  solvatochromism. Do not over-interpret differences between solvents.
- **State order:** near-degenerate bright states may swap relative to TD-DFT (flagged).
- Singlets of closed-shell molecules only.
- `stda` vectors are single precision: hole/electron integrals are 1 to about 1e-5.

## References

- T. Lu, F. Chen, *J. Comput. Chem.* **33**, 580 (2012) — Multiwfn; hole–electron and IFCT in its manual.
- S. Grimme, C. Bannwarth, *J. Chem. Phys.* **145**, 054103 (2016) — sTDA-xTB.
- std2: https://github.com/Theoretical-Chemistry-Group-UCLouvain/std2
- Q. Sun et al., *J. Chem. Phys.* **153**, 024109 (2020) — PySCF.

## License

MIT, see `LICENSE`.
