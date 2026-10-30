# HoleKing

HoleKing tells you where the hole and the electron go in each excited state of a molecule.
It works with sTDA-xTB and sTD-DFT-xTB, which are fast enough for large dye molecules. Multiwfn
can't read these results, so HoleKing does the hole–electron analysis itself.

For every excited state you get:

- energy, wavelength and oscillator strength;
- the usual hole–electron indices (Sr, D, H, t), same definitions as Multiwfn;
- how much of the hole and of the electron sits on each fragment, and how much charge moves
  from one fragment to another (IFCT).

Everything ends up in one CSV file, one row per state.

## Requirements

- Python 3.9 or newer
- `xtb4stda` and `stda` (or `std2`)
- rdkit, only if you want fragments split around a core automatically

## Tutorial: your first run

A printable version of this tutorial, with more detail, is in
[`docs/how-to-run.pdf`](docs/how-to-run.pdf).

**1. Tell HoleKing where xtb4stda and stda are.** Put these lines in your `~/.bashrc` so you
only do it once:

```bash
export XTB4STDA_BIN=/path/to/xtb4stda
export STDA_BIN=/path/to/stda
export XTB4STDAHOME=/path/to/folder/with/.param_stda1.xtb
```

The `XTB4STDAHOME` folder must hold `.param_stda1.xtb`, `.param_stda2.xtb` and, if you use a
solvent, `.param_gbsa_<solvent>`.

**2. Install HoleKing.** From the folder you downloaded:

```bash
pip install ".[core]"      # or: pip install .   (without rdkit)
```

**3. Run it.** Open a terminal in the folder with your `.xyz` files and type `holeking`.
It asks a few questions. Press Enter to keep the value in brackets:

```
$ holeking
xyz file(s) [dad.xyz]:
solvent (none, water, toluene, thf, ...) [none]: water
method: tda = sTDA-xTB, rpa = sTD-DFT-xTB [tda]:
number of excited states [10]:
molecular charge [0]:
fragments: 1 = BTD core, 2 = other core (SMARTS), 3 = fragment file, 4 = none [1]:
save results to [dad.csv]:
```

- **method:** `tda` (sTDA-xTB) is the fast default; `rpa` (sTD-DFT-xTB) is a bit slower.
- **fragments:** how to split the molecule into donor/acceptor parts. With a BTD core (or any
  core you describe with SMARTS), the core and its hydrogens are fragment 1 and each piece
  attached to it becomes its own fragment. Or write a fragment file the Multiwfn way: one
  fragment per line, atom numbers starting at 1, ranges allowed (`1,3,6-10`).

If you type something it can't use, it tells you why and asks again. Ctrl+C quits.

**4. Look at the results.** HoleKing makes a folder `dad_holeking/` with the xtb and stda
files and logs, and writes `dad.csv`, one row per excited state:

```
state,energy_ev,wavelength_nm,f,first_bright,sr,d_index,t_index,...
1,4.5825,270.56,0.2645,1,0.8997,0.0453,-1.711,...
```

| column | meaning |
|---|---|
| `energy_ev`, `wavelength_nm`, `f` | excitation energy, wavelength, oscillator strength |
| `first_bright` | 1 for the lowest state with real intensity |
| `sr` | hole–electron overlap: near 1 is a local excitation, near 0 is charge transfer |
| `d_index` | distance between hole and electron centres (Å) |
| `t_index` | positive means hole and electron are well separated |
| `hole_fragN_*`, `ele_fragN_*` | fraction of the hole / electron on fragment N |
| `ct_intrinsic_*` | fraction of the excitation that is real charge transfer |

**5. Next time, skip the questions.** At the end HoleKing prints the same job as a one-line
command. Use it in scripts or on a cluster:

```bash
holeking run dad.xyz --solvent h2o --core btd --csv dad.csv
holeking run *.xyz --jobs 8 --core btd --csv all.csv       # 8 molecules at a time
holeking analyze dad_holeking/ --fragments frags.txt       # redo the analysis only
holeking run --help                                         # every option
```

The results were checked against Multiwfn (Sr within 0.001, D within 0.03 Å) and the oscillator
strengths are recomputed and compared with `stda` for every state. See `tests/`.

## References

- T. Lu, F. Chen, *J. Comput. Chem.* **33**, 580 (2012). Multiwfn.
- S. Grimme, C. Bannwarth, *J. Chem. Phys.* **145**, 054103 (2016). sTDA-xTB.
- std2: https://github.com/Theoretical-Chemistry-Group-UCLouvain/std2
- Q. Sun et al., *J. Chem. Phys.* **153**, 024109 (2020). PySCF.

MIT license.
