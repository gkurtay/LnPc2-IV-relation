# Structure–spectrum relation for the intervalence band of lanthanide double-deckers: data, code and explorer

Data and software accompanying

> G. Kurtay, M. İnce, *Inter-ring distance and skew angle control the near-infrared intervalence band of lanthanide bis(phthalocyaninato) double-deckers* (in preparation; reference to be updated on publication).

Source code and Explorer: https://github.com/gkurtay/LnPc2-IV-relation (Explorer online: https://gkurtay.github.io/LnPc2-IV-relation/). Archived record, including all calculation files: Zenodo, https://doi.org/10.5281/zenodo.23839079.

## Contents

| Path | What it is |
|---|---|
| `explorer/LnPc2_Explorer_offline.html` | Interactive explorer. Opens in any modern browser, needs no network connection. Every computed and X-ray structure in 3D, descriptors recomputed in the browser, spectra, the relation, the experimental test against eight unsubstituted LnPc₂ and the NTO analysis. Your own ORCA outputs, Gaussian logs or XYZ files can be opened in it as well (see below). |
| `data/constants.json` | Constants of both forms of the relation (fitted to 11 donor-free structures), the donor correction, fit statistics and the experimental scale factors. |
| `data/structures.csv` | All 181 entries (170 computed structures and 11 X-ray entries): metal, set, functional, θ, d_N4, d_Cα, IV energy, oscillator strength, donor share D_sub, whether the structure enters the donor regression, deviation from both forms without and with the donor correction (computed structures above the node), ring spin asymmetry, ground-state ⟨S²⟩ and ⟨S²⟩ of the IV state (`S2_IV_state`). For the X-ray rows, IV_cm-1 is the measured solution band and the deviations are left empty. |
| `data/unsubstituted_LnPc2_test.csv` | Crystal descriptors of unsubstituted LnPc₂ (Pr–Yb) and the four IV components measured by Ostendorp and Homborg (CH₂Cl₂). |
| `data/geometry.csv` | Geometric parameters of every optimised structure, crystal structure and heteroleptic (Table S11): θ, d_N4, d_Cα, Ln–N mean and range, height of each N₄ plane above the metal, cis N–Ln–N angle, ring doming; computed with the Explorer descriptor routine. |
| `data/orbitals_a1u_LuPc2opt.txt` | Signed β HOMO and β LUMO values 0.75 Å above and below every C/N atom of the optimised LuPc₂ (Figure 1d), read from the ORCA cube files; `orbitals_a1u_LuPh8core.txt`, the same for the LuPc₂ skeleton of the Lu(Ph₈Pc)₂ crystal structure. |
| `data/aryl_optimised.csv` | Phenyl groups at the X-ray and at gas-phase optimised geometries (Table S9g): torsion, Pc–Ph bond, IV energy and core, shift from the core, per-group rule, shift from the relation, β gap. |
| `data/orbital_scan.csv` | Symmetry-labelled orbital scan of LuPc₂ (Figure S4): energy difference between the A and B1 combinations of the ring a₁u orbitals against the skew angle, for the anion and the neutral radical (α). |
| `data/heteroleptic.csv` | The 15 computed heteroleptic sandwiches: substituent and its class, descriptors, IV energy, ring spin asymmetry a, D_sub, share of the ground-state hole on the substituted ring, Δ_eff = |a|·ν, deviation from both forms raw, after the asymmetry correction 1/√(1 − a²) and after both corrections. |
| `data/aryl_series.csv` | Lu(Ph₈Pc)₂ at its X-ray geometry with the phenyl groups rotated (35–90°), halved, replaced by H (core) or by methyl: torsion, IV energy, shift from the core, α a₁u splitting, β HOMO–LUMO gap, electron–hole term K = ω − gap, and the phenyl share of the β HOMO and LUMO. |
| `data/solvent.csv` | IV energies with CPCM (chlorobenzene) and in the gas phase for seven structures, and their ratio (Figure S3). |
| `data/fragments.csv` | Substituent fragments (Ar–D on benzene, PBE0-D3BJ/def2-SVP): HOMO, vertical ionisation energy and the ring-HOMO: the highest occupied orbital with at least 25% of its weight on the benzene ring that stands for the benzo ring of Pc (energy, weight, its position below the HOMO), with the ring-weighted centroid of the eight highest occupied orbitals. D-fragments model the computed heteroleptics, A-fragments are candidate donors. |
| `data/nto_partition.csv` | Share of the IV hole and particle densities on the metal, each macrocycle, the substituents and the substituent heteroatoms, and D_sub. The row of the X-ray entry LuPh8_x comes from the NTO calculation at the same geometry (`calculations/nto/LuPh8_k_td2_nto`). |
| `data/LnPc2_IV_coordinates.xyz` | Cartesian coordinates (Å) of every computed structure as used in the TD-DFT calculations: the 170 of `data/structures.csv` and the 5 of `data/aryl_optimised.csv`; the comment line gives the structure code, label, functional, IV energy and descriptors. |
| `code/descriptors.py` | The descriptor routine used throughout the article (d_N4, d_Cα, θ). Reads CIF, XYZ or ORCA output. |
| `code/predict.py` | Predicts the IV band of a homoleptic double-decker from its geometry. |
| `code/unsubstituted_LnPc2_test.py` | Reproduces the external test of the article. |
| `code/explorer_images.js` | Exports the structure images of the article (Chart 1, Figures 1, 2 and 4–7, S6, S7) from the Explorer through its scripting interface, with the measured distances and angles; job list in `code/explorer_jobs.json`. Needs Node.js and Playwright. |
| `code/heteroleptic.py` | The relation with the asymmetry and donor corrections for heteroleptic sandwiches; `--table` reproduces the heteroleptic analysis. |
| `code/aryl.py` | Torsion dependence and gap/electron–hole decomposition of the aryl effect (`--series`); for a new structure, finds aryl rings on the benzo rings and estimates their shift additively. |
| `code/figstyle.py` | Common figure style used by `make_outputs.py` (palette, panel background, marker shading, image cards). |
| `code/make_outputs.py` | Regenerates Table 1, Figures 1–3, 6 and 7, Figures S1–S3 and S5–S8 of the Supporting Information and the analyses behind Tables S9c, S9f and S10b from the files in `data/` (writes `outputs/`). |
| `code/ring_asymmetry.py` | Ring spin asymmetry a from an ORCA output (Mulliken) or a Gaussian log with pop=hirshfeld; reproduces `ring_spin_asym_a` in `data/structures.csv` and the crystal-structure values in `calculations/hirshfeld/`. |
| `code/torsion_mode.py` | Inter-ring torsional mode from an ORCA Hessian: frequency, thermal amplitude and its effect on the angular factor (Section S3); run on `calculations/freq/*/*.hess`. |
| `code/nto_partition.py` | Partitions NTO cube files onto metal, rings and substituents. The cube files are not deposited; they are regenerated with orca_plot from the NTO calculations in `calculations/nto/`, and the results are in `data/nto_partition.csv`. |
| `outputs/` | The output of `make_outputs.py`, as deposited. |
| `calculations/` | ORCA 6.1.1 input and output files of all calculations (TD-DFT, natural transition orbitals, CASSCF/NEVPT2, aryl series, harmonic frequencies, basis-set check, orbital scan) and Gaussian 16 files of the substituent fragments and Hirshfeld calculations, with `index.csv` linking each job to its structure; see `calculations/README.txt`. In the Zenodo record only; not mirrored on GitHub because of its size. |

## Quick start

```
pip install numpy scipy matplotlib gemmi     # gemmi only for CIF input; scipy and matplotlib for make_outputs.py
python code/predict.py my_structure.cif
python code/predict.py my_structure.xyz --theta 45
python code/predict.py my_structure.xyz --dsub 0.02   # with heteroatom donors
python code/unsubstituted_LnPc2_test.py
python code/heteroleptic.py --table
python code/aryl.py --series
python code/make_outputs.py          # tables and figures of the article
```

`predict.py` prints both forms in the gas phase and multiplied by the experimental scale factors (Cα 0.710, N₄ 0.696), which convert the relation to solution band positions. For substituents with heteroatom donors pass the donor share D_sub of the IV transition (`--dsub`), obtained from the natural transition orbitals with `code/nto_partition.py`; the relation is then multiplied by (1 + c·D_sub) with c = -1.44 (Cα) or -0.91 (N₄).

## Reproducing the article

`python code/make_outputs.py` refits both forms to the 11 donor-free structures in `data/structures.csv` (with the bootstrap of the article, 2000 resamples, seed 7) and redraws Figures 1–3, 6, 7 and S1–S7, placing the deposited Explorer images of `outputs/explorer/`. Figures 4 and 5 combine Explorer images with literature spectra and are not regenerated. Panel letters only; captions are in the article.

## Your own structures in the Explorer

Open your own file in the library panel, or drop files anywhere on the page, to add a lanthanide double-decker of your own: an ORCA output (`.out`), a Gaussian log (`.log`) or an XYZ file. The file is read in the browser and is not sent anywhere. The last geometry in the file is used. The Explorer computes d_N4, d_Cα, θ and the ring doming with the routine that gives the published values, and both forms of the relation. When the file contains TD-DFT results it also shows the stick spectrum and picks the IV band: the strongest root between 1500 and 14000 cm⁻¹ polarised along the sandwich axis, or the strongest root in that window when the file has no transition dipoles; another root can be chosen. It reports the ring spin asymmetry a from the Mulliken spin populations. The structure can then be measured, overlaid on any published structure and exported, and it appears on the descriptor map and in the data table. Files are kept in the browser for the next visit until removed. For structures with heteroatom donors, the donor correction needs Δ_sub from `code/nto_partition.py`, which can be entered in the calculator.

## Scripting the Explorer

The Explorer exposes `window.LNPC` in the browser console (or to a headless browser):

```
LNPC.select('LuPc2bare_gas')          // structure by id (see data/structures.csv)
LNPC.ui({style:'ball', rings:true})    // display: style, h, rings (colour by ring), spinsph (spin spheres), half, num ...
LNPC.view('top'); LNPC.rotate('x',20); LNPC.zoom(0.9)
LNPC.meas('dist',[1,5])                // Ln–N distance, atom numbers 1-based as shown in the viewer
LNPC.meas('ang',[5,1,10]); LNPC.meas('dih',[196,193,141,145])
LNPC.desc()                            // theta, d_N4, d_Ca, metal and N8 atom numbers
LNPC.png('c2_600')                     // PNG as data URL
```

`node code/explorer_images.js` runs the job list `code/explorer_jobs.json` and writes `outputs/explorer/`; `make_outputs.py` places these images in the figures and reads the measured values from the accompanying .json files. The images used in the article are deposited in `outputs/explorer/`.

In the data files and code the donor share Δ_sub of the article is written D_sub, and the N₄ form is ν = (p + q·d_N4)(β − cos 4θ)/(1 + β).

## Descriptors

* **d_N4**: separation of the two N₄ centroids along the sandwich axis.
* **d_Cα**: d_N4 plus the mean height of the eight α-carbons of each ring above its own N₄ plane.
* **θ**: mean azimuthal offset between the coordinating nitrogens of the two rings, folded into 0–45°, averaged over the two choices of reference ring.

## Scope

The relation applies to homoleptic double-deckers with skew angles from about 25° (the lowest fitted point; the coupling node lies near 20°) to 45°. Heteroatom donors (thioether, amino) lower the band below the relation by up to about 7%, correlated with their share of the transition hole; the donor correction accounts for this. Heteroleptic sandwiches, whose rings are inequivalent, are outside its calibration; at computed geometries they follow it once the ring asymmetry (factor 1/√(1 − a²)) and the donor share are included (`code/heteroleptic.py`), but at the X-ray geometries of measured heteroleptics neither the relation nor TD-DFT reproduces the bands. Aryl groups lower the band by a route that D_sub does not register (`code/aryl.py`). In the test against unsubstituted LnPc2 the crystal structure of LaPc₂ (COD 7243588, 7243598) is not used, because its N₄–N₄ separation (2.82–2.86 Å) is shorter than that of NdPc₂, contrary to the lanthanide contraction.

## Methods in brief

Geometries: PBE0-D3BJ/def2-SVP. Excited states: CAM-B3LYP/def2-TZVP, RIJCOSX, unrestricted, ORCA 6.1.1. Natural transition orbitals at the same level, with the IV state identified by its energy (within 50 cm⁻¹ of the reference calculation, at least 12 roots). NTO cubes plotted from the transition orbitals without SCF iterations. CASSCF(3,2)/SC-NEVPT2 from quasi-restricted orbitals, with FIC-NEVPT2 and CAS(3,6) checks. Details are in the article and its Supporting Information.

## Licence

Code: MIT. Data and explorer: CC BY 4.0. 3Dmol.js (embedded in the explorer) is distributed under its own BSD licence.

## Contact

Gülbin Kurtay, Department of Chemistry, Hacettepe University, Ankara, Türkiye (gulbinkurtay@hacettepe.edu.tr; ORCID 0000-0003-0920-8409).
