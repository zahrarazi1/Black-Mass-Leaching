# Simulator: black-mass leaching campaign

The data in `data/raw/` are **synthetic**. They are produced by `simulate_leaching.py`, a simple physics-based model of a batch leach of Li-ion battery black mass (NMC-type cathode) in H₂SO₄ + H₂O₂, with measurement noise added. No real laboratory or plant data are used.

Regenerate the data with:

```bash
python simulator/simulate_leaching.py
```

The random seed is fixed, so you get exactly the same files every time.

## The process being simulated

Black mass is the powder left after batteries are discharged, shredded and sieved. It contains the cathode oxide (Li, Ni, Co, Mn), graphite, and pieces of the Al and Cu current-collector foils and the steel casing (Fe).

In a stirred 0.5 L reactor, black mass is leached with sulfuric acid. Hydrogen peroxide is added as a **reductant**, because:

- Co(III), Ni(III) and Mn(IV) in the cathode lattice dissolve poorly in acid alone. They must be reduced to the soluble M(II) state.
- Copper metal needs an **oxidant** to dissolve, and H₂O₂ also plays that role.

Example reaction:

2 LiCoO₂ + 3 H₂SO₄ + H₂O₂ → 2 CoSO₄ + Li₂SO₄ + 4 H₂O + O₂

## Model in five ideas

1. **Stoichiometry.** Every mole of metal dissolved consumes acid (Li 0.5, Ni/Co/Mn/Cu/Fe 1, Al 1.5 mol H₂SO₄ per mol). The reduced or oxidised fraction also consumes H₂O₂ (Co, Ni 0.5; Mn, Cu 1 mol per mol).
2. **Two fractions per metal.** One fraction dissolves with acid alone (`f_acid`). The rest dissolves only while H₂O₂ is present, at a rate scaled by c(H₂O₂) / (c(H₂O₂) + 0.02 mol/L).
3. **Shrinking-core kinetics, chemical-reaction control.** The conversion follows dX/dt = 3k(1 − X)^(2/3).
4. **Rate constant.** k = k_ref · exp(−Ea/R · (1/T − 1/333 K)) · (c_acid / 1.5)^0.8 · lot factor · run factor. So temperature (Arrhenius), the free acid that remains, the black-mass lot and run-to-run variability all matter.
5. **H₂O₂ is lost.** It is consumed by the reactions, and it also decomposes. Decomposition is first order, gets faster with temperature, and is catalysed by dissolved Co and Mn.

**Sensors:**
- **pH** follows the free acid that remains. As metals dissolve, acid is consumed and pH rises.
- **ORP** is high while H₂O₂ is present and drops once it is used up.

These two signals are the basis for a soft sensor.

The kinetic parameters (k_ref, Ea, f_acid) are rounded values chosen to give typical literature behaviour: >95 % Co/Ni extraction at 60–80 °C with enough H₂O₂, and poor Co/Ni and almost no Cu extraction without it. They are not fitted to one specific paper.

## Experimental design

Full factorial design, 3 × 3 × 3 × 2 = 54 conditions:

| Factor | Levels |
|---|---|
| Temperature | 40, 60, 80 °C |
| H₂SO₄ | 1.0, 1.5, 2.0 mol/L |
| H₂O₂ (30 wt-% solution) | 0, 2, 4 vol-% |
| Solid/liquid ratio | 50, 100 g/L |

- **Replicates and blocks:** each condition is run twice, once with black-mass lot A and once with lot B, so the **lot is a block**.
- **Run order:** randomised inside each block.
- **Size:** 108 runs × 12 samples (0, 5, 10, 15, 20, 30, 45, 60, 75, 90, 105, 120 min).
- **Calendar:** two runs per working day from March 2026.

## Files

### `data/raw/leaching_lab_log.csv`: one row per sample

| Column | Meaning | Unit |
|---|---|---|
| run_id | batch run identifier | – |
| sample_time | date and time the sample was taken | – |
| operator | person running the experiment | – |
| reactor | reactor used (R1 or R2) | – |
| lot | black-mass lot | – |
| temperature_set_C | temperature set-point | °C |
| h2so4_M | initial sulfuric acid concentration | mol/L |
| h2o2_vol_pct | H₂O₂ added (30 wt-% solution) | vol-% |
| solid_liquid_g_L | black mass per litre of leach solution | g/L |
| time_min | time since the black mass was added | min |
| temperature_meas | measured slurry temperature | °C |
| pH | online pH probe | – |
| orp_mV | online ORP probe, vs Ag/AgCl | mV |
| Li_mg_L … Fe_mg_L | metal concentration in the filtered sample, ICP-OES (detection limit 5 mg/L) | mg/L |

The file is written the way a real lab spreadsheet often looks, so it has to be checked and cleaned before use.

### `data/raw/black_mass_lots.csv`: assay of the two lots

Metal content of lots A and B in wt-%. You need it to turn concentrations into **extraction (%)**:

extraction_pct = c (mg/L) / (w (wt-%) / 100 × S/L (g/L) × 1000) × 100
