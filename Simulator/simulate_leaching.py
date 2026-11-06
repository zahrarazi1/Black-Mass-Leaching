"""
Simulated lab campaign: sulfuric-acid / hydrogen-peroxide leaching of
Li-ion battery black mass (NMC-type cathode material).

The data produced by this script are SYNTHETIC. They come from a simple
physics-based model (stoichiometry, shrinking-core kinetics, Arrhenius
temperature dependence, acid and peroxide consumption) plus measurement
noise. No real laboratory or plant data are used. See simulator/README.md.

Experimental design
-------------------
Full factorial design, 3 x 3 x 3 x 2 = 54 conditions:
    temperature      40 / 60 / 80 degC
    H2SO4            1.0 / 1.5 / 2.0 mol/L
    H2O2             0 / 2 / 4 vol-% (30 wt-% solution)
    solid/liquid     50 / 100 g/L
Each condition is run twice: replicate 1 with black-mass lot A, replicate 2
with lot B (lot = block). Run order is randomised inside each block.
Total 108 batch runs, 12 samples each (0-120 min).

Output files (data/raw/)
------------------------
    leaching_lab_log.csv   one row per sample, as typed in the lab (messy)
    black_mass_lots.csv    chemical assay of the two black-mass lots

Run:  python simulator/simulate_leaching.py
"""

from pathlib import Path

import numpy as np
import pandas as pd

SEED = 2026
R = 8.314                          # J/(mol K)
T_REF = 333.15                     # 60 degC, reference temperature for k
VOLUME_L = 0.5                     # reactor volume
SAMPLE_TIMES = [0, 5, 10, 15, 20, 30, 45, 60, 75, 90, 105, 120]   # min
DT = 0.25                          # integration step, min
H2O2_MOL_PER_VOLPCT = 0.098        # 1 vol-% of 30 wt-% H2O2 = 0.098 mol/L

# ------------------------------------------------------------ chemistry data
METALS = ["Li", "Co", "Ni", "Mn", "Al", "Cu", "Fe"]
MOLAR_MASS = {"Li": 6.94, "Co": 58.93, "Ni": 58.69, "Mn": 54.94,
              "Al": 26.98, "Cu": 63.55, "Fe": 55.85}

# black-mass assay, wt-%
LOTS = {
    "A": {"Li": 3.6, "Co": 6.4, "Ni": 19.5, "Mn": 6.1, "Al": 2.2, "Cu": 2.6, "Fe": 0.6},
    "B": {"Li": 3.3, "Co": 5.9, "Ni": 18.2, "Mn": 5.6, "Al": 3.4, "Cu": 3.1, "Fe": 0.9},
}
LOT_RATE_FACTOR = {"A": 1.0, "B": 0.85}     # lot B is coarser -> slower

# kinetic parameters (shrinking core, chemical reaction control)
#   k_ref  rate constant at 60 degC and 1.5 mol/L H2SO4, 1/min
#   Ea     apparent activation energy, J/mol
#   f_acid fraction that dissolves with acid alone; the rest needs H2O2
#          (Co(III), Ni(III), Mn(IV) must be reduced to M(II); Cu metal must be oxidised)
KIN = {
    "Li": {"k_ref": 0.020, "Ea": 30e3, "f_acid": 0.80},
    "Co": {"k_ref": 0.010, "Ea": 55e3, "f_acid": 0.40},
    "Ni": {"k_ref": 0.009, "Ea": 55e3, "f_acid": 0.55},
    "Mn": {"k_ref": 0.010, "Ea": 50e3, "f_acid": 0.30},
    "Al": {"k_ref": 0.012, "Ea": 45e3, "f_acid": 1.00},
    "Cu": {"k_ref": 0.015, "Ea": 35e3, "f_acid": 0.03},
    "Fe": {"k_ref": 0.004, "Ea": 40e3, "f_acid": 1.00},
}
# mol H2SO4 consumed per mol metal dissolved
ACID_STOICH = {"Li": 0.5, "Co": 1.0, "Ni": 1.0, "Mn": 1.0, "Al": 1.5, "Cu": 1.0, "Fe": 1.0}
# mol H2O2 consumed per mol metal dissolved in the reductive / oxidative fraction
H2O2_STOICH = {"Li": 0.0, "Co": 0.5, "Ni": 0.5, "Mn": 1.0, "Al": 0.0, "Cu": 1.0, "Fe": 0.0}

K_H2O2 = 0.02          # mol/L, half-saturation of the peroxide effect
K_DEC_REF = 0.002      # 1/min, H2O2 thermal decomposition at 60 degC
EA_DEC = 60e3          # J/mol
ACID_ORDER = 0.8
ICP_LOD_MG_L = 5.0


def arrhenius(ea, temp_k):
    """Factor by which a rate changes from 60 degC to temp_k."""
    return np.exp(-ea / R * (1 / temp_k - 1 / T_REF))


def simulate_run(cond, rng):
    """Integrate one batch leach. Returns list of true-state dicts at sample times."""
    temp_k = cond["temperature_C"] + 273.15
    sl = cond["solid_liquid_g_L"]
    lot = LOTS[cond["lot"]]
    run_factor = cond["run_factor"]

    # total moles of each metal per litre of slurry
    n_tot = {m: sl * lot[m] / 100 / MOLAR_MASS[m] for m in METALS}
    x_acid = {m: 0.0 for m in METALS}     # conversion of the acid-soluble fraction
    x_red = {m: 0.0 for m in METALS}      # conversion of the fraction needing H2O2
    c_acid = cond["h2so4_M"]
    c_h2o2 = cond["h2o2_vol_pct"] * H2O2_MOL_PER_VOLPCT

    out = []
    t = 0.0
    sample_idx = 0
    n_steps = int(SAMPLE_TIMES[-1] / DT)
    for step in range(n_steps + 1):
        t = step * DT
        if sample_idx < len(SAMPLE_TIMES) and abs(t - SAMPLE_TIMES[sample_idx]) < 1e-9:
            extraction = {m: KIN[m]["f_acid"] * x_acid[m] + (1 - KIN[m]["f_acid"]) * x_red[m]
                          for m in METALS}
            out.append({"time_min": t, "c_acid": c_acid, "c_h2o2": c_h2o2,
                        "extraction": extraction, "n_tot": n_tot})
            sample_idx += 1

        acid_term = (max(c_acid, 0.0) / 1.5) ** ACID_ORDER
        h2o2_term = c_h2o2 / (c_h2o2 + K_H2O2)
        d_acid = 0.0
        d_h2o2 = 0.0
        dissolved_tm = 0.0
        for m in METALS:
            k = (KIN[m]["k_ref"] * arrhenius(KIN[m]["Ea"], temp_k) * acid_term
                 * LOT_RATE_FACTOR[cond["lot"]] * run_factor)
            # shrinking core, chemical control: dX/dt = 3 k (1 - X)^(2/3)
            dxa = 3 * k * (1 - x_acid[m]) ** (2 / 3) * DT
            dxr = 3 * k * h2o2_term * (1 - x_red[m]) ** (2 / 3) * DT
            dxa = min(dxa, 1 - x_acid[m])
            dxr = min(dxr, 1 - x_red[m])
            x_acid[m] += dxa
            x_red[m] += dxr
            dn_acid = n_tot[m] * KIN[m]["f_acid"] * dxa
            dn_red = n_tot[m] * (1 - KIN[m]["f_acid"]) * dxr
            d_acid += (dn_acid + dn_red) * ACID_STOICH[m]
            d_h2o2 += dn_red * H2O2_STOICH[m]
            if m in ("Co", "Mn"):
                dissolved_tm += n_tot[m] * (KIN[m]["f_acid"] * x_acid[m]
                                            + (1 - KIN[m]["f_acid"]) * x_red[m])

        # peroxide also decomposes: faster when hot and when Co/Mn ions catalyse it
        k_dec = K_DEC_REF * arrhenius(EA_DEC, temp_k) * (1 + 5 * dissolved_tm)
        d_h2o2 += k_dec * c_h2o2 * DT

        c_acid = max(c_acid - d_acid, 0.0)
        c_h2o2 = max(c_h2o2 - d_h2o2, 0.0)
    return out


def design_campaign(rng):
    """Full factorial x 2 blocks (lots), randomised run order, lab calendar."""
    levels = pd.MultiIndex.from_product(
        [[40, 60, 80], [1.0, 1.5, 2.0], [0, 2, 4], [50, 100]],
        names=["temperature_C", "h2so4_M", "h2o2_vol_pct", "solid_liquid_g_L"],
    ).to_frame(index=False)

    blocks = []
    for lot in ["A", "B"]:
        b = levels.copy()
        b["lot"] = lot
        blocks.append(b.sample(frac=1, random_state=int(rng.integers(1e6))))   # randomise order
    runs = pd.concat(blocks, ignore_index=True)

    n = len(runs)
    runs["run_id"] = [f"L{i:03d}" for i in range(1, n + 1)]
    runs["run_factor"] = np.exp(rng.normal(0, 0.07, n))        # run-to-run variability
    runs["reactor"] = rng.choice(["R1", "R2"], size=n)
    runs["operator"] = rng.choice(["AK", "JS", "MK"], size=n, p=[0.4, 0.3, 0.3])
    runs["ph_offset"] = rng.normal(0, 0.04, n)                   # pH probe calibration
    runs["orp_offset"] = rng.normal(0, 8, n)                     # mV

    # two runs per working day from 2 March 2026, starting 08:30 and 13:00
    days = pd.bdate_range("2026-03-02", periods=n // 2 + 1)
    starts = []
    for i in range(n):
        day = days[i // 2]
        base = pd.Timedelta(hours=8, minutes=30) if i % 2 == 0 else pd.Timedelta(hours=13)
        starts.append(day + base + pd.Timedelta(minutes=int(rng.integers(0, 25))))
    runs["start"] = starts
    return runs


def measure(runs, rng):
    """Convert true states to what the instruments would report."""
    rows = []
    for _, c in runs.iterrows():
        states = simulate_run(c, rng)
        for s in states:
            t = s["time_min"]
            row = {
                "run_id": c["run_id"],
                "sample_time": c["start"] + pd.Timedelta(minutes=t),
                "operator": c["operator"],
                "reactor": c["reactor"],
                "lot": c["lot"],
                "temperature_set_C": c["temperature_C"],
                "h2so4_M": c["h2so4_M"],
                "h2o2_vol_pct": c["h2o2_vol_pct"],
                "solid_liquid_g_L": c["solid_liquid_g_L"],
                "time_min": t,
            }
            # small exotherm in the first 15 min, stronger at high pulp density
            exo = (c["solid_liquid_g_L"] / 100) * 2.5 * np.exp(-t / 12) if t > 0 else 0
            row["temperature_meas"] = c["temperature_C"] + exo + rng.normal(0, 0.4)

            h_plus = 1.1 * s["c_acid"] + 10 ** -2.8
            row["pH"] = -np.log10(h_plus) + c["ph_offset"] + rng.normal(0, 0.02)

            g = s["c_h2o2"] / (s["c_h2o2"] + K_H2O2)
            row["orp_mV"] = (420 + 180 * g + 25 * np.log10(s["c_acid"] + 0.05)
                             + c["orp_offset"] + rng.normal(0, 4))

            for m in METALS:
                conc = s["extraction"][m] * s["n_tot"][m] * MOLAR_MASS[m] * 1000    # mg/L
                conc = conc * (1 + rng.normal(0, 0.03)) + rng.normal(0, 1.0)
                row[f"{m}_mg_L"] = max(conc, 0.0)
            rows.append(row)
    return pd.DataFrame(rows)


def make_messy(df, runs, rng):
    """Write the log the way real lab spreadsheets often look."""
    df = df.copy()
    n = len(df)

    # numbers as a lab would round them
    df["temperature_meas"] = df["temperature_meas"].round(1)
    df["pH"] = df["pH"].round(2)
    df["orp_mV"] = df["orp_mV"].round(0)
    for m in METALS:
        df[f"{m}_mg_L"] = df[f"{m}_mg_L"].round(1)

    # 1. operator MK writes dates as dd/mm/yyyy HH:MM, everyone else ISO
    iso = df["sample_time"].dt.strftime("%Y-%m-%d %H:%M")
    eu = df["sample_time"].dt.strftime("%d/%m/%Y %H:%M")
    df["sample_time"] = np.where(df["operator"] == "MK", eu, iso)

    # 2. reactor R2's logger exports Kelvin after a firmware update on 1 May
    k_mask = (df["reactor"] == "R2") & (pd.to_datetime(iso) >= "2026-05-01")
    df.loc[k_mask, "temperature_meas"] = (df.loc[k_mask, "temperature_meas"] + 273.15).round(1)

    # 3. operator JS reads ORP against SHE instead of Ag/AgCl (+197 mV)
    js = df["operator"] == "JS"
    df.loc[js, "orp_mV"] = df.loc[js, "orp_mV"] + 197

    # 4. some runs have the acid typed in g/L instead of mol/L
    gl_runs = rng.choice(runs["run_id"], 9, replace=False)
    gl = df["run_id"].isin(gl_runs)
    df.loc[gl, "h2so4_M"] = (df.loc[gl, "h2so4_M"] * 98.08).round(0)

    # 5. lot code typed by hand
    spell = {"A": ["A", "A", "lot A", "LotA", "a"], "B": ["B", "B", "lot B", "LotB", "b"]}
    df["lot"] = [spell[v][rng.integers(0, 5)] if rng.uniform() < 0.25 else v for v in df["lot"]]

    # 6. pH probe dropouts logged as -99
    df.loc[rng.uniform(size=n) < 0.02, "pH"] = -99

    # 7. ICP: broken vials (whole sample missing) and values below LOD as "<LOD"
    broken = rng.uniform(size=n) < 0.03
    for m in METALS:
        col = f"{m}_mg_L"
        values = df[col].astype(object)
        values[df[col] < ICP_LOD_MG_L] = "<LOD"
        values[broken] = np.nan
        df[col] = values

    # 8. three typing slips in Co (decimal point lost)
    numeric_co = pd.to_numeric(df["Co_mg_L"], errors="coerce")
    idx = rng.choice(df.index[numeric_co > 500], 3, replace=False)
    df.loc[idx, "Co_mg_L"] = (numeric_co[idx] * 10).round(1)

    # 9. rows pasted twice when two spreadsheets were merged
    dup = df.sample(n=24, random_state=SEED)
    df = pd.concat([df, dup]).sort_index(kind="stable").reset_index(drop=True)

    cols = ["run_id", "sample_time", "operator", "reactor", "lot", "temperature_set_C",
            "h2so4_M", "h2o2_vol_pct", "solid_liquid_g_L", "time_min", "temperature_meas",
            "pH", "orp_mV"] + [f"{m}_mg_L" for m in METALS]
    return df[cols]


def main():
    rng = np.random.default_rng(SEED)
    runs = design_campaign(rng)
    clean = measure(runs, rng)
    messy = make_messy(clean, runs, rng)

    out = Path(__file__).resolve().parents[1] / "data" / "raw"
    out.mkdir(parents=True, exist_ok=True)
    messy.to_csv(out / "leaching_lab_log.csv", index=False)

    lots = pd.DataFrame(LOTS).T.rename_axis("lot").reset_index()
    lots.columns = ["lot"] + [f"{m}_wt_pct" for m in METALS]
    lots.to_csv(out / "black_mass_lots.csv", index=False)

    print(f"{len(runs)} runs, {len(messy)} rows -> {out / 'leaching_lab_log.csv'}")
    return runs, clean


if __name__ == "__main__":
    main()
