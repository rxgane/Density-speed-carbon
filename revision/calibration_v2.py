"""
v2 calibration: carbon intensity per dollar for each activity, derived from public sources,
and the classification of realistic products by the joint-design result (Proposition prop:classify).

Rule used throughout:  r = (carbon of the activity) / (cost of the activity), in kg CO2e per US dollar,
computed bottom-up as (energy use x emission factor) / (activity cost) and cross-checked against the
EPA supply-chain (USEEIO v1.3) industry factors, kg CO2e per 2022 USD.

Sources (all opened):
  [EPA-HUB]   EPA GHG Emission Factors Hub: diesel 10.21 kg CO2/gal; US grid 0.371 kg CO2e/kWh
              https://www.epa.gov/system/files/documents/2022-04/ghg_emission_factors_hub.pdf
  [EIA-DSL]   EIA weekly retail diesel, 2025 average 3.66 $/gal
              https://www.eia.gov/petroleum/supply/weekly/pdf/table14.pdf
  [EIA-ELEC]  EIA average electricity price, 8.1 (industrial) to 12.8 (commercial) cents/kWh
              https://www.eia.gov/electricity/annual/html/epa_02_10.html
  [MRV]       EU MRV: heavy fuel oil 3.114 t CO2 per t fuel  https://eur-lex.europa.eu/eli/reg_del/2016/2071/oj/eng
  [BUNKER]    VLSFO average 535-626 $/t (2025)
              https://shipandbunker.com/news/world/223373-analysis-average-vlsfo-bunker-prices-at-key-ports-fall-20-in-2025
  [EEIO]      USEEIO v1.3 supply-chain factors (kg CO2e/2022 USD, with margins): trucking 484121 0.595;
              local trucking 484110 0.595; couriers 492110 0.303; rail 482111 0.466; general warehousing
              493110 0.244; refrigerated warehousing 493120 0.244; nonresidential lessors 531120 0.246;
              animal slaughtering 311611 1.12; poultry processing 311615 0.391; fluid milk 311511 0.803;
              cheese 311513 0.969; vegetable farming 111219 0.631
              https://pasteur.epa.gov/uploads/10.23719/1531143/SupplyChainGHGEmissionFactors_v1.3.0_NAICS_CO2e_USD2022.csv
  [ATRI]      ATRI operational costs: truck 2.26 $/mile, of which fuel 0.48 $/mile
              https://www.fleetmaintenance.com/equipment/article/55301363/american-transportation-research-institute-atri-breakdown-of-atri-2025-operational-costs-report
  [IMX]       Intermodal about 60% less carbon per load than truckload; intermodal 10-30% cheaper
              https://talkbusiness.net/2023/10/intermodal-shippers-see-emissions-cost-benefits-survey-shows
              https://www.freightwaves.com/news/massive-opportunities-for-j-b-hunt-in-intermodal-shift
  [AIR]       Air about 4-5x road cost (GCCA); air 1.10 kg/t-km incl. radiative forcing (Climatiq/DEFRA)
              https://www.gcca.org/sites/default/files/Trends%20in%20the%20Cold%20Chain%20%E2%80%93%20From%20Air%20to%20Sea.pdf
  [TEAM]      Team drivers cover 4,500-5,000 vs 2,400-2,800 miles/week; team pay premium about
              0.25-0.40 $/mile; same tractor and fuel
              https://www.freightwaves.com/checkpoint/team-truck-driving/
  [ESTAR]     ENERGY STAR US medians: refrigerated warehouse 84.1, non-refrigerated 22.7 kBtu/ft2 (site)
              https://portfoliomanager.energystar.gov/pdf/reference/US%20National%20Median%20Table.pdf
  [CW]        Cushman & Wakefield Q2 2025: industrial asking rent 10.12 $/ft2/yr
              https://www.cushmanwakefield.com/en/united-states/news/2025/07/q2-2025-us-industrial-market-report
  [LIN]       Lineage 2024 Scope 1 + market-based Scope 2: 1.068 Mt CO2e (2024 Sustainability Report,
              published 2025, pp. 15, 39). Dividing by the March 2024 capacity of 9.8M pallet positions
              (prospectus) gives about 109 kg CO2e per pallet position-year. This company-wide proxy
              includes mobile combustion and refrigerants. Ambient 30 kg is 109*(22.7/84.1), rounded.
              The 0.28 kg/$ cross-check divides this total by about $3.86B in 2023 warehousing revenue;
              its periods and operating boundaries differ, so it is only an approximate comparison.
              https://www.onelineage.com/sites/default/files/2025-06/Lineage-2024-Sustainability-Report_6_25_25.pdf
              https://www.sec.gov/Archives/edgar/data/1868159/000119312524185331/d577649d424b4.htm
  [PAL]       Ambient pallet storage about 20 $/pallet-month (range 12-25); refrigerated 1.5-3x ambient
              https://olimpwarehousing.com/pallet-storage-cost-per-month-usa/
              https://racklify.com/encyclopedia/temperature-controlled-storage-fee-pricing-how-rates-are-calculated/
  [CARRY]     Inventory carrying cost 15-25% of value per year in total; capital is one component
              https://www.cleverence.com/articles/for-business/carrying-cost-inventory-5982/
  [EVAN]      Electric van 0.37-0.39 kWh/mile (AFDC)  https://afdc.energy.gov/files/u/publication/cost-ev-vans-pickups-us-2040-jan22.pdf
"""
import os, sys, itertools
import numpy as np

out = []
def say(s=""): print(s); out.append(s)

# ------------------------------------------------------------------ energy primitives
DIESEL_KG_PER_GAL, DIESEL_USD = 10.21, 3.66          # [EPA-HUB], [EIA-DSL]
GRID_KG_PER_KWH = 0.371                              # [EPA-HUB]
ELEC_USD = (0.081, 0.128)                            # [EIA-ELEC]
kg_per_diesel_usd = DIESEL_KG_PER_GAL / DIESEL_USD   # 2.79
kg_per_elec_usd = tuple(GRID_KG_PER_KWH / p for p in ELEC_USD)

say("ENERGY PRIMITIVES (kg CO2e per dollar of energy)")
say(f"  diesel {kg_per_diesel_usd:.2f}   grid electricity {kg_per_elec_usd[1]:.2f}-{kg_per_elec_usd[0]:.2f}"
    f"   heavy fuel oil {3114/626:.1f}-{3114/535:.1f}")

# ------------------------------------------------------------------ inbound speed levers: r_t = marginal kg / marginal $
TL_KG_PER_MI = 1.52   # kg/mile, EPA SmartWay 2023 truckload dry-van median (fuel cost basis gives 1.34)
TL_USD_PER_MI = 2.26
speed_levers = {
    # name: (low, central, high, how)
    "ocean vessel speed": (3114/626, 3114/580, 3114/535, "marginal cost of speed is bunker fuel [MRV][BUNKER]"),
    "intermodal -> truckload": (0.60*TL_KG_PER_MI/(0.30*TL_USD_PER_MI), 0.60*TL_KG_PER_MI/(0.20*TL_USD_PER_MI),
                                0.60*TL_KG_PER_MI/(0.10*TL_USD_PER_MI), "60% carbon saving, 10-30% cost saving [ATRI][IMX]"),
    "truckload -> air": (0.6, 1.1, 1.9, "air 4-5x road cost, about 5-10x road carbon per t-km [AIR]"),
    "solo -> team drivers": (0.0, 0.03, 0.10, "same tractor and fuel; premium is driver labor [TEAM]"),
}
say("\nINBOUND SPEED LEVERS: r_t = marginal carbon / marginal cost of faster inbound (kg/$)")
for k, (lo, c, hi, how) in speed_levers.items():
    say(f"  {k:<26} {c:5.2f}  (range {lo:.2f}-{hi:.2f})   {how}")

# ------------------------------------------------------------------ last mile: r_L
last_mile = {
    "diesel van / courier": (0.15, 0.30, 0.45, "EEIO couriers 0.303; fuel share x 2.79 kg/$"),
    "diesel local truck":   (0.30, 0.50, 0.60, "EEIO local trucking 0.595; ATRI fuel share 21% -> 0.59"),
    "electric van":         (0.05, 0.075, 0.10, "0.37-0.39 kWh/mile x 0.371 kg/kWh over van cost per mile [EVAN]"),
    "electric truck":       (0.12, 0.16, 0.20, "about 1.2 kWh/mile x 0.371 kg/kWh over about 2.75 $/mile; driver-dominated cost"),
}
say("\nLAST MILE: r_L (kg/$)")
for k, (lo, c, hi, how) in last_mile.items():
    say(f"  {k:<26} {c:5.3f}  (range {lo:.2f}-{hi:.2f})   {how}")

# ------------------------------------------------------------------ facilities: r_F
kwh_per_ft2 = {"ambient": 22.7/3.412, "refrigerated": 84.1/3.412}   # [ESTAR], site energy, all-electric bound
facility = {}
for k, kwh in kwh_per_ft2.items():
    kg_ft2 = kwh * GRID_KG_PER_KWH
    rent = 10.12 if k == "ambient" else 2.0 * 10.12                 # [CW]; refrigerated rent premium assumed 2x
    rent_only = kg_ft2 / rent
    full_fixed = kg_ft2 / (rent / 0.4)                              # rent assumed 40% of a DC's fixed operating cost
    # the model's F is the full fixed operating cost of a DC, so the central value uses full fixed cost;
    # the range runs from 70% of that value up to the rent-only ratio (capped at Lineage's 0.28)
    facility[k] = (0.7 * full_fixed, full_fixed, min(rent_only, 0.28))
say("\nFACILITIES: r_F (kg/$)  [ESTAR site energy x grid factor, over fixed cost]")
for k, (lo, c, hi) in facility.items():
    say(f"  {k:<26} {c:5.3f}  (range {lo:.2f}-{hi:.2f})   cross-check: EEIO warehousing 0.244, Lineage 0.28")

# ------------------------------------------------------------------ holding: r_h per product (per pallet-year)
CAP = (0.08, 0.10, 0.15)                       # capital charge on inventory value, per year [CARRY]
storage_kg = {"ambient": (20, 30, 50), "refrigerated": (80, 109, 140)}   # kg CO2e / pallet-yr [LIN], [ESTAR] 3.7x
storage_usd = {"ambient": (144, 240, 300), "refrigerated": (360, 480, 720)}  # $/pallet-yr [PAL]
products = {
    # name: (pallet value $ (lo,c,hi), storage type, embodied kg/$ of product [EEIO] or None)
    "Ambient packaged goods":   ((1000, 1500, 2500), "ambient", None),
    "Bulky low-value goods":    ((300, 500, 800), "ambient", None),
    "Refrigerated dairy":       ((1000, 1250, 1500), "refrigerated", 0.803),
    "Fresh beef":               ((4000, 6000, 8000), "refrigerated", 1.12),
    "Fresh poultry":            ((2000, 3000, 4000), "refrigerated", 0.391),
    "Cold-chain pharmaceuticals": ((15000, 22050, 30000), "refrigerated", None),
}
def r_h(value, stype, emb=None, s_week=0.0, which=1):
    """carbon per dollar of holding; s_week = fraction of stock value lost to spoilage per week held"""
    E = storage_kg[stype][which]; H = CAP[which]*value + storage_usd[stype][which]
    spoil_usd = 52*s_week*value
    return (E + spoil_usd*(emb or 0.0)) / (H + spoil_usd)
def r_h_range(value, stype, emb=None, s_week=0.0):
    vals = [ (storage_kg[stype][a] + 52*s_week*v*(emb or 0)) / (CAP[b]*v + storage_usd[stype][c] + 52*s_week*v)
             for a, b, c in itertools.product(range(3), range(3), range(3)) for v in value ]
    return min(vals), max(vals)

say("\nHOLDING: r_h (kg/$), storage energy only (no spoilage)")
for k, (val, st, emb) in products.items():
    lo, hi = r_h_range(val, st)
    say(f"  {k:<26} {r_h(val[1], st):5.3f}  (range {lo:.3f}-{hi:.3f})")

# ------------------------------------------------------------------ classification by Proposition prop:classify
# cost side (pooling weight, speed bands) from the calibrated network of application_v2.py
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import contextlib, io
with contextlib.redirect_stdout(io.StringIO()):
    import application_v2 as app
_, _, OMEGA, RLOW, RHIGH = app.cost_side(dict(app.base, **app.lane("fuel")))
OMEGAS = (OMEGA, 0.1, 0.3, 0.6)
say(f"\nCOST SIDE (calibrated network): omega = {OMEGA:.3f}; speed bands {RLOW:.3f} / {RHIGH:.3f}; "
    f"robustness also over omega in {OMEGAS[1:]}")

def classify(rh, rt, rL, rF, omega=OMEGA, lo=RLOW, hi=RHIGH):
    thr = (1-omega)*rF + omega*(rh**0.8)*(max(rt, 1e-12)**0.2)
    disperse = rL > thr
    Rs = rh/max(rt, 1e-12)
    if Rs > hi or (Rs >= 1 and disperse): spd = "faster"
    elif Rs < lo or (Rs <= 1 and not disperse): spd = "slower"
    else: spd = "depends"
    return spd, ("more DCs" if disperse else "fewer DCs"), Rs, thr

say("\nTABLE: classification; each cell = direction at most range corners (share of corners agreeing)")
say(f"  {'product':<28}{'fuel speed':>14}{'labor speed':>14}{'diesel truck':>15}{'diesel van':>13}{'elec. truck':>13}{'elec. van':>12}")
for pk, (val, st, emb) in products.items():
    rhs = r_h_range(val, st); rFs = (facility[st][0], facility[st][2]); cells = []
    for lk in ("intermodal -> truckload", "solo -> team drivers"):
        rts = (speed_levers[lk][0], speed_levers[lk][2]); cnt = {}; t = 0
        for rh_, rt_, rL_, rF_, o in itertools.product(rhs, rts, (0.15, 0.6), rFs, OMEGAS):
            s2, *_ = classify(rh_, rt_, rL_, rF_, omega=o); cnt[s2] = cnt.get(s2, 0)+1; t += 1
        best = max(cnt, key=cnt.get); cells.append(f"{best} {round(100*cnt[best]/t)}%")
    for mk in ("diesel local truck", "diesel van / courier", "electric truck", "electric van"):
        rLs = (last_mile[mk][0], last_mile[mk][2]); cnt = {}; t = 0
        for rh_, rt_, rL_, rF_, o in itertools.product(rhs, (0.03, 1.1, 1.78, 5.4), rLs, rFs, OMEGAS):
            _, d2, *_ = classify(rh_, rt_, rL_, rF_, omega=o); cnt[d2] = cnt.get(d2, 0)+1; t += 1
        best = max(cnt, key=cnt.get); cells.append(f"{best.split()[0]} {round(100*cnt[best]/t)}%")
    say(f"  {pk:<28}" + "".join(f"{x:>14}" for x in cells))

# ------------------------------------------------------------------ spoilage break-even
say("\nSPOILAGE BREAK-EVEN: weekly spoilage (share of stock value lost per week held) at which faster")
say("replenishment becomes carbon-minimizing whatever the network does (R_speed > upper band)")
say("Holding carbon then includes the embodied carbon of spoiled food [EEIO]; holding cost includes its value.")
for pk in ("Refrigerated dairy", "Fresh beef", "Fresh poultry"):
    val, st, emb = products[pk]
    V = val[1]; E = storage_kg[st][1]; H0 = CAP[1]*V + storage_usd[st][1]
    for lk in ("intermodal -> truckload", "truckload -> air", "ocean vessel speed"):
        k = speed_levers[lk][1]*RHIGH
        if emb <= k:
            say(f"  {pk:<22} vs {lk:<24}: never (embodied {emb:.2f} kg/$ <= {k:.2f} needed)")
        else:
            s = (k*H0 - E)/(52*V*(emb - k))
            say(f"  {pk:<22} vs {lk:<24}: above {100*max(s,0):.2f}% of stock value per week")

# ------------------------------------------------------------------ carbon share of priced cost (Proposition prop:path)
say("\nCARBON SHARE OF PRICED COST, s = tau r / (1 + tau r), at central ratios; and the price at which s = 1/2")
share_rows = [("holding, ambient packaged goods", r_h(products["Ambient packaged goods"][0][1], "ambient")),
              ("facilities, ambient DC", facility["ambient"][1]),
              ("last mile, diesel local truck", last_mile["diesel local truck"][1]),
              ("last mile, diesel van", last_mile["diesel van / courier"][1]),
              ("line haul, truck instead of intermodal", speed_levers["intermodal -> truckload"][1]),
              ("line haul, faster vessel speed", speed_levers["ocean vessel speed"][1])]
say(f"  {'activity':<40}{'$50/t':>8}{'$100/t':>8}{'$200/t':>8}{'s=1/2 at $/t':>14}")
for nm_, r in share_rows:
    sh = [100*t*r/(1+t*r) for t in (0.05, 0.10, 0.20)]
    say(f"  {nm_:<40}{sh[0]:7.1f}%{sh[1]:7.1f}%{sh[2]:7.1f}%{1000/r:14,.0f}")

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "calibration_v2_output.txt"), "w") as fh:
    fh.write("\n".join(out) + "\n")
