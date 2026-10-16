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
  [LIN]       Lineage 2023: 1.07 Mt CO2e over 3.86 B$ warehousing revenue (0.28 kg/$); about 109 kg CO2e
              and 393 $ revenue per pallet position per year
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
TL_KG_PER_MI = 0.48 / DIESEL_USD * DIESEL_KG_PER_GAL   # 1.34 kg/mile  [ATRI]
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
    facility[k] = (full_fixed, (full_fixed + min(rent_only, 0.28)) / 2, min(rent_only, 0.28))
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
# cost-side quantities from the baseline cost structure (numerics.Params): omega and the speed bands
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numerics as nm
from scipy.optimize import minimize_scalar
p0 = nm.Params()
def vopt(n, p, tau): return minimize_scalar(lambda v: nm.J(n, v, p, tau), bounds=(p.vmin, p.vmax), method="bounded",
                                            options={"xatol": 1e-12}).x
nK = minimize_scalar(lambda n: nm.J(n, vopt(n, p0, 0.0), p0, 0.0), bounds=(1, p0.N), method="bounded").x
vK = vopt(nK, p0, 0.0)
G = p0.z*p0.sigma*np.sqrt(p0.D/vK)*np.sqrt(p0.N)*(1-p0.rho)/(2*np.sqrt(nm.Psi(nK, p0)))
OMEGA = p0.h*G/(p0.F + p0.h*G)
RLOW, RHIGH = np.sqrt(nm.Psi(nK, p0)/p0.N), np.sqrt(nm.Psi(nK, p0)/nm.Psi(1.0, p0))
say(f"\nCOST SIDE (baseline cost structure): omega = {OMEGA:.3f}; speed bands {RLOW:.3f} / {RHIGH:.3f}")

def classify(rh, rt, rL, rF, omega=OMEGA, lo=RLOW, hi=RHIGH):
    thr = (1-omega)*rF + omega*(rh**0.8)*(max(rt, 1e-12)**0.2)
    disperse = rL > thr
    Rs = rh/max(rt, 1e-12)
    if Rs > hi or (Rs >= 1 and disperse): spd = "faster"
    elif Rs < lo or (Rs <= 1 and not disperse): spd = "slower"
    else: spd = "speed depends on size of density move"
    return spd, ("more dispersed" if disperse else "more consolidated"), Rs, thr

say("\nCLASSIFICATION AT CENTRAL VALUES (Proposition prop:classify)")
say(f"  {'product':<28}{'speed lever':<25}{'last mile':<22}{'R_speed':>8}{'thr':>7}{'r_L':>7}  result")
rows = []
for (pk, (val, st, emb)), (lk, lv), (mk, mv) in itertools.product(products.items(), speed_levers.items(), last_mile.items()):
    if pk == "Bulky low-value goods" and mk == "diesel van / courier": continue
    if pk != "Bulky low-value goods" and mk == "diesel local truck": continue
    rh = r_h(val[1], st); rt = lv[1]; rL = mv[1]; rF = facility[st][1]
    spd, den, Rs, thr = classify(rh, rt, rL, rF)
    rows.append((pk, lk, mk, Rs, thr, rL, spd, den))
    say(f"  {pk:<28}{lk:<25}{mk:<22}{Rs:8.3f}{thr:7.3f}{rL:7.3f}  {spd}, {den}")

# robustness of each cell over the full ranges and over omega in [0.1, 0.6]
say("\nROBUSTNESS: share of range corners (r_h, r_t, r_L, r_F, omega) giving the central-value classification")
for (pk, lk, mk, Rs, thr, rL, spd, den) in rows:
    val, st, emb = products[pk]
    rhs = r_h_range(val, st); rts = (speed_levers[lk][0], speed_levers[lk][2]); rLs = (last_mile[mk][0], last_mile[mk][2])
    rFs = (facility[st][0], facility[st][2])
    agree = tot = 0
    for rh_, rt_, rL_, rF_, om in itertools.product(rhs, rts, rLs, rFs, (0.1, OMEGA, 0.6)):
        s2, d2, *_ = classify(rh_, rt_, rL_, rF_, omega=om)
        tot += 1; agree += (s2 == spd and d2 == den)
    say(f"  {pk:<28}{lk:<25}{mk:<22}{100*agree/tot:5.0f}%")

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

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "calibration_v2_output.txt"), "w") as fh:
    fh.write("\n".join(out) + "\n")
