# config_utils.py
import re
import cobra
from cobra.flux_analysis import pfba

NGAM_FLUX = 0.398

# metabolite IDs to be filled according to each SBML model
CPD_O2    = ""  # e.g. oxygen metabolite id (base id / unique substring)
CPD_PI    = ""  # e.g. inorganic phosphate metabolite id
CPD_POLYP = ""  # e.g. polyphosphate metabolite id

ATP_MET_ID = ""  # intracellular ATP metabolite id used in the NGAM reaction

# reaction IDs (or unique substrings) for polyphosphate hydrolysis / synthesis
RXN_POLY_HYD = ""  # PolyP -> Pi
RXN_POLY_SYN = ""  # Pi -> PolyP

AA_MODE    = "auto"
AUXO_TABLE = "auxotrophy_summary.csv"
AA_UP_MAIN = 1.0
AA_UP_LOW  = 0.5

# optional amino acid exchange mapping; fill according to your model
AA_TO_SEED = {
    # "Ala": "your_ala_cpd_id",
    # "Val": "your_val_cpd_id",
    # "Leu": "your_leu_cpd_id",
}

# carbon sources to scan; keys are labels, values are metabolite ids or substrings
CARBON_CPDS = {
    # "Glucose": "your_glucose_cpd_id",
    # "Acetate": "your_acetate_cpd_id",
}


def sanitize(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]", "_", s)


def rxn_by_id(model, rid):
    if not rid:
        return None
    try:
        return model.reactions.get_by_id(rid)
    except KeyError:
        return None


def rxn_by_id_or_like(model, base_id):
    if not base_id:
        return None
    r = rxn_by_id(model, base_id)
    if r:
        return r
    candidates = [
        f"{base_id}_c0", f"{base_id}_c",
        f"R_{base_id}", f"R_{base_id}_c0", f"R_{base_id}_c",
        base_id.upper(), base_id.lower()
    ]
    for rid in candidates:
        r = rxn_by_id(model, rid)
        if r:
            return r
    low = base_id.lower()
    for r in model.reactions:
        if low and low in r.id.lower():
            return r
    return None


def met_like(model, seed, comps=("c0", "c")):
    seed = (seed or "").lower()
    if not seed:
        return None
    for m in model.metabolites:
        if seed in m.id.lower():
            comp = (m.compartment or "").lower()
            if comp.endswith(comps):
                return m
    return None


def find_exchange_for_cpd(model, seed_id):
    if not seed_id:
        return None
    key = seed_id.lower()
    exact, fuzzy = [], []
    for ex in model.exchanges:
        mets = list(ex.metabolites.keys())
        if len(mets) == 1 and key in mets[0].id.lower():
            exact.append(ex.id)
        if key in ex.id.lower():
            fuzzy.append(ex.id)
    for pool in (exact, fuzzy):
        if pool:
            pref = [r for r in pool if r.lower().endswith(("_e0", "_e"))]
            starts = [r for r in pref if r.startswith(("R_EX_", "EX_"))]
            if starts:
                return starts[0]
            if pref:
                return pref[0]
            starts2 = [r for r in pool if r.startswith(("R_EX_", "EX_"))]
            return starts2[0] if starts2 else pool[0]
    return None


def safe_box(r, lb, ub):
    r.upper_bound = 1000.0
    r.lower_bound = float(lb)
    r.upper_bound = float(ub)


def safe_fix(r, val):
    v = float(val)
    r.bounds = (v, v)


def solve_pfba(model):
    try:
        return pfba(model)
    except Exception:
        return model.optimize()


def find_o2_consumers(model):
    out = []
    pattern = (CPD_O2 or "").lower()
    if not pattern:
        return out
    for r in model.reactions:
        for met, coef in r.metabolites.items():
            if coef < 0 and pattern in met.id.lower():
                out.append(r.id)
                break
    return out


def uptake_sign(ex_rxn):
    if ex_rxn.lower_bound < 0 <= ex_rxn.upper_bound:
        return "neg"
    if ex_rxn.lower_bound >= 0 and ex_rxn.upper_bound > 0:
        return "pos"
    return "neg"


def force_uptake(ex_rxn, req_rate, forbid_secretion=True):
    s = uptake_sign(ex_rxn)
    if s == "neg":
        lb = -float(req_rate)
        ub = 0.0 if forbid_secretion else 1000.0
        safe_box(ex_rxn, lb, ub)
    else:
        lb = float(req_rate)
        safe_box(ex_rxn, lb, 1000.0)


def flux_to_uptake(ex_rxn, v):
    if uptake_sign(ex_rxn) == "neg":
        return max(0.0, -v)
    return max(0.0, v)


def ensure_ngam_reaction(model, ngam_flux):
    rid = "rxn00062_c0"
    r = model.reactions.get_by_id(rid)
    v = float(ngam_flux)
    try:
        atp_id = ATP_MET_ID or ""
        if atp_id:
            atp = model.metabolites.get_by_id(atp_id)
            coef = r.metabolites.get(atp, None)
            if coef is not None and coef > 0:
                v = -v
    except KeyError:
        pass
    safe_fix(r, v)


def add_internal_storage_reactions(model):
    stin_ids = {}
    stout_ids = {}
    for met in model.metabolites:
        comp = (met.compartment or "").lower()
        if not (comp.endswith("c") or comp.endswith("c0")):
            continue
        mid = met.id
        if "biomass" in mid.lower() or mid.lower().startswith("m_biomass"):
            continue
        stin_id = f"STIN_{sanitize(mid)}"
        stout_id = f"STOUT_{sanitize(mid)}"
        if stin_id not in model.reactions:
            rin = cobra.Reaction(stin_id)
            rin.lower_bound = 0.0
            rin.upper_bound = 0.0
            rin.add_metabolites({met: -1.0})
            model.add_reactions([rin])
        if stout_id not in model.reactions:
            rout = cobra.Reaction(stout_id)
            rout.lower_bound = 0.0
            rout.upper_bound = 0.0
            rout.add_metabolites({met: +1.0})
            model.add_reactions([rout])
        stin_ids[mid] = stin_id
        stout_ids[mid] = stout_id
    return stin_ids, stout_ids


def robust_biomass_internal(model):
    def is_boundary(rid):
        rid = rid.lower()
        return rid.startswith(("ex_", "r_ex_", "dm_", "sink_"))

    cands = []
    for r in model.reactions:
        rid = r.id.lower()
        nm = (r.name or "").lower()
        if is_boundary(rid):
            continue
        if ("biomass" in rid) or ("growth" in rid) or ("biomass" in nm) or ("growth" in nm) or (r.id == "bio1"):
            if len(r.metabolites) >= 2:
                cands.append(r.id)
    if not cands:
        raise RuntimeError("No biomass reaction found.")
    cands.sort(
        key=lambda x: (
            0
            if x == "bio1"
            else (0 if "biomass" in x.lower() else 1),
            -len(rxn_by_id(model, x).metabolites),
        )
    )
    return cands[0]


def open_amino_acids(model, mode, auxo_df=None, model_name=None):
    opened = []
    if mode == "none":
        return opened
    if mode == "all_low":
        for aa, seed in AA_TO_SEED.items():
            exid = find_exchange_for_cpd(model, seed)
            r = rxn_by_id(model, exid)
            if r:
                safe_box(r, -AA_UP_LOW, 1000.0)
                opened.append(r.id)
        return opened
    if auxo_df is None or model_name not in auxo_df.index:
        return opened
    row = auxo_df.loc[model_name]
    cols_norm = {c.strip().lower(): c for c in row.index}

    def marked(seed):
        key = seed.strip().lower()
        if key not in cols_norm:
            return False
        val = str(row[cols_norm[key]]).strip().lower()
        return val in {"1", "true", "yes", "y"}

    for aa, seed in AA_TO_SEED.items():
        if marked(seed):
            exid = find_exchange_for_cpd(model, seed)
            r = rxn_by_id(model, exid)
            if r:
                safe_box(r, -AA_UP_MAIN, 1000.0)
                opened.append(r.id)
    return opened


def carbon_number_from_ex(model, ex_id):
    ex = rxn_by_id(model, ex_id)
    if not ex:
        return None
    mets = list(ex.metabolites.keys())
    if len(mets) != 1:
        return None
    met = mets[0]
    formula = getattr(met, "formula", None)
    if not formula:
        return None
    m = re.search(r"C(\d+)", formula)
    if m:
        return int(m.group(1))
    if "C" in formula and not re.search(r"[A-Za-z]1", formula):
        return 1
    return None


def get_group_reactions(model, target_ids=("subsys_GLYCOLYSIS",), target_names=("GLYCOLYSIS",)):
    rxns = []
    for g in getattr(model, "groups", []):
        gid = (getattr(g, "id", "") or "").strip()
        gname = (getattr(g, "name", "") or "").strip().upper()
        if gid in target_ids or any(n.upper() in gname for n in target_names):
            for member in g.members:
                r = member if hasattr(member, "id") else rxn_by_id(model, str(member))
                if r is not None:
                    rxns.append(r)
    seen = set()
    uniq = []
    for r in rxns:
        if r.id not in seen:
            uniq.append(r)
            seen.add(r.id)
    return uniq


def find_like_reactions(model, keywords):
    keys = [k.lower() for k in keywords]
    out = []
    for r in model.reactions:
        txt = (r.id + " " + (r.name or "")).lower()
        if any(k in txt for k in keys):
            out.append(r)
    return out
