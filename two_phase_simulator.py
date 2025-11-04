# two_phase_simulator.py
from config_utils import (
    NGAM_FLUX,
    CPD_O2,
    CPD_PI,
    CPD_POLYP,
    RXN_POLY_HYD,
    RXN_POLY_SYN,
    sanitize,
    rxn_by_id,
    rxn_by_id_or_like,
    met_like,
    find_exchange_for_cpd,
    safe_box,
    safe_fix,
    solve_pfba,
    find_o2_consumers,
    force_uptake,
    flux_to_uptake,
    ensure_ngam_reaction,
    add_internal_storage_reactions,
    get_group_reactions,
    find_like_reactions,
)

# ---- time / dose configuration (typical ranges; tune per system) ----
ANA_TIME    = 1.0   # ~1–2 h
AER_TIME    = 3.0   # ~2–4 h
DT          = 0.2   # ~0.1–0.5 h
FEED_CUTOFF = 1.0   # within anaerobic window
DOSE_MOLAR  = 2.0   # ~1–3 mmol per gDW

SUB_MAX_ANA  =   # large enough upper bound
REUP_MAX_AER =
INT_REL_MAX  =

O2_SUPPLY    =   # typical O2 uptake cap, e.g. 5–20 (model units)
INIT_BIOMASS = 1.0   # reference biomass unit

INIT_PI_POOL =   # effective external Pi pool (arbitrary model units)

# PolyP-related pools and flux bounds (order-of-magnitude choices)
INIT_POLYP   =    # initial internal PolyP pool
POLY_REQ_LB  =    # lower bound for hydrolysis in anaerobic phase
POLY_SYN_EPS =   # minimal synthesis flux in aerobic phase

SOFT_SUPPRESS_GLY_ATP = True
SOFT_EPS              = 0.0
GLY_ATP_CAND_IDS      = ["PGK", "PYK"]  # matched by substring in reaction id/name


def integrate_two_phase(
    model,
    ex_sub_id,
    ex_o2_id,
    bio_id,
    ana_time,
    aer_time,
    dt,
    dose,
    feed_cut,
    sub_max_ana,
    reup_max_aer,
    int_rel_max,
    o2_supply,
    seed_pi_ex=None,
    init_pi_amount=0.0,
    sub_max_ANA=None,
    reup_max_AER=None,
    int_rel_MAX=None,
):
    if sub_max_ANA is not None:
        sub_max_ana = sub_max_ANA
    if reup_max_AER is not None:
        reup_max_aer = reup_max_AER
    if int_rel_MAX is not None:
        int_rel_max = int_rel_MAX

    tA = int(round(ana_time / dt))
    tO = int(round(aer_time / dt))
    tW = int(round(feed_cut / dt))

    biomass = rxn_by_id(model, bio_id)
    ex_sub  = rxn_by_id(model, ex_sub_id) if ex_sub_id else None
    ex_o2   = rxn_by_id(model, ex_o2_id)  if ex_o2_id  else None
    ex_pi   = rxn_by_id(model, seed_pi_ex) if seed_pi_ex else None

    r_hyd = rxn_by_id_or_like(model, RXN_POLY_HYD)
    r_syn = rxn_by_id_or_like(model, RXN_POLY_SYN)
    met_polyp = met_like(model, CPD_POLYP)
    met_pi    = met_like(model, CPD_PI)

    ensure_ngam_reaction(model, NGAM_FLUX)

    stin_map, stout_map = add_internal_storage_reactions(model)
    ext_pool = {}
    int_pool = {m: 0.0 for m in stin_map.keys()}

    if seed_pi_ex and init_pi_amount > 0:
        ext_pool[seed_pi_ex] = ext_pool.get(seed_pi_ex, 0.0) + float(init_pi_amount)

    if met_polyp is not None:
        int_pool[met_polyp.id] = int_pool.get(met_polyp.id, 0.0) + float(INIT_POLYP)

    stout_polyp = (
        rxn_by_id(model, f"STOUT_{sanitize(met_polyp.id)}") if met_polyp is not None else None
    )
    stin_pi = (
        rxn_by_id(model, f"STIN_{sanitize(met_pi.id)}") if met_pi is not None else None
    )

    X = float(INIT_BIOMASS)
    Ms = float(dose) if ex_sub is not None else 0.0
    sub_used = 0.0

    pi_rel_ana     = 0.0
    pi_upt_aer     = 0.0
    polyp_hyd_ana  = 0.0
    polyp_syn_aer  = 0.0

    if ex_o2:
        safe_fix(ex_o2, 0.0)
    o2_cons = find_o2_consumers(model)
    saved_bounds = {}
    for rid in o2_cons:
        r = rxn_by_id(model, rid)
        if r:
            saved_bounds[rid] = (r.lower_bound, r.upper_bound)
            safe_box(r, 0.0, 0.0)
    saved_bio_bounds = (biomass.lower_bound, biomass.upper_bound)
    safe_fix(biomass, 0.0)

    if ex_sub is not None:
        obj = ex_sub.flux_expression
        if SOFT_SUPPRESS_GLY_ATP:
            for r in find_like_reactions(model, GLY_ATP_CAND_IDS):
                obj = obj + SOFT_EPS * r.flux_expression
        model.objective = model.problem.Objective(obj, direction="min")
    else:
        glc_rxns = get_group_reactions(model, ("subsys_GLYCOLYSIS",), ("GLYCOLYSIS",))
        if glc_rxns:
            expr = 0
            for r in glc_rxns:
                expr = expr + r.flux_expression
            model.objective = model.problem.Objective(expr, direction="max")

    if r_hyd is not None:
        lb = max(POLY_REQ_LB, r_hyd.lower_bound)
        safe_box(r_hyd, lb, 1000.0)
    if r_syn is not None:
        safe_box(r_syn, 0.0, 0.0)

    rec = []
    sum_mu_dt = 0.0
    sum_dt = 0.0

    # anaerobic phase
    for k in range(tA):
        if ex_sub:
            if k < tW and Ms > 1e-12:
                time_left = max((tW - k) * dt, 1e-9)
                need = min(sub_max_ana, Ms / time_left)
                force_uptake(ex_sub, need, forbid_secretion=True)
            else:
                safe_box(ex_sub, 0.0, 0.0)

        for stin_id in stin_map.values():
            r = rxn_by_id(model, stin_id)
            if r:
                safe_box(r, 0.0, 1000.0)
        for stout_id in stout_map.values():
            r = rxn_by_id(model, stout_id)
            if r:
                safe_box(r, 0.0, 0.0)

        if stout_polyp is not None and met_polyp is not None:
            amt = float(int_pool.get(met_polyp.id, 0.0))
            ub = min(INT_REL_MAX, amt / dt) if amt > 1e-12 else 0.0
            safe_box(stout_polyp, 0.0, ub)

        if stin_pi is not None:
            safe_box(stin_pi, 0.0, 0.0)

        if ex_pi is not None:
            safe_box(ex_pi, 0.0, 1000.0)

        sol = solve_pfba(model)
        status = getattr(sol, "status", "unknown")
        v_bio = float(sol.fluxes.get(biomass.id, 0.0)) if sol is not None else 0.0

        v_sub = float(sol.fluxes.get(ex_sub.id, 0.0)) if (sol is not None and ex_sub) else 0.0
        if status == "optimal" and ex_sub:
            take = flux_to_uptake(ex_sub, v_sub) * dt
            take = min(take, Ms)
            Ms -= take
            sub_used += take

        if status == "optimal":
            if ex_pi is not None:
                vpi = float(sol.fluxes.get(ex_pi.id, 0.0))
                if vpi > 1e-12:
                    pi_rel_ana += vpi * dt
            if r_hyd is not None:
                vh = float(sol.fluxes.get(r_hyd.id, 0.0))
                if vh > 1e-12:
                    polyp_hyd_ana += vh * dt

            for ex in model.exchanges:
                if ex_sub and ex.id == ex_sub.id:
                    continue
                if ex_o2 and ex.id == ex_o2.id:
                    continue
                v = float(sol.fluxes.get(ex.id, 0.0))
                if v > 1e-12:
                    ext_pool[ex.id] = ext_pool.get(ex.id, 0.0) + v * dt

            for mid, stin_id in stin_map.items():
                if stin_pi is not None and stin_id == stin_pi.id:
                    continue
                v = float(sol.fluxes.get(stin_id, 0.0))
                if v > 1e-12:
                    int_pool[mid] = int_pool.get(mid, 0.0) + v * dt

            if stout_polyp is not None and met_polyp is not None:
                vrel = float(sol.fluxes.get(stout_polyp.id, 0.0))
                if vrel > 1e-12:
                    int_pool[met_polyp.id] = max(
                        0.0, int_pool.get(met_polyp.id, 0.0) - vrel * dt
                    )

        sum_dt += dt
        rec.append(
            {
                "t": k * dt,
                "phase": "ana",
                "X": X,
                "S_left": Ms,
                "Sub_v": v_sub,
                "status": status,
            }
        )

    # aerobic phase
    for rid, (lb, ub) in saved_bounds.items():
        r = rxn_by_id(model, rid)
        if r:
            safe_box(r, lb, ub)
    if ex_o2:
        safe_fix(ex_o2, -abs(float(o2_supply)))
    biomass.lower_bound, biomass.upper_bound = saved_bio_bounds
    if ex_sub:
        safe_box(ex_sub, 0.0, 1000.0)
    model.objective = biomass
    model.objective_direction = "maximize"

    if r_syn is not None:
        lb = max(POLY_SYN_EPS, r_syn.lower_bound)
        safe_box(r_syn, lb, 1000.0)
    if r_hyd is not None:
        safe_box(r_hyd, 0.0, 0.0)
    if ex_pi is not None:
        safe_box(ex_pi, -1000.0, 1000.0)

    for k in range(tO):
        for ex_id, amt in ext_pool.items():
            r = rxn_by_id(model, ex_id)
            if not r:
                continue
            lb = -min(reup_max_aer, amt / dt) if amt > 1e-12 else 0.0
            safe_box(r, lb, 1000.0)
        for mid, amt in int_pool.items():
            r = rxn_by_id(model, f"STOUT_{sanitize(mid)}")
            if not r:
                continue
            ub = min(int_rel_max, amt / dt) if amt > 1e-12 else 0.0
            safe_box(r, 0.0, ub)

        sol_max = model.optimize()
        if getattr(sol_max, "status", "") == "optimal":
            sol = pfba(model)
        else:
            sol = sol_max

        status = getattr(sol, "status", "unknown")
        v_bio = float(sol.fluxes.get(biomass.id, 0.0)) if sol is not None else 0.0
        if status == "optimal":
            if ex_pi is not None:
                vpi = float(sol.fluxes.get(ex_pi.id, 0.0))
                upt = max(0.0, -vpi) * dt
                if upt > 1e-12:
                    pi_upt_aer += upt
            if r_syn is not None:
                vs = float(sol.fluxes.get(r_syn.id, 0.0))
                if vs > 1e-12:
                    polyp_syn_aer += vs * dt

            for ex_id, amt in list(ext_pool.items()):
                v = float(sol.fluxes.get(ex_id, 0.0))
                take = max(0.0, -v) * dt
                if take > 0:
                    ext_pool[ex_id] = max(0.0, amt - take)
            for mid, amt in list(int_pool.items()):
                v = float(sol.fluxes.get(f"STOUT_{sanitize(mid)}", 0.0))
                rel = max(0.0, v) * dt
                if rel > 1e-12:
                    int_pool[mid] = max(0.0, amt - rel)

            if v_bio > 0:
                X = max(0.0, X * (1.0 + v_bio * dt))
                sum_mu_dt += v_bio * dt

        sum_dt += dt
        rec.append(
            {
                "t": (tA + k) * dt,
                "phase": "aer",
                "X": X,
                "S_left": Ms,
                "Sub_v": 0.0,
                "status": status,
            }
        )

    mu_avg = (sum_mu_dt / sum_dt) if sum_dt > 0 else 0.0
    dX = max(0.0, X - INIT_BIOMASS)
    yxs = (dX / sub_used) if sub_used > 1e-12 else 0.0

    return (
        rec,
        X,
        mu_avg,
        sub_used,
        dX,
        yxs,
        pi_rel_ana,
        pi_upt_aer,
        polyp_hyd_ana,
        polyp_syn_aer,
    )
