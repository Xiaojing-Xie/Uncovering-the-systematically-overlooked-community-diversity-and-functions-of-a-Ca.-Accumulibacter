# run_batch.py

from config_utils import (
    NGAM_FLUX,
    CPD_O2,
    CPD_PI,
    CARBON_CPDS,
    AA_MODE,
    AUXO_TABLE,
    ensure_ngam_reaction,
    robust_biomass_internal,
    open_amino_acids,
    find_exchange_for_cpd,
    carbon_number_from_ex,
)
from two_phase_simulator import (
    ANA_TIME,
    AER_TIME,
    DT,
    FEED_CUTOFF,
    DOSE_MOLAR,
    SUB_MAX_ANA,
    REUP_MAX_AER,
    INT_REL_MAX,
    O2_SUPPLY,
    INIT_PI_POOL,
    integrate_two_phase,
)

SBML_ROOT = ""
MODELS = None


def main():
    warnings.filterwarnings("ignore")

    models = MODELS or [
        d
        for d in sorted(os.listdir(SBML_ROOT))
        if os.path.exists(os.path.join(SBML_ROOT, d, f"{d}.xml"))
    ]

    auxo_df = None
    if AA_MODE == "auto" and os.path.exists(AUXO_TABLE):
        try:
            auxo_df = pd.read_csv(AUXO_TABLE, index_col=0)
        except Exception as e:
            print(f"Warning: cannot read {AUXO_TABLE}: {e}")

    all_rows = []

    for m in models:
        xml = os.path.join(SBML_ROOT, m, f"{m}.xml")
        try:
            base = cobra.io.read_sbml_model(xml)
        except Exception as e:
            print(f"[{m}] SBML load failed: {e}")
            continue

        outdir = os.path.join("out", m)
        os.makedirs(outdir, exist_ok=True)

        try:
            bio_id = robust_biomass_internal(base)
        except Exception as e:
            print(f"[{m}] biomass detection failed: {e}")
            continue

        ex_o2 = find_exchange_for_cpd(base, CPD_O2)
        ex_pi = find_exchange_for_cpd(base, CPD_PI)

        # background only (no external carbon)
        for cname, seed in CARBON_CPDS.items():
            model_bg = cobra.io.read_sbml_model(xml)
            ensure_ngam_reaction(model_bg, NGAM_FLUX)
            open_amino_acids(model_bg, AA_MODE, auxo_df=auxo_df, model_name=m)

            try:
                recB, XB, muB, subB, dXB, yxsB, *_ = integrate_two_phase(
                    model_bg,
                    ex_sub_id=None,
                    ex_o2_id=ex_o2,
                    bio_id=robust_biomass_internal(model_bg),
                    ana_time=ANA_TIME,
                    aer_time=AER_TIME,
                    dt=DT,
                    dose=0.0,
                    feed_cut=0.0,
                    sub_max_ana=SUB_MAX_ANA,
                    reup_max_aer=REUP_MAX_AER,
                    int_rel_max=INT_REL_MAX,
                    o2_supply=O2_SUPPLY,
                    seed_pi_ex=None,
                    init_pi_amount=0.0,
                    sub_max_ANA=SUB_MAX_ANA,
                    reup_max_AER=REUP_MAX_AER,
                    int_rel_MAX=INT_REL_MAX,
                )
            except Exception as e:
                print(f"[{m} | {cname}] background run failed: {e}")
                dXB = 0.0
                recB = None

            if recB is not None:
                try:
                    pd.DataFrame(recB).to_csv(
                        os.path.join(outdir, f"{m}_{cname}_BACKGROUND_timecourse.csv"),
                        index=False,
                    )
                except Exception:
                    pass

            # main trajectory with carbon source
            model = cobra.io.read_sbml_model(xml)
            ensure_ngam_reaction(model, NGAM_FLUX)
            open_amino_acids(model, AA_MODE, auxo_df=auxo_df, model_name=m)

            ex_sub_id = find_exchange_for_cpd(model, seed)
            if ex_sub_id is None:
                row = {
                    "Model": m,
                    "Carbon": cname,
                    "Mu_avg_h-1": 0.0,
                    "DeltaX": 0.0,
                    "Sub_used": 0.0,
                    "Yx_per_mmol": 0.0,
                    "Yx_per_Cmmol": 0.0,
                    "DeltaX_net": 0.0,
                    "YxS_net_mmol": 0.0,
                    "YxS_net_Cmmol": 0.0,
                    "Note": "NO carbon EX",
                }
                all_rows.append(row)
                print(f"[{m} | {cname}] no carbon exchange identified")
                continue

            try:
                (
                    rec,
                    X_end,
                    mu_avg,
                    sub_used,
                    dX,
                    yxs,
                    _pi_rel_ana,
                    _pi_upt_aer,
                    _polyp_hyd_ana,
                    _polyp_syn_aer,
                ) = integrate_two_phase(
                    model,
                    ex_sub_id=ex_sub_id,
                    ex_o2_id=ex_o2,
                    bio_id=bio_id,
                    ana_time=ANA_TIME,
                    aer_time=AER_TIME,
                    dt=DT,
                    dose=DOSE_MOLAR,
                    feed_cut=FEED_CUTOFF,
                    sub_max_ana=SUB_MAX_ANA,
                    reup_max_aer=REUP_MAX_AER,
                    int_rel_max=INT_REL_MAX,
                    o2_supply=O2_SUPPLY,
                    seed_pi_ex=ex_pi,
                    init_pi_amount=INIT_PI_POOL,
                    sub_max_ANA=SUB_MAX_ANA,
                    reup_max_AER=REUP_MAX_AER,
                    int_rel_MAX=INT_REL_MAX,
                )
            except Exception as e:
                row = {
                    "Model": m,
                    "Carbon": cname,
                    "Mu_avg_h-1": 0.0,
                    "DeltaX": 0.0,
                    "Sub_used": 0.0,
                    "Yx_per_mmol": 0.0,
                    "Yx_per_Cmmol": 0.0,
                    "DeltaX_net": 0.0,
                    "YxS_net_mmol": 0.0,
                    "YxS_net_Cmmol": 0.0,
                    "Note": f"Error:{e}",
                }
                all_rows.append(row)
                print(f"[{m} | {cname}] error during main run")
                continue

            dX_net = max(0.0, dX - max(0.0, dXB))
            yxs_net = (dX_net / sub_used) if sub_used > 1e-12 else 0.0

            Cn = carbon_number_from_ex(model, ex_sub_id)
            if Cn and sub_used > 1e-12:
                sub_Cmmol = sub_used * Cn
                yxs_C = dX / sub_Cmmol
                yxs_net_C = dX_net / sub_Cmmol
            else:
                yxs_C = 0.0
                yxs_net_C = 0.0

            try:
                pd.DataFrame(rec).to_csv(
                    os.path.join(outdir, f"{m}_{cname}_timecourse.csv"),
                    index=False,
                )
            except Exception:
                pass

            note = "OK" if sub_used > 1e-12 else "NO anaerobic uptake"
            row = {
                "Model": m,
                "Carbon": cname,
                "Mu_avg_h-1": mu_avg,
                "DeltaX": dX,
                "Sub_used": sub_used,
                "Yx_per_mmol": yxs,
                "Yx_per_Cmmol": yxs_C,
                "DeltaX_net": dX_net,
                "YxS_net_mmol": yxs_net,
                "YxS_net_Cmmol": yxs_net_C,
                "Note": note,
            }
            all_rows.append(row)

            print(
                f"[{m} | {cname}] growth_avg={mu_avg:.3f}, "
                f"DeltaX={dX:.3f}, S_used={sub_used:.2f}, "
                f"yield={yxs:.2f}/mmol ({yxs_C:.2f}/Cmmol), status={note}"
            )

        # per-model summary
        dfm = pd.DataFrame([r for r in all_rows if r["Model"] == m])
        if not dfm.empty:
            dfm.to_csv(os.path.join(outdir, f"{m}_summary.csv"), index=False)

    # global wide tables
    if all_rows:
        df_all = pd.DataFrame(all_rows)
        df_all.to_csv("ebpr_summary_all.csv", index=False)
        carb_order = list(CARBON_CPDS.keys())

        wide_m = (
            df_all.pivot(index="Model", columns="Carbon", values="Yx_per_mmol")
            .reindex(columns=carb_order)
            .fillna(0.0)
        )
        wide_c = (
            df_all.pivot(index="Model", columns="Carbon", values="Yx_per_Cmmol")
            .reindex(columns=carb_order)
            .fillna(0.0)
        )
        wide_nm = (
            df_all.pivot(index="Model", columns="Carbon", values="YxS_net_mmol")
            .reindex(columns=carb_order)
            .fillna(0.0)
        )
        wide_nc = (
            df_all.pivot(index="Model", columns="Carbon", values="YxS_net_Cmmol")
            .reindex(columns=carb_order)
            .fillna(0.0)
        )

        wide_m.to_csv("ebpr_yxs_wide.csv")
        wide_c.to_csv("ebpr_yxsC_wide.csv")
        wide_nm.to_csv("ebpr_yxs_net_wide.csv")
        wide_nc.to_csv("ebpr_yxsC_net_wide.csv")


if __name__ == "__main__":
    main()
