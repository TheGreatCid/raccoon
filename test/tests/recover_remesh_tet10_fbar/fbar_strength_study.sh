#!/usr/bin/env bash
# F-bar averaging strength sweep at fixed mesh refinement.
#
# Holds mesh refinement at n=1 (24 TET10 elements) and varies the Poisson
# ratio nu, which controls the bulk-to-shear modulus ratio K/G:
#
#   nu=0.30  -> K/G ~  2.2   (compressible, F-bar barely active)
#   nu=0.40  -> K/G ~  3.3
#   nu=0.45  -> K/G ~  4.8
#   nu=0.48  -> K/G ~  8.7
#   nu=0.49  -> K/G ~ 16.5   (near-incompressible, F-bar critical)
#   nu=0.495 -> K/G ~ 33     (very near-incompressible)
#
# Higher nu -> more volumetric locking -> F-bar correction does more work per
# element.  The reference's fbar_correction_max/int diagnostic is the actual
# measured "amount of averaging correction" for each nu (read straight from
# the reference CSV at target_time).  For each nu the script:
#
#   1. runs reference.i ONCE with that nu, n=1, extended by N_RESTART_STEPS past
#      ref_end_time so the file has a TIMESTEP_END row at target_time, and the
#      displaced recovery dump (sync_times = ref_end_time) contains BOTH the
#      raw stretch_tensor AND the pre-averaged stretch_tensor_fbar.
#   2. runs TWO restart legs from that same reference dump:
#        - method "new"  (approach A): consume the pre-averaged U_bar
#          (stretch_tensor_fbar) directly from the recovery file; no initial
#          volume averaging is applied on restart.
#        - method "old"  (approach B): consume the raw recovered U
#          (stretch_tensor) and apply the cfb539fe1-era F-bar volumetric
#          averaging operator at INITIAL to rebuild F_bar on the new mesh.
#      Each method produces its own labeled .csv / .e files
#      (output_tag = ${nu_tag}_${method}).
#   3. records reference fbar_correction_int / _max / pressure variants
#      (x-axis: "amount of F-bar work") AND the psie recovery error
#      (y-axis: "recovery accuracy") in TWO separate result files -- one row
#      per nu in each:
#        fbar_strength_results_new.csv  (approach A, recover U_bar)
#        fbar_strength_results_old.csv  (approach B, recover raw U + init avg)
#      Same column layout in both; same nu rows appear in both for direct
#      side-by-side plotting.
#
# Plot rel_err_target vs fbar_pressure_correction_max_ref from the results
# CSV grouped by method to see (a) whether recovery accuracy degrades as the
# volumetric correction work grows, and (b) whether OLD (recover-then-average)
# and NEW (recover-pre-averaged) behave differently in that limit.
#
# Usage:
#   ./fbar_strength_study.sh                                  # plastic (default)
#   MODE=elastic ./fbar_strength_study.sh                     # elastic-only suite
#   NU_VALUES="0.30 0.45 0.49" ./fbar_strength_study.sh
#   N_RESTART_STEPS=10 ./fbar_strength_study.sh
#   N=2 NP=4 ./fbar_strength_study.sh                         # different mesh
#
# Modes:
#   MODE=plastic (default)  reference.i / restart.i           full elasto-plastic
#                                                             test (JC hardening,
#                                                             T-coupled).  OUT_DIR
#                                                             defaults to
#                                                             outputs_fbar/
#   MODE=elastic            reference_elastic.i /             elastic-only test,
#                           restart_elastic.i                 no plasticity, no T,
#                                                             no JC.  OUT_DIR
#                                                             defaults to
#                                                             outputs_fbar_elastic/
#   MODE=elastic_qs         reference_elastic_qs.i /          elastic + quasi-static
#                           restart_elastic_qs.i              (no inertia, no
#                                                             Newmark, no HHT, no T,
#                                                             no plasticity).
#                                                             OUT_DIR defaults to
#                                                             outputs_fbar_elastic_qs/
#   MODE=plastic_qs         reference_plastic_qs.i /          plastic + quasi-static
#                           restart_plastic_qs.i              (J2 + JC + T-coupled,
#                                                             no inertia, no Newmark).
#                                                             OUT_DIR defaults to
#                                                             outputs_fbar_plastic_qs/
# Override OUT_DIR= to change any default.
#
# Output:
#   ${OUT_DIR}/fbar_strength_results_new.csv   (one row per nu, method=new)
#   ${OUT_DIR}/fbar_strength_results_old.csv   (one row per nu, method=old)

set -euo pipefail

cd "$(dirname "$0")"

RACCOON="${RACCOON:-/home/det12/projects/raccoon/raccoon-opt}"
read -r -a NU_VALUES <<< "${NU_VALUES:-0.30 0.40 0.45 0.48 0.49 0.495}"

# MODE selects which physics suite to sweep:
#   plastic (default)  -> reference.i / restart.i           (elastoplastic + JC + T-aware)
#   elastic            -> reference_elastic.i / restart_elastic.i  (elastic-only, no T)
# Each mode writes to its own default OUT_DIR so the two sweeps don't collide.
MODE="${MODE:-plastic}"
case "$MODE" in
  plastic)
    REF_INPUT="reference.i"
    RES_INPUT="restart.i"
    DEFAULT_OUT_DIR="outputs_fbar"
    ;;
  elastic)
    REF_INPUT="reference_elastic.i"
    RES_INPUT="restart_elastic.i"
    DEFAULT_OUT_DIR="outputs_fbar_elastic"
    ;;
  elastic_qs)
    REF_INPUT="reference_elastic_qs.i"
    RES_INPUT="restart_elastic_qs.i"
    DEFAULT_OUT_DIR="outputs_fbar_elastic_qs"
    ;;
  plastic_qs)
    REF_INPUT="reference_plastic_qs.i"
    RES_INPUT="restart_plastic_qs.i"
    DEFAULT_OUT_DIR="outputs_fbar_plastic_qs"
    ;;
  *)
    echo "ERROR: MODE='$MODE' is not recognized.  Use MODE=plastic, MODE=elastic, MODE=elastic_qs, or MODE=plastic_qs." >&2
    exit 1
    ;;
esac
OUT_DIR="${OUT_DIR:-$DEFAULT_OUT_DIR}"

# Fixed mesh refinement and rank count.  Defaults match the n=1 case in the
# convergence study.  Override with N=2 NP=2 etc. if you want to repeat the
# nu sweep at a different mesh size.
N="${N:-1}"
NP="${NP:-2}"

# Steps to march past ref_end_time in the restart (same semantics as
# mesh_convergence.sh).  Defaults to 5.
N_RESTART_STEPS="${N_RESTART_STEPS:-5}"

# MPI launcher autodetect.
if [[ -n "${MPIRUN:-}" ]]; then
  : # honor override
elif command -v mpirun >/dev/null 2>&1; then
  MPIRUN=mpirun
else
  MPIRUN="conda run -n moose mpirun"
fi

# Pull end_time and dt out of reference.i.
REF_END_TIME=$(awk -F'=' '/^[[:space:]]*end_time[[:space:]]*=/{gsub(/[[:space:]]/,"",$2); print $2; exit}' "$REF_INPUT")
REF_END_TIME="${REF_END_TIME:-1.5}"
DT=$(awk -F'=' '/^[[:space:]]*dt[[:space:]]*=/{gsub(/[[:space:]]/,"",$2); print $2; exit}' "$REF_INPUT")
DT="${DT:-0.1}"

TARGET_TIME=$(python3 -c "print(${REF_END_TIME} + ${N_RESTART_STEPS}*${DT})")
FIRST_STEP_TIME=$(python3 -c "print(${REF_END_TIME} + ${DT})")
# reference.i's exodusqp now uses sync_times=dump_time + sync_only, so the
# displaced recovery file holds exactly one slice at dump_time (= REF_END_TIME
# below).  LATEST always picks the right one.
RECOVER_TIMESTEP="LATEST"

echo "F-bar strength study config:"
echo "  raccoon-opt      = $RACCOON"
echo "  MPIRUN           = $MPIRUN"
echo "  MODE             = $MODE   ($REF_INPUT + $RES_INPUT)"
echo "  OUT_DIR          = $OUT_DIR"
echo "  N (fixed mesh)   = $N"
echo "  NP               = $NP"
echo "  NU_VALUES        = ${NU_VALUES[*]}"
echo "  ref_end_time     = $REF_END_TIME"
echo "  dt               = $DT"
echo "  N_RESTART_STEPS  = $N_RESTART_STEPS"
echo "  target_time      = $TARGET_TIME"
echo "  recover_timestep = $RECOVER_TIMESTEP"

if [[ ! -x "$RACCOON" ]]; then
  echo "ERROR: raccoon-opt not found or not executable at $RACCOON" >&2
  exit 1
fi

mkdir -p "$OUT_DIR"

# ---------------------------------------------------------------------------
# extract_col <csv_file> <column_name> <target_time>
#   Returns the value of <column_name> at the row whose "time" is closest to
#   <target_time>.  Same locator logic as in mesh_convergence.sh, but with
#   the column name passed in.
# ---------------------------------------------------------------------------
extract_col() {
  local csvfile="$1"
  local col_name="$2"
  local target="$3"
  if [[ ! -f "$csvfile" ]]; then
    echo "NaN"
    return 1
  fi
  awk -F, -v target="$target" -v want="$col_name" '
    NR==1 {
      for (i=1; i<=NF; i++) {
        h=$i; gsub(/[[:space:]\r"]/,"",h)
        col[h] = i
      }
      next
    }
    {
      t = $col["time"] + 0.0
      v = $col[want] + 0.0
      diff = (t > target) ? t - target : target - t
      if (NR==2 || diff < best_diff) {
        best_diff = diff
        best_v    = v
      }
    }
    END { printf "%.16e\n", best_v }
  ' "$csvfile"
}

RESULTS_NEW="${OUT_DIR}/fbar_strength_results_new.csv"
RESULTS_OLD="${OUT_DIR}/fbar_strength_results_old.csv"
# Each nu produces TWO rows -- one per recovery method:
#   method = "old"  (approach B):
#       Restart consumes the recovered RAW U (`stretch_tensor`) and applies the
#       cfb539fe1-era F-bar volumetric averaging operator at INITIAL to rebuild
#       F_bar on the new mesh.  Toggled via:
#           Materials/defgrad/recover_apply_fbar_to_U = true
#           UserObjects/epsol/tensor_materials      =
#               'stress be_bar stretch_tensor rotation_tensor'
#   method = "new"  (approach A; default in restart.i):
#       Restart consumes the pre-averaged U_bar (`stretch_tensor_fbar`) directly
#       from the recovery file -- no initial volume averaging is performed; the
#       averaging that produced U_bar was already done on the reference side.
#       Toggled via:
#           Materials/defgrad/recover_apply_fbar_to_U = false (default)
#           UserObjects/epsol/tensor_materials      =
#               'stress be_bar stretch_tensor_fbar rotation_tensor' (default)
# Both methods recover from the SAME reference dump (the reference's
# RecoverVariables/rec block dumps both stretch_tensor and stretch_tensor_fbar),
# so the reference is run only once per nu.
RESULTS_HEADER="nu,K_over_G,fbar_correction_int_ref,fbar_correction_max_ref,fbar_pressure_correction_int_ref,fbar_pressure_correction_max_ref,psie_ref_first,psie_restart_first,abs_err_first,rel_err_first,psie_ref_target,psie_restart_target,abs_err_target,rel_err_target"
echo "$RESULTS_HEADER" > "$RESULTS_NEW"
echo "$RESULTS_HEADER" > "$RESULTS_OLD"

# ---------------------------------------------------------------------------
# run_method <nu> <nu_tag> <method> <recover_apply_fbar_to_U> <restart_tensor_mats> <ref_recover_tensor_mats>
#   Runs ONE reference + ONE restart for the given method.  Each method gets
#   its own reference run because libmesh's Exodus reader silently shadows
#   `stretch_tensor_*` when `stretch_tensor_fbar_*` is present in the same
#   file -- so the reference must dump ONLY the variant the restart needs.
#   The reference writes to tag=${nu_tag}_${method}; the restart writes its
#   own outputs to output_tag=${nu_tag}_${method} and reads recover_file
#   tagged with the same suffix.  Appends one row to RESULTS.
# ---------------------------------------------------------------------------
run_method() {
  local nu="$1"
  local nu_tag="$2"
  local method="$3"
  local apply_fbar="$4"
  local restart_tensor_mats="$5"
  local ref_recover_tensor_mats="$6"
  local pair_tag="${nu_tag}_${method}"
  local ref_csv="${OUT_DIR}/reference_fbar_out_${N}${pair_tag}.csv"
  local res_csv="${OUT_DIR}/restart_fbar_out_${N}${pair_tag}.csv"

  echo "[nu=$nu method=$method] $REF_INPUT  (tag=$pair_tag)"
  # Build the reference command as an array so we can echo it before running.
  # `printf '%q '` quotes each token so the printed line is copy-paste-runnable
  # and shows the exact argv that gets dispatched (including spaces inside
  # multi-word parameter values).
  local ref_cmd=(
    $MPIRUN -n "$NP" "$RACCOON" -i $REF_INPUT
    n=$N out_dir=$OUT_DIR tag=$pair_tag nu=$nu
    end_time=$TARGET_TIME dump_time=$REF_END_TIME
    "RecoverVariables/rec/tensor_materials=$ref_recover_tensor_mats"
    -w
  )
  printf '  $ '; printf '%q ' "${ref_cmd[@]}"; printf '\n'
  "${ref_cmd[@]}" | tail -n 4

  if [[ ! -f "$ref_csv" ]]; then
    echo "ERROR: $ref_csv not produced -- reference failed at nu=$nu method=$method" >&2
    exit 1
  fi

  echo "[nu=$nu method=$method] $RES_INPUT    (tag=$pair_tag, recover_apply_fbar_to_U=$apply_fbar)"
  local res_cmd=(
    $MPIRUN -n "$NP" "$RACCOON" -i $RES_INPUT
    n=$N out_dir=$OUT_DIR tag=$pair_tag output_tag=$pair_tag nu=$nu
    start_time=$REF_END_TIME end_time=$TARGET_TIME
    recover_timestep=$RECOVER_TIMESTEP
    Materials/defgrad/recover_apply_fbar_to_U=$apply_fbar
    "UserObjects/epsol/tensor_materials=$restart_tensor_mats"
    -w
  )
  printf '  $ '; printf '%q ' "${res_cmd[@]}"; printf '\n'
  "${res_cmd[@]}" | tail -n 4

  if [[ ! -f "$res_csv" ]]; then
    echo "ERROR: $res_csv not produced -- restart failed at nu=$nu method=$method" >&2
    exit 1
  fi

  local fbar_int fbar_max fbar_p_int fbar_p_max
  fbar_int=$(extract_col "$ref_csv" "fbar_correction_int" "$TARGET_TIME")
  fbar_max=$(extract_col "$ref_csv" "fbar_correction_max" "$TARGET_TIME")
  fbar_p_int=$(extract_col "$ref_csv" "fbar_pressure_correction_int" "$TARGET_TIME")
  fbar_p_max=$(extract_col "$ref_csv" "fbar_pressure_correction_max" "$TARGET_TIME")

  local psie_ref_first psie_res_first psie_ref_targ psie_res_targ
  psie_ref_first=$(extract_col "$ref_csv" "psie_corr_active_int" "$FIRST_STEP_TIME")
  psie_res_first=$(extract_col "$res_csv" "psie_corr_active_int" "$FIRST_STEP_TIME")
  psie_ref_targ=$(extract_col "$ref_csv" "psie_corr_active_int" "$TARGET_TIME")
  psie_res_targ=$(extract_col "$res_csv" "psie_corr_active_int" "$TARGET_TIME")

  local abs_first rel_first abs_targ rel_targ
  read abs_first rel_first abs_targ rel_targ < <(python3 -c "
rf, rsf = ${psie_ref_first}, ${psie_res_first}
rt, rst = ${psie_ref_targ},  ${psie_res_targ}
def pair(ref, res):
    abs_e = abs(res - ref)
    rel_e = abs_e / max(abs(ref), 1e-30)
    return f'{abs_e:.16e} {rel_e:.16e}'
print(pair(rf, rsf), pair(rt, rst))
")

  local k_over_g
  k_over_g=$(python3 -c "nu=${nu}; print(2*(1+nu)/(3*(1-2*nu)))")

  printf '    first-step (t=%s):  abs_err=%s  rel_err=%s\n' \
         "$FIRST_STEP_TIME" "$abs_first" "$rel_first"
  printf '    target     (t=%s):  abs_err=%s  rel_err=%s\n' \
         "$TARGET_TIME"     "$abs_targ"  "$rel_targ"

  local results_file
  if [[ "$method" == "new" ]]; then
    results_file="$RESULTS_NEW"
  else
    results_file="$RESULTS_OLD"
  fi
  echo "$nu,$k_over_g,$fbar_int,$fbar_max,$fbar_p_int,$fbar_p_max,$psie_ref_first,$psie_res_first,$abs_first,$rel_first,$psie_ref_targ,$psie_res_targ,$abs_targ,$rel_targ" >> "$results_file"
}

for nu in "${NU_VALUES[@]}"; do
  echo
  echo "================================================================"
  echo "  nu = $nu   (n = $N, np = $NP)"
  echo "================================================================"

  k_over_g=$(python3 -c "nu=${nu}; print(2*(1+nu)/(3*(1-2*nu)))")
  echo "  K/G = $k_over_g"

  NU_TAG="_nu_${nu//./p}"

  # be_bar is only available when plasticity is active.  Drop it from the
  # recovery-variable lists in any elastic-only mode.
  case "$MODE" in
    elastic|elastic_qs)   be_bar_tok="" ;;
    plastic|plastic_qs)   be_bar_tok="be_bar " ;;
    *)                    be_bar_tok="be_bar " ;;
  esac

  # Method "new": pre-averaged U_bar carried over, no init averaging.
  #   Reference dumps only `stretch_tensor_fbar` (not `stretch_tensor`) so the
  #   restart's epsol UO can load the fbar variant without prefix shadowing.
  run_method "$nu" "$NU_TAG" "new" "false" \
    "stress ${be_bar_tok}stretch_tensor_fbar rotation_tensor" \
    "${be_bar_tok}stress rotation_tensor stretch_tensor_fbar"

  # Method "old": recover raw U, apply F-bar averaging at init.
  #   Reference dumps only `stretch_tensor` (not `stretch_tensor_fbar`).
  run_method "$nu" "$NU_TAG" "old" "true" \
    "stress ${be_bar_tok}stretch_tensor rotation_tensor" \
    "${be_bar_tok}stress rotation_tensor stretch_tensor"
done

echo
echo "================================================================"
echo "  F-bar strength summary -- NEW method (approach A)"
echo "================================================================"
column -ts, "$RESULTS_NEW"
echo
echo "================================================================"
echo "  F-bar strength summary -- OLD method (approach B)"
echo "================================================================"
column -ts, "$RESULTS_OLD"
echo
echo "Saved: $RESULTS_NEW"
echo "       $RESULTS_OLD"
echo
echo "Suggested plots:"
echo "  rel_err_target vs fbar_pressure_correction_max_ref, grouped by method"
echo "    --> does recovery accuracy degrade with K-scaled correction work,"
echo "        and is the trend different for OLD (init-averaging) vs"
echo "        NEW (pre-averaged-U_bar) recovery?"
