#!/usr/bin/env bash
# Nu sweep for the single-HEX8 incompressible bending recovery test.
#
# For each nu in NU_VALUES:
#   1. Run pull_test.i as REFERENCE end_time = TARGET_TIME, with dump_time
#      held at REF_END_TIME so the displaced exodusqp recovery dump captures
#      the kinematic state at REF_END_TIME.  Two separate reference runs per
#      nu, one per recovery method:
#        method "new" (approach A):  dumps stretch_tensor_fbar only
#        method "old" (approach B):  dumps stretch_tensor only
#      (libmesh's exodus reader shadows the shorter prefix when both
#       stretch_tensor_* and stretch_tensor_fbar_* are in the same file, so
#       each method gets its own reference run.)
#   2. Run pull_test_restart.i from REF_END_TIME to TARGET_TIME, consuming
#      that dump.  NEW restart uses Materials/defgrad/recover_apply_fbar_to_U
#      = false (consume U_bar directly).  OLD restart uses = true (rebuild
#      F_bar from raw U by averaging at INITIAL).
#
# The script emits two CSVs -- one per method -- with one row per nu containing
# F-bar diagnostic values from the reference at TARGET_TIME (x-axis: "amount
# of F-bar correction work") AND psie / J recovery error vs the reference
# (y-axis: "recovery accuracy").
#
# MODE selects which reference+restart pair the sweep drives:
#   elastic_qs (default):  pull_test.i              + pull_test_restart.i
#   dynamic            :   pull_test_dynamic.i      + pull_test_dynamic_restart.i
#   plastic            :   pull_test_plastic.i      + pull_test_plastic_restart.i
#   plastic_qs         :   pull_test_plastic_qs.i   + pull_test_plastic_qs_restart.i
# OUT_DIR defaults to outputs_nu_sweep_${MODE} so the four modes don't clobber.
#
# Usage:
#   ./nu_sweep_study.sh                                  # default sweep (elastic_qs)
#   MODE=dynamic ./nu_sweep_study.sh
#   MODE=plastic ./nu_sweep_study.sh
#   MODE=plastic_qs ./nu_sweep_study.sh
#   NU_VALUES="0.49 0.499 0.4999" ./nu_sweep_study.sh
#   N_RESTART_STEPS=10 ./nu_sweep_study.sh
#   OUT_DIR=runs/sweep1 ./nu_sweep_study.sh
#
# Output:
#   ${OUT_DIR}/nu_sweep_results_new.csv
#   ${OUT_DIR}/nu_sweep_results_old.csv

set -euo pipefail

cd "$(dirname "$0")"

RACCOON="${RACCOON:-/home/det12/projects/raccoon/raccoon-opt}"
MODE="${MODE:-elastic_qs}"
read -r -a NU_VALUES <<< "${NU_VALUES:-0.30 0.45 0.48 0.49 0.495 0.499 0.4999}"

NP="${NP:-1}"
# Default chosen so the restart's post-recovery time window stays 0.25
# (TARGET = 1.25, FIRST = 1.05) regardless of the dt halving in the .i files.
N_RESTART_STEPS="${N_RESTART_STEPS:-20}"

case "$MODE" in
  elastic_qs)
    REF_INPUT="pull_test.i"
    RES_INPUT="pull_test_restart.i"
    ;;
  dynamic)
    REF_INPUT="pull_test_dynamic.i"
    RES_INPUT="pull_test_dynamic_restart.i"
    ;;
  plastic)
    REF_INPUT="pull_test_plastic.i"
    RES_INPUT="pull_test_plastic_restart.i"
    ;;
  plastic_qs)
    REF_INPUT="pull_test_plastic_qs.i"
    RES_INPUT="pull_test_plastic_qs_restart.i"
    ;;
  *)
    echo "ERROR: unknown MODE='$MODE'; expected elastic_qs|dynamic|plastic|plastic_qs" >&2
    exit 1
    ;;
esac

OUT_DIR="${OUT_DIR:-outputs_nu_sweep_${MODE}}"

# MPI launcher autodetect.
if [[ -n "${MPIRUN:-}" ]]; then
  :
elif command -v mpirun >/dev/null 2>&1; then
  MPIRUN=mpirun
else
  MPIRUN="conda run -n moose mpirun"
fi

# Pull end_time and dt from pull_test.i.  These are the reference's natural
# loading window; REF_END_TIME (= the original end_time) is the time at which
# the displaced recovery dump should be written, and TARGET_TIME is the time
# the restart marches to (REF_END_TIME + N_RESTART_STEPS * dt).
REF_END_TIME=$(awk -F'=' '/^[[:space:]]*end_time[[:space:]]*=/{gsub(/[[:space:]]/,"",$2); print $2; exit}' "$REF_INPUT")
REF_END_TIME="${REF_END_TIME:-1.0}"
DT=$(awk -F'=' '/^[[:space:]]*dt[[:space:]]*=/{gsub(/[[:space:]]/,"",$2); print $2; exit}' "$REF_INPUT")
DT="${DT:-0.05}"

TARGET_TIME=$(python3 -c "print(${REF_END_TIME} + ${N_RESTART_STEPS}*${DT})")
FIRST_STEP_TIME=$(python3 -c "print(${REF_END_TIME} + ${DT})")
RECOVER_TIMESTEP="LATEST"

echo "Single-HEX8 nu-sweep config:"
echo "  raccoon-opt      = $RACCOON"
echo "  MPIRUN           = $MPIRUN"
echo "  MODE             = $MODE"
echo "  inputs           = $REF_INPUT + $RES_INPUT"
echo "  OUT_DIR          = $OUT_DIR"
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

# extract_col <csv_file> <column_name> <target_time>
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

RESULTS_NEW="${OUT_DIR}/nu_sweep_results_new.csv"
RESULTS_OLD="${OUT_DIR}/nu_sweep_results_old.csv"
RESULTS_HEADER="nu,K_over_G,fbar_correction_int_ref,fbar_correction_max_ref,fbar_pressure_correction_int_ref,fbar_pressure_correction_max_ref,psie_ref_first,psie_restart_first,abs_err_first,rel_err_first,psie_ref_target,psie_restart_target,abs_err_target,rel_err_target"
echo "$RESULTS_HEADER" > "$RESULTS_NEW"
echo "$RESULTS_HEADER" > "$RESULTS_OLD"

# run_method <nu> <nu_tag> <method> <recover_apply_fbar_to_U> <restart_tensor_mats> <ref_recover_tensor_mats>
#   Runs ONE reference + ONE restart for the given method.
run_method() {
  local nu="$1"
  local nu_tag="$2"
  local method="$3"
  local apply_fbar="$4"
  local restart_tensor_mats="$5"
  local ref_recover_tensor_mats="$6"
  local pair_tag="${nu_tag}_${method}"
  local ref_csv="${OUT_DIR}/pull_test_out${pair_tag}.csv"
  local res_csv="${OUT_DIR}/pull_test_restart_out${pair_tag}.csv"

  echo "[nu=$nu method=$method] $REF_INPUT  (tag=$pair_tag, dump_time=$REF_END_TIME, end_time=$TARGET_TIME)"
  local ref_cmd=(
    $MPIRUN -n "$NP" "$RACCOON" -i "$REF_INPUT"
    out_dir="$OUT_DIR" tag="$pair_tag" nu=$nu
    end_time=$TARGET_TIME dump_time=$REF_END_TIME
    "RecoverVariables/rec/tensor_materials=$ref_recover_tensor_mats"
  )
  printf '  $ '; printf '%q ' "${ref_cmd[@]}"; printf '\n'
  "${ref_cmd[@]}" | tail -n 4

  if [[ ! -f "$ref_csv" ]]; then
    echo "ERROR: $ref_csv not produced -- reference failed at nu=$nu method=$method" >&2
    exit 1
  fi

  echo "[nu=$nu method=$method] $RES_INPUT  (tag=$pair_tag, recover_apply_fbar_to_U=$apply_fbar)"
  local res_cmd=(
    $MPIRUN -n "$NP" "$RACCOON" -i "$RES_INPUT"
    out_dir="$OUT_DIR" tag="$pair_tag" output_tag="$pair_tag" nu=$nu
    start_time=$REF_END_TIME end_time=$TARGET_TIME end_time_for_ramp=$TARGET_TIME
    recover_timestep=$RECOVER_TIMESTEP
    Materials/defgrad/recover_apply_fbar_to_U=$apply_fbar
    "UserObjects/epsol/tensor_materials=$restart_tensor_mats" -w
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
  psie_ref_first=$(extract_col "$ref_csv" "psie_active_int" "$FIRST_STEP_TIME")
  psie_res_first=$(extract_col "$res_csv" "psie_active_int" "$FIRST_STEP_TIME")
  psie_ref_targ=$(extract_col "$ref_csv" "psie_active_int" "$TARGET_TIME")
  psie_res_targ=$(extract_col "$res_csv" "psie_active_int" "$TARGET_TIME")

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

  printf '    K/G                              = %s\n' "$k_over_g"
  printf '    fbar_pressure_correction_max_ref = %s\n' "$fbar_p_max"
  printf '    first-step (t=%s):  abs_err=%s  rel_err=%s\n' "$FIRST_STEP_TIME" "$abs_first" "$rel_first"
  printf '    target     (t=%s):  abs_err=%s  rel_err=%s\n' "$TARGET_TIME"     "$abs_targ"  "$rel_targ"

  local results_file
  if [[ "$method" == "new" ]]; then
    results_file="$RESULTS_NEW"
  else
    results_file="$RESULTS_OLD"
  fi
  echo "$nu,$k_over_g,$fbar_int,$fbar_max,$fbar_p_int,$fbar_p_max,$psie_ref_first,$psie_res_first,$abs_first,$rel_first,$psie_ref_targ,$psie_res_targ,$abs_targ,$rel_targ" >> "$results_file"
}

# Plastic modes need be_bar carried through both legs so the isochoric plastic
# state survives the dump/recover roundtrip.  Elastic modes don't have a
# be_bar material so adding it would error.
case "$MODE" in
  plastic|plastic_qs)
    EXTRA_TENSOR_MATS="be_bar"
    EXTRA_TENSOR_MATS_REF="be_bar"
    ;;
  *)
    EXTRA_TENSOR_MATS=""
    EXTRA_TENSOR_MATS_REF=""
    ;;
esac

for nu in "${NU_VALUES[@]}"; do
  echo
  echo "================================================================"
  echo "  nu = $nu   (NP=$NP)"
  echo "================================================================"

  NU_TAG="_nu_${nu//./p}"

  # Method "new": pre-averaged U_bar carried over, no init averaging.
  #   Reference dumps only stretch_tensor_fbar.
  run_method "$nu" "$NU_TAG" "new" "false" \
    "stress ${EXTRA_TENSOR_MATS} stretch_tensor_fbar rotation_tensor" \
    "stress ${EXTRA_TENSOR_MATS_REF} rotation_tensor stretch_tensor_fbar"

  # Method "old": recover raw U, apply F-bar averaging at init.
  #   Reference dumps only stretch_tensor.
  run_method "$nu" "$NU_TAG" "old" "true" \
    "stress ${EXTRA_TENSOR_MATS} stretch_tensor rotation_tensor" \
    "stress ${EXTRA_TENSOR_MATS_REF} rotation_tensor stretch_tensor"
done

echo
echo "================================================================"
echo "  Nu sweep summary -- NEW method (approach A)"
echo "================================================================"
column -ts, "$RESULTS_NEW"
echo
echo "================================================================"
echo "  Nu sweep summary -- OLD method (approach B)"
echo "================================================================"
column -ts, "$RESULTS_OLD"
echo
echo "Saved: $RESULTS_NEW"
echo "       $RESULTS_OLD"
