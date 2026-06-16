#!/usr/bin/env bash
# BC-constraint sweep for the F-bar recovery comparison.
#
# Companion to fbar_strength_study.sh, which varies nu to change how much
# F-bar correction work is needed.  This script holds nu fixed and instead
# sweeps the boundary-constraint pattern -- moving from a heavily-locked
# domain (lots of faces pinned, big F-bar correction needed) toward a
# uniaxial-like setup (only symmetry planes pinned, lateral contraction
# free, much smaller F-bar correction).
#
# Four physics suites are supported via the MODE env var (default elastic_qs):
#   MODE=plastic     reference.i / restart.i                  (elasto-plastic + JC + T)
#   MODE=elastic     reference_elastic.i / restart_elastic.i  (elastic-only, dynamic)
#   MODE=elastic_qs  reference_elastic_qs.i /                 (elastic + quasi-static;
#                    restart_elastic_qs.i                      cleanest BC isolation
#                                                              -- no inertia or
#                                                              plasticity blurring
#                                                              the F-bar signal)
#   MODE=plastic_qs  reference_plastic_qs.i /                 (plastic + quasi-static;
#                    restart_plastic_qs.i                      J2 + JC + T-coupled but
#                                                              no inertia/Newmark)
# Each mode writes to its own default OUT_DIR (outputs_bc_sweep_<mode>/) so
# the sweeps don't collide.  Override OUT_DIR= to redirect any default.
#
# Each BC case is a named tuple (xfix_bnd, yfix_bnd, zfix_bnd) of boundary
# strings that the reference and restart both honor.  For each case the
# script runs:
#   reference + NEW restart (approach A, recover U_bar)
#   reference + OLD restart (approach B, recover raw U + init averaging)
# Same per-method recovery-dump strategy as fbar_strength_study.sh.
#
# Usage:
#   ./bc_sweep_study.sh                                       # elastic_qs default sweep
#   MODE=plastic ./bc_sweep_study.sh                          # elasto-plastic suite
#   MODE=elastic ./bc_sweep_study.sh                          # elastic-dynamic suite
#   MODE=plastic_qs ./bc_sweep_study.sh                       # plastic + quasi-static
#   BC_CASES="uniaxial constrained" ./bc_sweep_study.sh
#   NU=0.30 NP=4 N=1 ./bc_sweep_study.sh                      # different nu / mesh / ranks
#
# Output:
#   ${OUT_DIR}/bc_sweep_results_new.csv     (one row per BC case, method=new)
#   ${OUT_DIR}/bc_sweep_results_old.csv     (one row per BC case, method=old)
#   ${OUT_DIR}/<files>...                   (per-case reference + restart e/csv files)

set -euo pipefail

cd "$(dirname "$0")"

RACCOON="${RACCOON:-/home/det12/projects/raccoon/raccoon-opt}"
NU="${NU:-0.49}"
N="${N:-1}"
NP="${NP:-2}"
N_RESTART_STEPS="${N_RESTART_STEPS:-5}"

# MODE selects which physics suite to sweep the BC pattern through:
#   plastic     -> reference.i / restart.i              (elasto-plastic + JC + T-aware)
#   elastic     -> reference_elastic.i / restart_elastic.i  (elastic-only, dynamic)
#   elastic_qs  -> reference_elastic_qs.i / restart_elastic_qs.i  (elastic + quasi-static)
# Each mode picks the matching pair of input files and writes to its own
# default OUT_DIR so the three sweeps don't collide.
MODE="${MODE:-elastic_qs}"
case "$MODE" in
  plastic)
    REF_INPUT="reference.i"
    RES_INPUT="restart.i"
    DEFAULT_OUT_DIR="outputs_bc_sweep_plastic"
    ;;
  elastic)
    REF_INPUT="reference_elastic.i"
    RES_INPUT="restart_elastic.i"
    DEFAULT_OUT_DIR="outputs_bc_sweep_elastic"
    ;;
  elastic_qs)
    REF_INPUT="reference_elastic_qs.i"
    RES_INPUT="restart_elastic_qs.i"
    DEFAULT_OUT_DIR="outputs_bc_sweep_elastic_qs"
    ;;
  plastic_qs)
    REF_INPUT="reference_plastic_qs.i"
    RES_INPUT="restart_plastic_qs.i"
    DEFAULT_OUT_DIR="outputs_bc_sweep_plastic_qs"
    ;;
  *)
    echo "ERROR: MODE='$MODE' not recognized.  Use MODE=plastic, MODE=elastic, MODE=elastic_qs, or MODE=plastic_qs." >&2
    exit 1
    ;;
esac
OUT_DIR="${OUT_DIR:-$DEFAULT_OUT_DIR}"

# be_bar is only available when plasticity is active.  Drop it from the
# recovery-variable lists for any elastic-only mode.
case "$MODE" in
  elastic|elastic_qs)   BE_BAR_TOK="" ;;
  plastic|plastic_qs)   BE_BAR_TOK="be_bar " ;;
  *)                    BE_BAR_TOK="be_bar " ;;
esac

# Hardcoded BC cases.  Order matters only for the output table.
# - heavily      : all 6 face-sides clamped except the loaded one; even tighter
#                  than the original setup (right face also fixed in x).
# - constrained  : original reference setup (xfix on left+top, zfix front+back).
# - uniaxial     : symmetry planes only -- left (x=0), bottom (y=0), back (z=0).
#                  Right/front faces free to contract laterally.
read -r -a BC_CASES <<< "${BC_CASES:-heavily constrained uniaxial}"

declare -A XFIX_OF YFIX_OF ZFIX_OF DESC_OF
XFIX_OF[heavily]="left right top"
YFIX_OF[heavily]="bottom"
ZFIX_OF[heavily]="front back"
DESC_OF[heavily]="left+right+top  / bottom / front+back  (most constrained)"

XFIX_OF[constrained]="left top"
YFIX_OF[constrained]="bottom"
ZFIX_OF[constrained]="front back"
DESC_OF[constrained]="left+top        / bottom / front+back  (reference default)"

XFIX_OF[uniaxial]="left"
YFIX_OF[uniaxial]="bottom"
ZFIX_OF[uniaxial]="back"
DESC_OF[uniaxial]="left            / bottom / back        (symmetry only)"

# Input pair was picked above based on MODE.  All four pairs carry the
# xfix_bnd / yfix_bnd / zfix_bnd CLI hooks the BC sweep needs.

# MPI launcher autodetect.
if [[ -n "${MPIRUN:-}" ]]; then
  :
elif command -v mpirun >/dev/null 2>&1; then
  MPIRUN=mpirun
else
  MPIRUN="conda run -n moose mpirun"
fi

# Pull end_time / dt from the reference input so the comparison time stays
# correct if you ever re-tune the reference's loading.
REF_END_TIME=$(awk -F'=' '/^[[:space:]]*end_time[[:space:]]*=/{gsub(/[[:space:]]/,"",$2); print $2; exit}' "$REF_INPUT")
REF_END_TIME="${REF_END_TIME:-1.5}"
DT=$(awk -F'=' '/^[[:space:]]*dt[[:space:]]*=/{gsub(/[[:space:]]/,"",$2); print $2; exit}' "$REF_INPUT")
DT="${DT:-0.1}"

TARGET_TIME=$(python3 -c "print(${REF_END_TIME} + ${N_RESTART_STEPS}*${DT})")
FIRST_STEP_TIME=$(python3 -c "print(${REF_END_TIME} + ${DT})")
RECOVER_TIMESTEP="LATEST"

echo "BC-sweep study config:"
echo "  MODE             = $MODE"
echo "  raccoon-opt      = $RACCOON"
echo "  MPIRUN           = $MPIRUN"
echo "  inputs           = $REF_INPUT + $RES_INPUT"
echo "  OUT_DIR          = $OUT_DIR"
echo "  N (mesh)         = $N"
echo "  NP               = $NP"
echo "  NU (fixed)       = $NU"
echo "  BC_CASES         = ${BC_CASES[*]}"
echo "  ref_end_time     = $REF_END_TIME"
echo "  dt               = $DT"
echo "  N_RESTART_STEPS  = $N_RESTART_STEPS"
echo "  target_time      = $TARGET_TIME"
for case in "${BC_CASES[@]}"; do
  echo "    case '$case'   xfix='${XFIX_OF[$case]:-?}'  yfix='${YFIX_OF[$case]:-?}'  zfix='${ZFIX_OF[$case]:-?}'"
done

if [[ ! -x "$RACCOON" ]]; then
  echo "ERROR: raccoon-opt not found or not executable at $RACCOON" >&2
  exit 1
fi
mkdir -p "$OUT_DIR"

# ---------------------------------------------------------------------------
# extract_col -- same as in fbar_strength_study.sh
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

RESULTS_NEW="${OUT_DIR}/bc_sweep_results_new.csv"
RESULTS_OLD="${OUT_DIR}/bc_sweep_results_old.csv"
RESULTS_HEADER="bc_case,xfix,yfix,zfix,fbar_correction_int_ref,fbar_correction_max_ref,fbar_pressure_correction_int_ref,fbar_pressure_correction_max_ref,psie_ref_first,psie_restart_first,abs_err_first,rel_err_first,psie_ref_target,psie_restart_target,abs_err_target,rel_err_target"
echo "$RESULTS_HEADER" > "$RESULTS_NEW"
echo "$RESULTS_HEADER" > "$RESULTS_OLD"

# ---------------------------------------------------------------------------
# run_method <bc_case> <case_tag> <method> <recover_apply_fbar_to_U> <restart_tensor_mats> <ref_recover_tensor_mats>
#   Runs reference + restart for the given (bc_case, method).  be_bar is
#   absent in the elastic-QS variant so it never appears in the lists.
# ---------------------------------------------------------------------------
run_method() {
  local bc_case="$1"
  local case_tag="$2"
  local method="$3"
  local apply_fbar="$4"
  local restart_tensor_mats="$5"
  local ref_recover_tensor_mats="$6"
  local pair_tag="${case_tag}_${method}"
  local ref_csv="${OUT_DIR}/reference_fbar_out_${N}${pair_tag}.csv"
  local res_csv="${OUT_DIR}/restart_fbar_out_${N}${pair_tag}.csv"

  local xfix="${XFIX_OF[$bc_case]}"
  local yfix="${YFIX_OF[$bc_case]}"
  local zfix="${ZFIX_OF[$bc_case]}"

  echo "[bc=$bc_case method=$method] $REF_INPUT  (tag=$pair_tag)"
  local ref_cmd=(
    $MPIRUN -n "$NP" "$RACCOON" -i "$REF_INPUT"
    n=$N out_dir="$OUT_DIR" tag="$pair_tag" nu=$NU
    end_time=$TARGET_TIME dump_time=$REF_END_TIME
    "xfix_bnd=$xfix" "yfix_bnd=$yfix" "zfix_bnd=$zfix"
    "RecoverVariables/rec/tensor_materials=$ref_recover_tensor_mats"
    -w
  )
  printf '  $ '; printf '%q ' "${ref_cmd[@]}"; printf '\n'
  "${ref_cmd[@]}" | tail -n 4

  if [[ ! -f "$ref_csv" ]]; then
    echo "ERROR: $ref_csv not produced -- reference failed at bc=$bc_case method=$method" >&2
    exit 1
  fi

  echo "[bc=$bc_case method=$method] $RES_INPUT    (tag=$pair_tag, recover_apply_fbar_to_U=$apply_fbar)"
  local res_cmd=(
    $MPIRUN -n "$NP" "$RACCOON" -i "$RES_INPUT"
    n=$N out_dir="$OUT_DIR" tag="$pair_tag" output_tag="$pair_tag" nu=$NU
    start_time=$REF_END_TIME end_time=$TARGET_TIME
    recover_timestep=$RECOVER_TIMESTEP
    "xfix_bnd=$xfix" "yfix_bnd=$yfix" "zfix_bnd=$zfix"
    Materials/defgrad/recover_apply_fbar_to_U=$apply_fbar
    "UserObjects/epsol/tensor_materials=$restart_tensor_mats"
    -w
  )
  printf '  $ '; printf '%q ' "${res_cmd[@]}"; printf '\n'
  "${res_cmd[@]}" | tail -n 4

  if [[ ! -f "$res_csv" ]]; then
    echo "ERROR: $res_csv not produced -- restart failed at bc=$bc_case method=$method" >&2
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

  printf '    fbar_pressure_correction_max_ref = %s\n' "$fbar_p_max"
  printf '    first-step (t=%s):  abs_err=%s  rel_err=%s\n' "$FIRST_STEP_TIME" "$abs_first" "$rel_first"
  printf '    target     (t=%s):  abs_err=%s  rel_err=%s\n' "$TARGET_TIME"     "$abs_targ"  "$rel_targ"

  local results_file
  if [[ "$method" == "new" ]]; then
    results_file="$RESULTS_NEW"
  else
    results_file="$RESULTS_OLD"
  fi
  echo "$bc_case,\"$xfix\",\"$yfix\",\"$zfix\",$fbar_int,$fbar_max,$fbar_p_int,$fbar_p_max,$psie_ref_first,$psie_res_first,$abs_first,$rel_first,$psie_ref_targ,$psie_res_targ,$abs_targ,$rel_targ" >> "$results_file"
}

for bc_case in "${BC_CASES[@]}"; do
  if [[ -z "${XFIX_OF[$bc_case]:-}" ]]; then
    echo "ERROR: unknown BC case '$bc_case'.  Available: ${!XFIX_OF[*]}" >&2
    exit 1
  fi
  echo
  echo "================================================================"
  echo "  bc_case = $bc_case   (N=$N, NP=$NP, nu=$NU)"
  echo "  ${DESC_OF[$bc_case]}"
  echo "================================================================"

  CASE_TAG="_bc_${bc_case}"

  # ${BE_BAR_TOK} is "be_bar " for plastic mode and empty for elastic/elastic_qs.
  run_method "$bc_case" "$CASE_TAG" "new" "false" \
    "stress ${BE_BAR_TOK}stretch_tensor_fbar rotation_tensor" \
    "${BE_BAR_TOK}stress rotation_tensor stretch_tensor_fbar"

  run_method "$bc_case" "$CASE_TAG" "old" "true" \
    "stress ${BE_BAR_TOK}stretch_tensor rotation_tensor" \
    "${BE_BAR_TOK}stress rotation_tensor stretch_tensor"
done

echo
echo "================================================================"
echo "  BC-sweep summary -- NEW method (approach A)"
echo "================================================================"
column -ts, "$RESULTS_NEW"
echo
echo "================================================================"
echo "  BC-sweep summary -- OLD method (approach B)"
echo "================================================================"
column -ts, "$RESULTS_OLD"
echo
echo "Saved: $RESULTS_NEW"
echo "       $RESULTS_OLD"
echo
echo "Suggested reading: rel_err_target vs fbar_pressure_correction_max_ref"
echo "  ordered across bc_cases shows how recovery accuracy scales with BC-driven"
echo "  F-bar correction load, with nu held fixed at $NU."
