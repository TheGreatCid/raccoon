#!/usr/bin/env bash
# Dump-time sweep for the single-HEX8 incompressible bending recovery test.
#
# Sister script to nu_sweep_study.sh.  Instead of varying nu (which scales K
# and thus the F-bar pressure-correction magnitude), this sweep fixes nu and
# varies the *physical time at which the dump is written*.  Because the BC
# ramps the bending displacement linearly in t, a later dump_time means more
# accumulated deformation, which means more pointwise volumetric work for
# F-bar to do.  The x-axis of the resulting plots is the F-bar pressure
# correction magnitude *at the dump moment*, swept by dump time rather than
# by Poisson ratio.
#
# To keep the loading rate (and thus the pre-dump trajectory) identical across
# every dump_time, the reference is run with end_time_for_ramp pinned to a
# constant (default 1.0) -- that's what the BC formula
#   target_uy = pull_amount * (t / end_time_for_ramp) * X0
# uses.  end_time itself varies per dump_time to give the reference enough
# room to walk past the dump and compare to the restart at TARGET_TIME.
#
# MODE selects which reference+restart pair the sweep drives, same set as
# nu_sweep_study.sh:
#   elastic_qs (default), dynamic, plastic, plastic_qs
# OUT_DIR defaults to outputs_dump_time_sweep_${MODE}.
#
# Usage:
#   ./dump_time_sweep_study.sh                                  # default
#   MODE=plastic ./dump_time_sweep_study.sh
#   DUMP_TIMES="0.2 0.5 0.8" ./dump_time_sweep_study.sh
#   MAX_DUMP_TIME=3.6 ./dump_time_sweep_study.sh                # 9 evenly-spaced dumps up to 3.6
#   MAX_DUMP_TIME=3.6 N_DUMP_TIMES=15 ./dump_time_sweep_study.sh# 15 evenly-spaced dumps up to 3.6
#   NU=0.49 ./dump_time_sweep_study.sh                          # change fixed nu
#   END_TIME_FOR_RAMP=2.0 ./dump_time_sweep_study.sh            # gentler loading rate
#   N_RESTART_STEPS=20 ./dump_time_sweep_study.sh
#   DT=0.00625 ./dump_time_sweep_study.sh                       # halve the timestep
#   PULSE_AMPLITUDE_BEND=0.05 ./dump_time_sweep_study.sh        # stronger bending pulse
#   PULSE_AMPLITUDE_COMPRESS=0.05 ./dump_time_sweep_study.sh    # stronger compression pulse (more locking)
#   PULSE_AMPLITUDE_BEND=0 PULSE_AMPLITUDE_COMPRESS=0 ./dump_time_sweep_study.sh  # no pulse
#   PULSE_OFFSET_STEPS=4 ./dump_time_sweep_study.sh             # pulse_center = dump_time - 4*dt
#   NY=8 ./dump_time_sweep_study.sh                             # 1x8x1 mesh so a vertical
#                                                               # wave can travel through
#                                                               # interior y-nodes
#
# Output:
#   ${OUT_DIR}/dump_time_sweep_results_new.csv
#   ${OUT_DIR}/dump_time_sweep_results_old.csv

set -euo pipefail

cd "$(dirname "$0")"

RACCOON="${RACCOON:-/home/det12/projects/raccoon/raccoon-opt}"
MODE="${MODE:-elastic_qs}"
NU="${NU:-0.4}"
END_TIME_FOR_RAMP="${END_TIME_FOR_RAMP:-1.0}"

# Three ways to set the dump-time grid (in priority order):
#
#  1) DUMP_TIMES="0.2 0.5 0.8 ..."   explicit list, wins over everything
#  2) MAX_DUMP_TIME=2.5               generates N_DUMP_TIMES evenly-spaced
#                                     values from MAX_DUMP_TIME/N to MAX_DUMP_TIME
#                                     (so the spacing scales with the max while
#                                     the count stays fixed); N_DUMP_TIMES
#                                     defaults to 9 to match the historical
#                                     default sweep length
#  3) Neither set                    fall back to the default explicit list
#                                     "0.2 0.4 ... 1.8" (9 values, max 1.8)
#
# Option (2) is the typical way to push the sweep into a different deformation
# regime without changing how many simulations you run.
N_DUMP_TIMES="${N_DUMP_TIMES:-9}"
if [[ -n "${DUMP_TIMES:-}" ]]; then
  read -r -a DUMP_TIMES <<< "$DUMP_TIMES"
elif [[ -n "${MAX_DUMP_TIME:-}" ]]; then
  read -r -a DUMP_TIMES <<< "$(python3 -c "
import sys
N = ${N_DUMP_TIMES}
maxv = ${MAX_DUMP_TIME}
step = maxv / N
print(' '.join(f'{step*(i+1):.6g}' for i in range(N)))
")"
else
  read -r -a DUMP_TIMES <<< "0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8"
fi

NP="${NP:-1}"
N_RESTART_STEPS="${N_RESTART_STEPS:-10}"

# Gaussian BC pulse parameters.  The pulse is split into two components -- a
# bending pulse (amplitude scales with X0 and adds to bending kinematics) and a
# compression pulse (uniform across the top face, drives F-bar locking work by
# forcing the top up uniformly while the bottom is pinned).  Set either env var
# to 0 to disable that component.  Defaults match the .i defaults so an
# unspecified env var is a no-op forwarding.
PULSE_AMPLITUDE_BEND="${PULSE_AMPLITUDE_BEND:-0.01}"
PULSE_AMPLITUDE_COMPRESS="${PULSE_AMPLITUDE_COMPRESS:-0.01}"
PULSE_CENTER="${PULSE_CENTER:-0.05}"
PULSE_WIDTH="${PULSE_WIDTH:-0.025}"

# Optional: pin pulse_center relative to dump_time instead of using an absolute
# value.  If PULSE_OFFSET_STEPS is set, pulse_center is computed per-dump_time
# as `dump_time - PULSE_OFFSET_STEPS * dt`, so the pulse always fires the same
# number of timesteps before the dump regardless of when the dump occurs.
# That keeps the bouncing-wave state at dump_time at the same "wave age" across
# every point of the sweep, isolating dump_time as the only variable.  Leave
# PULSE_OFFSET_STEPS unset to use the absolute PULSE_CENTER above.
PULSE_OFFSET_STEPS="${PULSE_OFFSET_STEPS:-}"

# Reference-mesh resolution (overrides Mesh/gmg/n{x,y,z} on the reference run
# only -- the restart uses FileMeshGenerator and picks up whatever the
# reference dumped).  Defaults reproduce the single-element behaviour; set
# NY > 1 to thicken in the bending direction so a vertical stress wave can
# propagate through interior nodes instead of being absorbed by the BC-pinned
# top and bottom of the single element.
NX="${NX:-1}"
NY="${NY:-1}"
NZ="${NZ:-1}"

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

# Default OUT_DIR includes the mesh dimensions when they differ from the
# single-element default so multi-element sweeps don't clobber the 1x1x1 runs
# (and so it's easy to tell at a glance what mesh produced a given results
# folder).  CLI override still wins.
if [[ "$NX" == "1" && "$NY" == "1" && "$NZ" == "1" ]]; then
  OUT_DIR="${OUT_DIR:-outputs_dump_time_sweep_${MODE}}"
else
  OUT_DIR="${OUT_DIR:-outputs_dump_time_sweep_${MODE}_${NX}x${NY}x${NZ}}"
fi

# MPI launcher autodetect.
if [[ -n "${MPIRUN:-}" ]]; then
  :
elif command -v mpirun >/dev/null 2>&1; then
  MPIRUN=mpirun
else
  MPIRUN="conda run -n moose mpirun"
fi

# dt is forwarded to the MOOSE runs and is also used locally for the
# target_time / first_time calculation so they stay in sync.  If DT is set in
# the environment, that wins; otherwise auto-extract it from the reference's
# input file as a sensible default.
if [[ -n "${DT:-}" ]]; then
  DT_SOURCE="env override"
else
  DT=$(awk -F'=' '/^[[:space:]]*dt[[:space:]]*=/{gsub(/[[:space:]]/,"",$2); print $2; exit}' "$REF_INPUT")
  DT_SOURCE="auto-extracted from $REF_INPUT"
fi
DT="${DT:-0.0125}"

echo "Single-HEX8 dump_time-sweep config:"
echo "  raccoon-opt        = $RACCOON"
echo "  MPIRUN             = $MPIRUN"
echo "  MODE               = $MODE"
echo "  inputs             = $REF_INPUT + $RES_INPUT"
echo "  OUT_DIR            = $OUT_DIR"
echo "  NP                 = $NP"
echo "  NU (fixed)         = $NU"
echo "  END_TIME_FOR_RAMP  = $END_TIME_FOR_RAMP"
echo "  MAX_DUMP_TIME      = ${MAX_DUMP_TIME:-<unset>}   N_DUMP_TIMES = $N_DUMP_TIMES"
echo "  DUMP_TIMES         = ${DUMP_TIMES[*]}"
echo "  dt                 = $DT  ($DT_SOURCE)"
echo "  N_RESTART_STEPS    = $N_RESTART_STEPS"
echo "  PULSE_AMPLITUDE_BEND     = $PULSE_AMPLITUDE_BEND"
echo "  PULSE_AMPLITUDE_COMPRESS = $PULSE_AMPLITUDE_COMPRESS"
echo "  PULSE_CENTER       = $PULSE_CENTER"
echo "  PULSE_WIDTH        = $PULSE_WIDTH"
echo "  PULSE_OFFSET_STEPS = ${PULSE_OFFSET_STEPS:-<unset>}"
echo "  NX x NY x NZ       = ${NX} x ${NY} x ${NZ}"

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

# Output layout:
#   $OUT_DIR/
#     new/                          one subfolder per recovery method
#       dump_time_sweep_results.csv   the per-method summary CSV
#       csv/                          per-pair timeseries CSVs, named with dump_tag
#         pull_test_out_dt_XpY.csv
#         pull_test_restart_out_dt_XpY.csv
#       exodus/                       per-pair exodus subfolders
#         dt_XpY/                       reference + restart .e files for one dump
#           pull_test_out.e
#           pull_test_out_disp.e
#           pull_test_restart_out.e
#     old/                           same layout for the OLD method
mkdir -p "${OUT_DIR}/new" "${OUT_DIR}/old"
RESULTS_NEW="${OUT_DIR}/new/dump_time_sweep_results.csv"
RESULTS_OLD="${OUT_DIR}/old/dump_time_sweep_results.csv"
RESULTS_HEADER="dump_time,target_time,fbar_correction_int_at_dump,fbar_correction_max_at_dump,fbar_pressure_correction_int_at_dump,fbar_pressure_correction_max_at_dump,psie_ref_first,psie_restart_first,abs_err_first,rel_err_first,psie_ref_target,psie_restart_target,abs_err_target,rel_err_target"
echo "$RESULTS_HEADER" > "$RESULTS_NEW"
echo "$RESULTS_HEADER" > "$RESULTS_OLD"

# Plastic modes need be_bar carried through both legs so the isochoric plastic
# state survives the dump/recover roundtrip.  Elastic modes don't have a be_bar
# material so adding it would error.
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

# run_method <dump_time> <dump_tag> <method> <recover_apply_fbar_to_U> <restart_tensor_mats> <ref_recover_tensor_mats>
run_method() {
  local dump_time="$1"
  local dump_tag="$2"
  local method="$3"
  local apply_fbar="$4"
  local restart_tensor_mats="$5"
  local ref_recover_tensor_mats="$6"
  # Final destination layout (see header comment).  dump_tag is like "_dt_0p2";
  # strip the leading underscore for the exodus subfolder name.
  local dump_tag_bare="${dump_tag#_}"
  local method_dir="$OUT_DIR/$method"
  local csv_dir="$method_dir/csv"
  local exodus_pair_dir="$method_dir/exodus/$dump_tag_bare"
  local workdir="$method_dir/_work_$dump_tag_bare"
  mkdir -p "$csv_dir" "$exodus_pair_dir" "$workdir"

  local ref_csv="$csv_dir/pull_test_out${dump_tag}.csv"
  local res_csv="$csv_dir/pull_test_restart_out${dump_tag}.csv"

  local target_time first_time effective_pulse_center
  target_time=$(python3 -c "print(${dump_time} + ${N_RESTART_STEPS}*${DT})")
  first_time=$(python3 -c "print(${dump_time} + ${DT})")

  # Resolve pulse_center: either absolute (PULSE_CENTER) or relative to
  # dump_time (PULSE_OFFSET_STEPS).
  if [[ -n "$PULSE_OFFSET_STEPS" ]]; then
    effective_pulse_center=$(python3 -c "print(${dump_time} - ${PULSE_OFFSET_STEPS} * ${DT})")
  else
    effective_pulse_center=$PULSE_CENTER
  fi

  echo "[dump=$dump_time method=$method] $REF_INPUT  (workdir=$workdir, end_time=$target_time, end_time_for_ramp=$END_TIME_FOR_RAMP)"
  local ref_cmd=(
    $MPIRUN -n "$NP" "$RACCOON" -i "$REF_INPUT"
    out_dir="$workdir" tag= nu=$NU
    dt=$DT
    end_time=$target_time dump_time=$dump_time end_time_for_ramp=$END_TIME_FOR_RAMP
    pulse_amplitude_bend=$PULSE_AMPLITUDE_BEND pulse_amplitude_compress=$PULSE_AMPLITUDE_COMPRESS pulse_center=$effective_pulse_center pulse_width=$PULSE_WIDTH
    Mesh/gmg/nx=$NX Mesh/gmg/ny=$NY Mesh/gmg/nz=$NZ
    "RecoverVariables/rec/tensor_materials=$ref_recover_tensor_mats"
  )
  printf '  $ '; printf '%q ' "${ref_cmd[@]}"; printf '\n'
  "${ref_cmd[@]}" | tail -n 4

  if [[ ! -f "$workdir/pull_test_out.csv" ]]; then
    echo "ERROR: $workdir/pull_test_out.csv not produced -- reference failed at dump_time=$dump_time method=$method" >&2
    exit 1
  fi

  echo "[dump=$dump_time method=$method] $RES_INPUT  (workdir=$workdir, recover_apply_fbar_to_U=$apply_fbar)"
  local res_cmd=(
    $MPIRUN -n "$NP" "$RACCOON" -i "$RES_INPUT"
    out_dir="$workdir" tag= output_tag= nu=$NU
    dt=$DT
    start_time=$dump_time end_time=$target_time end_time_for_ramp=$END_TIME_FOR_RAMP
    pulse_amplitude_bend=$PULSE_AMPLITUDE_BEND pulse_amplitude_compress=$PULSE_AMPLITUDE_COMPRESS pulse_center=$effective_pulse_center pulse_width=$PULSE_WIDTH
    recover_timestep=LATEST
    Materials/defgrad/recover_apply_fbar_to_U=$apply_fbar
    "UserObjects/epsol/tensor_materials=$restart_tensor_mats" -w
  )
  printf '  $ '; printf '%q ' "${res_cmd[@]}"; printf '\n'
  "${res_cmd[@]}" | tail -n 4

  if [[ ! -f "$workdir/pull_test_restart_out.csv" ]]; then
    echo "ERROR: $workdir/pull_test_restart_out.csv not produced -- restart failed at dump_time=$dump_time method=$method" >&2
    exit 1
  fi

  # Move outputs out of the per-pair workdir into the structured layout.
  # CSVs go to csv/ with the dump_tag in the filename so they don't collide
  # across the sweep.  Exodus files go to exodus/$dump_tag/ so the (ref, ref-
  # displaced-dump, restart) trio sits in its own subfolder.
  mv "$workdir/pull_test_out.csv"            "$ref_csv"
  mv "$workdir/pull_test_restart_out.csv"    "$res_csv"
  mv "$workdir/pull_test_out.e"              "$exodus_pair_dir/pull_test_out.e"
  mv "$workdir/pull_test_out_disp.e"         "$exodus_pair_dir/pull_test_out_disp.e"
  mv "$workdir/pull_test_restart_out.e"      "$exodus_pair_dir/pull_test_restart_out.e"
  # Catch any extra files (per-step exodus, checkpoints, etc.) before rmdir.
  shopt -s nullglob
  for stray in "$workdir"/*; do
    mv "$stray" "$exodus_pair_dir/"
  done
  shopt -u nullglob
  rmdir "$workdir"

  # F-bar diagnostics sampled AT the dump time -- this is the x-axis of the
  # sweep (how much volumetric correction work was present at the moment we
  # decided to dump).
  local fbar_int fbar_max fbar_p_int fbar_p_max
  fbar_int=$(extract_col "$ref_csv" "fbar_correction_int" "$dump_time")
  fbar_max=$(extract_col "$ref_csv" "fbar_correction_max" "$dump_time")
  fbar_p_int=$(extract_col "$ref_csv" "fbar_pressure_correction_int" "$dump_time")
  fbar_p_max=$(extract_col "$ref_csv" "fbar_pressure_correction_max" "$dump_time")

  local psie_ref_first psie_res_first psie_ref_targ psie_res_targ
  psie_ref_first=$(extract_col "$ref_csv" "psie_active_int" "$first_time")
  psie_res_first=$(extract_col "$res_csv" "psie_active_int" "$first_time")
  psie_ref_targ=$(extract_col "$ref_csv" "psie_active_int" "$target_time")
  psie_res_targ=$(extract_col "$res_csv" "psie_active_int" "$target_time")

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

  printf '    fbar_pressure_correction_max_at_dump = %s\n' "$fbar_p_max"
  printf '    first-step (t=%s):  abs_err=%s  rel_err=%s\n' "$first_time" "$abs_first" "$rel_first"
  printf '    target     (t=%s):  abs_err=%s  rel_err=%s\n' "$target_time" "$abs_targ" "$rel_targ"

  local results_file
  if [[ "$method" == "new" ]]; then
    results_file="$RESULTS_NEW"
  else
    results_file="$RESULTS_OLD"
  fi
  echo "$dump_time,$target_time,$fbar_int,$fbar_max,$fbar_p_int,$fbar_p_max,$psie_ref_first,$psie_res_first,$abs_first,$rel_first,$psie_ref_targ,$psie_res_targ,$abs_targ,$rel_targ" >> "$results_file"
}

for dump_time in "${DUMP_TIMES[@]}"; do
  echo
  echo "================================================================"
  echo "  dump_time = $dump_time   (NU=$NU, end_time_for_ramp=$END_TIME_FOR_RAMP, NP=$NP)"
  echo "================================================================"

  # tag: prefix _dt_ so the per-pair output files don't collide with nu sweep
  # outputs that may live in the same directory tree.
  DUMP_TAG="_dt_${dump_time//./p}"

  run_method "$dump_time" "$DUMP_TAG" "new" "false" \
    "stress ${EXTRA_TENSOR_MATS} stretch_tensor_fbar rotation_tensor" \
    "stress ${EXTRA_TENSOR_MATS_REF} rotation_tensor stretch_tensor_fbar"

  run_method "$dump_time" "$DUMP_TAG" "old" "true" \
    "stress ${EXTRA_TENSOR_MATS} stretch_tensor rotation_tensor" \
    "stress ${EXTRA_TENSOR_MATS_REF} rotation_tensor stretch_tensor"
done

echo
echo "================================================================"
echo "  dump_time sweep summary -- NEW method (approach A)"
echo "================================================================"
column -ts, "$RESULTS_NEW"
echo
echo "================================================================"
echo "  dump_time sweep summary -- OLD method (approach B)"
echo "================================================================"
column -ts, "$RESULTS_OLD"
echo
echo "Saved: $RESULTS_NEW"
echo "       $RESULTS_OLD"
