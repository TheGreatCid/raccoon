#!/usr/bin/env bash
# Mesh-convergence study for recover/remesh initialization error (F-bar variant).
#
# For each n in N_VALUES:
#   1. run reference.i with end_time = ref_end_time + dt (ONE extra step past
#      ref_end_time so its CSV has a row at target_time = ref_end_time + dt).
#      reference_fbar_out_disp_${n}.e is written at every TIMESTEP_END so the
#      file contains a slice at ref_end_time as well as the extra slice.
#   2. run restart.i  with start_time = ref_end_time, end_time = target_time
#      AND recover_timestep = (slice-index-of-ref_end_time) so the restart
#      recovers from the SAME physical state the reference had at ref_end_time
#      and then advances by one dt.
#
# Both CSVs now have a TIMESTEP_END row at target_time, so
#   abs_err = |psie_restart(target_time) - psie_reference(target_time)|
# is the recovery error after one forward step.  On the same mesh with perfect
# recovery this is zero; otherwise it scales with mesh refinement.
#
# Why not compare at t = ref_end_time directly?  At EXEC_INITIAL (or
# TIMESTEP_BEGIN) MOOSE calls Material::initStatefulProperties but NOT the
# regular Material::computeQpProperties, so material-property integral
# postprocessors are zero on the recovered state at that flag.  TIMESTEP_END
# is the first stage at which the recovered materials are evaluated -- hence
# the extra-step trick.
#
# Usage:
#   ./mesh_convergence.sh                              # N_VALUES = NP_VALUES = 1 2 4 8 16
#   N_VALUES="2 4 8" NP_VALUES="2 4 8" ./mesh_convergence.sh
#   OUT_DIR=runs RACCOON=/path/to/raccoon-opt ./mesh_convergence.sh
#
# Cores-per-sim scaling:
#   NP_VALUES is a positional list of mpirun rank counts; element i is the
#   rank count used for N_VALUES[i].  Defaults to the same list as N_VALUES,
#   so each n uses n ranks (n=1 -> 1 rank, n=2 -> 2, ..., n=16 -> 16).  If
#   NP_VALUES is shorter than N_VALUES, missing entries fall back to using
#   n itself as the rank count.
#
# Output:
#   ${OUT_DIR}/convergence_results.csv
#       columns: n, np, n_elem, psie_ref_last, psie_restart_init, abs_err, rel_err
#   ${OUT_DIR}/*.e, *.csv  (all per-n exodus and csv files)
#   stdout: progress + final aligned table

set -euo pipefail

cd "$(dirname "$0")"

RACCOON="${RACCOON:-/home/det12/projects/raccoon/raccoon-opt}"
OUT_DIR="${OUT_DIR:-outputs}"
read -r -a N_VALUES  <<< "${N_VALUES:-1 2 4 8 16}"
read -r -a NP_VALUES <<< "${NP_VALUES:-1 2 4 8 16}"

# Number of timesteps the restart leg runs PAST ref_end_time.  The reference
# is extended by the same number so both have TIMESTEP_END rows at every
# intermediate time, and at the final comparison time.  Higher = more forward
# drift folded into the error metric; lower = closer to pure initialization
# error.  Defaults to 5.
N_RESTART_STEPS="${N_RESTART_STEPS:-5}"

# MPI launcher.  If mpirun isn't on PATH, fall back to `conda run -n moose
# mpirun` since the moose conda env ships its own libmpi.  Override with
# MPIRUN=... to force a specific launcher.
if [[ -n "${MPIRUN:-}" ]]; then
  : # honor override as-is
elif command -v mpirun >/dev/null 2>&1; then
  MPIRUN=mpirun
else
  MPIRUN="conda run -n moose mpirun"
fi

# Pull end_time and dt out of reference.i so we always extract the right row,
# even if someone re-tunes them upstream.  awk picks the FIRST top-level
# `end_time = ...` and `dt = ...` assignments in the file.
REF_END_TIME=$(awk -F'=' '/^[[:space:]]*end_time[[:space:]]*=/{gsub(/[[:space:]]/,"",$2); print $2; exit}' reference.i)
REF_END_TIME="${REF_END_TIME:-2.5}"
DT=$(awk -F'=' '/^[[:space:]]*dt[[:space:]]*=/{gsub(/[[:space:]]/,"",$2); print $2; exit}' reference.i)
DT="${DT:-0.1}"

# target_time = ref_end_time + N_RESTART_STEPS * dt.  Reference is extended
# by N_RESTART_STEPS steps past ref_end_time so it has a TIMESTEP_END row at
# target_time (and at every intermediate restart step) for direct
# row-by-row comparison.
TARGET_TIME=$(python3 -c "print(${REF_END_TIME} + ${N_RESTART_STEPS}*${DT})")

# With reference.i's exodusqp gated to sync_times = dump_time + sync_only,
# the displaced recovery dump contains exactly one slice (at dump_time, which
# this script sets to REF_END_TIME).  recover_timestep = LATEST therefore
# always picks the right one.
RECOVER_TIMESTEP="LATEST"

echo "Convergence study config:"
echo "  raccoon-opt      = $RACCOON"
echo "  MPIRUN           = $MPIRUN"
echo "  OUT_DIR          = $OUT_DIR"
echo "  N_VALUES         = ${N_VALUES[*]}"
echo "  NP_VALUES        = ${NP_VALUES[*]}"
echo "  ref_end_time     = $REF_END_TIME"
echo "  dt               = $DT"
echo "  N_RESTART_STEPS  = $N_RESTART_STEPS  (steps after ref_end_time)"
echo "  target_time      = $TARGET_TIME   (= ref_end_time + N_RESTART_STEPS*dt)"
echo "  dump_time        = $REF_END_TIME  (displaced recovery dump time)"
echo "  recover_timestep = $RECOVER_TIMESTEP  (single slice in dump)"

if [[ ! -x "$RACCOON" ]]; then
  echo "ERROR: raccoon-opt not found or not executable at $RACCOON" >&2
  exit 1
fi

# Make sure the output subdirectory exists.  MOOSE will write into it via the
# file_base = ${out_dir}/... pattern in the input files.
mkdir -p "$OUT_DIR"

# ---------------------------------------------------------------------------
# extract_psie <csv_file> <target_time>
#   Locates the column named "psie_corr_active_int" and the row whose "time"
#   column is closest to target_time, then prints that value.
# ---------------------------------------------------------------------------
extract_psie() {
  local csvfile="$1"
  local target="$2"
  if [[ ! -f "$csvfile" ]]; then
    echo "NaN"
    return 1
  fi
  awk -F, -v target="$target" '
    NR==1 {
      for (i=1; i<=NF; i++) {
        h=$i; gsub(/[[:space:]\r"]/,"",h)
        col[h] = i
      }
      next
    }
    {
      t = $col["time"] + 0.0
      v = $col["psie_corr_active_int"] + 0.0
      diff = (t > target) ? t - target : target - t
      if (NR==2 || diff < best_diff) {
        best_diff = diff
        best_v    = v
      }
    }
    END { printf "%.16e\n", best_v }
  ' "$csvfile"
}

RESULTS="${OUT_DIR}/convergence_results.csv"
# Two comparison points are emitted:
#   first_step   = ref_end_time + dt        (initialization-error indicator,
#                                            minimal forward drift)
#   target_time  = ref_end_time + N*dt      (long-horizon drift indicator)
echo "n,np,n_elem,psie_ref_first,psie_restart_first,abs_err_first,rel_err_first,psie_ref_target,psie_restart_target,abs_err_target,rel_err_target" > "$RESULTS"

FIRST_STEP_TIME=$(python3 -c "print(${REF_END_TIME} + ${DT})")

for i in "${!N_VALUES[@]}"; do
  n="${N_VALUES[$i]}"
  # Default rank count = matching NP_VALUES entry; if NP_VALUES is shorter
  # than N_VALUES, fall back to using n itself as the rank count.
  np="${NP_VALUES[$i]:-$n}"

  echo
  echo "================================================================"
  echo "  n = $n   (np = $np)"
  echo "================================================================"

  # Each hex of an nxnxn grid is split into 24 TET10 by libMesh, so
  # n_elem = 24 * n^3.  Reported just for context.
  n_elem=$(( 24 * n * n * n ))

  REF_CSV="${OUT_DIR}/reference_fbar_out_${n}.csv"
  RES_CSV="${OUT_DIR}/restart_fbar_out_${n}.csv"

  echo "[$n] reference.i (np=$np, end_time=$TARGET_TIME, dump_time=$REF_END_TIME) ..."
  # Extend the reference past ref_end_time so its CSV has TIMESTEP_END rows
  # at target_time AND all intermediate restart-comparison times.  But hold
  # dump_time at ref_end_time so the displaced exodusqp recovery dump still
  # captures the kinematic state at ref_end_time, which is what the restart
  # is meant to recover from.  -w because end_time/dump_time overrides may
  # leave top-level symbols unused.
  $MPIRUN -n "$np" "$RACCOON" -i reference.i n=$n out_dir=$OUT_DIR \
    end_time=$TARGET_TIME dump_time=$REF_END_TIME -w \
    | tail -n 20

  if [[ ! -f "$REF_CSV" ]]; then
    echo "ERROR: $REF_CSV not produced — reference run failed" >&2
    exit 1
  fi

  echo "[$n] restart.i (np=$np, recover_timestep=$RECOVER_TIMESTEP) ..."
  # recover_timestep points the SolutionUserObjectQP at the exodus slice for
  # t = ref_end_time, NOT the LATEST slice (which is at target_time).
  $MPIRUN -n "$np" "$RACCOON" -i restart.i n=$n out_dir=$OUT_DIR \
    start_time=$REF_END_TIME end_time=$TARGET_TIME \
    recover_timestep=$RECOVER_TIMESTEP -w \
    | tail -n 20

  if [[ ! -f "$RES_CSV" ]]; then
    echo "ERROR: $RES_CSV not produced — restart run failed" >&2
    exit 1
  fi

  # First-step values (initialization-error indicator) and target_time
  # values (after N_RESTART_STEPS forward steps).
  psie_ref_first=$(extract_psie "$REF_CSV"  "$FIRST_STEP_TIME")
  psie_res_first=$(extract_psie "$RES_CSV"  "$FIRST_STEP_TIME")
  psie_ref_targ=$(extract_psie "$REF_CSV"  "$TARGET_TIME")
  psie_res_targ=$(extract_psie "$RES_CSV"  "$TARGET_TIME")

  read abs_first rel_first abs_targ rel_targ < <(python3 -c "
import sys
rf, rsf = ${psie_ref_first}, ${psie_res_first}
rt, rst = ${psie_ref_targ},  ${psie_res_targ}
def pair(ref, res):
    abs_e = abs(res - ref)
    rel_e = abs_e / max(abs(ref), 1e-30)
    return f'{abs_e:.16e} {rel_e:.16e}'
print(pair(rf, rsf), pair(rt, rst))
")

  printf '  first-step (t=%s):\n' "$FIRST_STEP_TIME"
  printf '    psie_reference  = %s\n' "$psie_ref_first"
  printf '    psie_restart    = %s\n' "$psie_res_first"
  printf '    abs_err         = %s\n' "$abs_first"
  printf '    rel_err         = %s\n' "$rel_first"
  printf '  target (t=%s, after %d restart steps):\n' "$TARGET_TIME" "$N_RESTART_STEPS"
  printf '    psie_reference  = %s\n' "$psie_ref_targ"
  printf '    psie_restart    = %s\n' "$psie_res_targ"
  printf '    abs_err         = %s\n' "$abs_targ"
  printf '    rel_err         = %s\n' "$rel_targ"

  echo "$n,$np,$n_elem,$psie_ref_first,$psie_res_first,$abs_first,$rel_first,$psie_ref_targ,$psie_res_targ,$abs_targ,$rel_targ" >> "$RESULTS"
done

echo
echo "================================================================"
echo "  Convergence summary"
echo "================================================================"
column -ts, "$RESULTS"
echo
echo "Saved: $RESULTS"
echo "All exodus / csv outputs are in: $OUT_DIR/"
