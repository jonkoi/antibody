#!/usr/bin/env bash
# Batch CDRH3 redesign against HER2: design -> refold-prepare -> fold -> score.
#
# Generalises the upstream 4-design example (pipeline/run-her2.upstream.sh) into a
# parameterised batch driver, and adds a fold-model sweep so designs can be scored
# by more than one predictor.
#
#   bash pipeline/run_design.sh all
#   N=100 TOPN=20 SEED=7 bash pipeline/run_design.sh all
#   WEIGHTS=her2-2026-11-07 TSRK=1 N=200 bash pipeline/run_design.sh design
#   FOLDER=protenix-v2 bash pipeline/run_design.sh predict score
#
# Run from the folder holding Data/ and Settings/ (on the VM: ~/rfd4).
set -euo pipefail

Q=quay.io/apheris/oslo-hackathon

# --- knobs ---------------------------------------------------------------
WEIGHTS=${WEIGHTS:-her2-novel-2026-11-09}  # her2-novel | her2-2026-11-07 | antibody-2026-10-30
N=${N:-4}                                   # designs to generate (max 200)
TOPN=${TOPN:-2}                             # designs carried into refolding
SEED=${SEED:-1}
TSRK=${TSRK:-}                              # sequence randomness: 1 exploratory, 2 balanced
RECYCLE=${RECYCLE:-high/front_loaded}       # recommended for antibodies
INPUT=${INPUT:-Data/12-her2/her2_cdrh3.yaml}
FOLDER=${FOLDER:-boltz2}                    # scoring model name
OUT=${OUT:-runs/her2-${WEIGHTS}-n${N}-seed${SEED}}

# map scoring model name -> image tag; template-capable models fold HER2 in context
case "$FOLDER" in
  boltz2)            FOLD_IMAGE=apheris-boltz-v1.3.0-boltz2;           TEMPLATED=1 ;;
  openfold3-of3ob)   FOLD_IMAGE=apheris-openfold3-v1.3.0-of3ob;        TEMPLATED=1 ;;
  openfold3-of3p2)   FOLD_IMAGE=apheris-openfold3-v1.3.0-of3p2;        TEMPLATED=1 ;;
  protenix-v2)       FOLD_IMAGE=apheris-protenix-v1.3.0-v2;            TEMPLATED=1 ;;
  mol1-base)         FOLD_IMAGE=apheris-mol1-v1.3.0-base;              TEMPLATED=1 ;;
  esmfold2-release)  FOLD_IMAGE=apheris-esmfold2-v1.4.0-release;       TEMPLATED=0 ;;
  *) echo "unknown FOLDER=$FOLDER" >&2; exit 2 ;;
esac

RFD=$Q:apheris-rfd4-proteina-v0.3.0-$WEIGHTS
FOLD=$Q:$FOLD_IMAGE

run=(docker run --rm --gpus all --shm-size 16g --user "$(id -u):$(id -g)"
     -e HOME=/tmp -v /etc/passwd:/etc/passwd:ro -v "$PWD:/work" -w /work)

design() {
  echo ">> design: $N designs, weights=$WEIGHTS seed=$SEED -> $OUT/design"
  rm -rf "$OUT/design"
  local extra=()
  [ -n "$TSRK" ] && extra+=(--sequence-tsr-k "$TSRK")
  "${run[@]}" "$RFD" module-run design \
    --input "$INPUT" --output "$OUT/design" \
    --num-replicates "$N" --seed "$SEED" \
    --recycle-strategy "$RECYCLE" "${extra[@]}"
}

prepare() {
  echo ">> refold-prepare: top $TOPN -> $OUT/prepared"
  rm -rf "$OUT/prepared"
  "${run[@]}" "$RFD" module-run refold-prepare \
    --input "$OUT/design" --output "$OUT/prepared" --top-n "$TOPN"
}

predict() {
  echo ">> fold with $FOLDER -> $OUT/$FOLDER"
  rm -rf "$OUT/$FOLDER"
  local cmd template=()
  cmd=$(docker inspect -f '{{join .Config.Cmd " "}}' "$FOLD")
  [ "$TEMPLATED" = 1 ] && template=(--template-paths "$OUT/prepared/holo/template.cif")
  # shellcheck disable=SC2086
  "${run[@]}" "$FOLD" $cmd predict \
    --input "$OUT/prepared/holo/queries.json" --output "$OUT/$FOLDER" \
    "${template[@]}" --write-detailed-scores true \
    --num-diffusion-samples 1 --seeds 42 --query-failures skip
}

score() {
  echo ">> refold-score ($FOLDER) -> $OUT/scored-$FOLDER"
  rm -rf "$OUT/scored-$FOLDER"
  "${run[@]}" "$RFD" module-run refold-score \
    --input "$OUT/design" --output "$OUT/scored-$FOLDER" \
    --predictions "$OUT/$FOLDER" --model "$FOLDER"
  column -s, -t < "$OUT/scored-$FOLDER/refold_metrics.csv" | cut -c1-200
}

[ $# -gt 0 ] || { sed -n '2,14p' "$0"; exit 1; }
for stage in "$@"; do
  case "$stage" in
    design|prepare|predict|score) "$stage" ;;
    all) design; prepare; predict; score ;;
    *) echo "unknown stage: $stage" >&2; exit 2 ;;
  esac
done
echo ">> done: $OUT"
