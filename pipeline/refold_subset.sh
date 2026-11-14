#!/usr/bin/env bash
# Refold a chosen SUBSET of designs, instead of refold-prepare's top-N ranking.
#
# Why this exists: refold-prepare ranks on interface H-bonds and hotspot
# obedience, which in our n=100 batch selected almost entirely high-substitution
# designs -- and 30% of what it picked was later rejected by the novelty gate,
# so that GPU time was wasted. Worse, every design it picked sat in the
# substitution range where measured affinity collapses. This lets us fold the
# designs we actually want: novelty-accepted, and conservative enough to bind.
#
#   IDS="her2_cdrh3_4 her2_cdrh3_88" bash refold_subset.sh
#
# Run from the folder holding Data/ and Settings/ (on the VM: ~/rfd4).
set -euo pipefail

Q=quay.io/apheris/oslo-hackathon
WEIGHTS=${WEIGHTS:-her2-novel-2026-11-09}
RUN=${RUN:-runs/team15-claude-n100-seed1}
OUT=${OUT:-$RUN/subset}
FOLDER=${FOLDER:-boltz2}
FOLD_IMAGE=${FOLD_IMAGE:-apheris-boltz-v1.3.0-boltz2}
IDS=${IDS:?set IDS to a space-separated list of query_ids}

RFD=$Q:apheris-rfd4-proteina-v0.3.0-$WEIGHTS
FOLD=$Q:$FOLD_IMAGE
run=(docker run --rm --gpus all --shm-size 16g --user "$(id -u):$(id -g)"
     -e HOME=/tmp -v /etc/passwd:/etc/passwd:ro -v "$PWD:/work" -w /work)

echo ">> refold-prepare over ALL designs (ranking ignored; we filter ourselves)"
rm -rf "$OUT/prepared"
"${run[@]}" "$RFD" module-run refold-prepare \
  --input "$RUN/design" --output "$OUT/prepared" --top-n 100

echo ">> filtering queries to: $IDS"
IDS="$IDS" python3 -I - "$OUT/prepared/holo/queries.json" <<'PY'
import json, os, sys
path = sys.argv[1]
keep = set(os.environ["IDS"].split())
doc = json.load(open(path))
before = len(doc["queries"])
doc["queries"] = [q for q in doc["queries"] if q["query_id"] in keep]
missing = keep - {q["query_id"] for q in doc["queries"]}
json.dump(doc, open(path, "w"), indent=1)
print(f"   kept {len(doc['queries'])} of {before} queries")
if missing:
    print(f"   WARNING not found: {sorted(missing)}", file=sys.stderr)
    if not doc["queries"]:
        sys.exit("   no queries left -- aborting")
PY

echo ">> folding with $FOLDER"
rm -rf "$OUT/$FOLDER"
cmd=$(docker inspect -f '{{join .Config.Cmd " "}}' "$FOLD")
# shellcheck disable=SC2086
"${run[@]}" "$FOLD" $cmd predict \
  --input "$OUT/prepared/holo/queries.json" --output "$OUT/$FOLDER" \
  --template-paths "$OUT/prepared/holo/template.cif" \
  --write-detailed-scores true --num-diffusion-samples 1 --seeds 42 \
  --query-failures skip

echo ">> scoring"
rm -rf "$OUT/scored-$FOLDER"
"${run[@]}" "$RFD" module-run refold-score \
  --input "$RUN/design" --output "$OUT/scored-$FOLDER" \
  --predictions "$OUT/$FOLDER" --model "$FOLDER"

column -s, -t < "$OUT/scored-$FOLDER/refold_metrics.csv" | cut -c1-160
echo ">> done: $OUT"
