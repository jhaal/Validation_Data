#!/bin/sh
# Worker: claims decks via mkdir under $CLAIM, runs each with rdragon (start lock avoids /tmp/rundir collisions).
CLAIM=/tmp/vrf_claims; LOCK=/tmp/vrf_rdragon_start.lock; mkdir -p $CLAIM
cd /Users/alokjha/DRAGON/5.1/Dragon
for c in $(cat "$1"); do
  mkdir $CLAIM/$c 2>/dev/null || continue
  until mkdir $LOCK 2>/dev/null; do sleep 1; done
  ./rdragon $c.x2m > /dev/null 2>&1 &
  pid=$!; sleep 4; rmdir $LOCK; wait $pid
  echo "$(date +%H:%M:%S) $c branches=$(grep -a -c '^>|VRF_BRANCH_B' Darwin_arm64/$c.result) done=$(grep -a -c '^>|VRF_COMPLETED' Darwin_arm64/$c.result)"
done
