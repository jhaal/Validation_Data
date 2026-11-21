#!/bin/sh
# Run all VRF_FULL_<lib>_* decks for one library, serializing DRAGON start-up with a lock
# so two queues never pick the same /tmp/rundirN.
lib=$1; LOCK=/tmp/vrf_rdragon_start.lock; cd /Users/alokjha/DRAGON/5.1/Dragon
for f in data/VRF_${PFX:-FULL}_${lib}_BU*.x2m; do
  c=$(basename $f .x2m)
  until mkdir $LOCK 2>/dev/null; do sleep 1; done
  ./rdragon $c.x2m > /dev/null 2>&1 &
  pid=$!; sleep 4; rmdir $LOCK; wait $pid
  echo "$(date +%H:%M:%S) $c branches=$(grep -a -c '^>|VRF_BRANCH_B' Darwin_arm64/$c.result) done=$(grep -a -c '^>|VRF_COMPLETED' Darwin_arm64/$c.result)"
done
