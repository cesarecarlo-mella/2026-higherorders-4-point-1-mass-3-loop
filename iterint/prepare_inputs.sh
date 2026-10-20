#!/bin/bash
# Copies the Atilde of every family into iterint/input/atilde/<fam>_atilde.m (from
# all.ints/a_tildes, <fam>_atilde.m or <fam>_Atilde.m), so that a clone of the repo has every
# input of the cluster run. Only needed if an Atilde changes; the copies are committed.
# The BC are read directly from families/results/<fam>_bc.m.
set -e
here="$(cd "$(dirname "$0")" && pwd)"
ALLINTS="/Users/cesarecarlomella/Desktop/Physics/PhD/XjetProject/INTEGRALS/4P1M/3L1M/all.ints"
mkdir -p "$here/input/atilde"
while read -r f; do
  [ -z "$f" ] && continue
  src=$(ls "$ALLINTS/a_tildes/${f}_atilde.m" "$ALLINTS/a_tildes/${f}_Atilde.m" 2>/dev/null | head -1)
  [ -z "$src" ] && { echo "no atilde for $f"; exit 1; }
  cp "$src" "$here/input/atilde/${f}_atilde.m"
done < "$here/families.txt"
echo "atilde for $(grep -c . "$here/families.txt") families in $here/input/atilde"
