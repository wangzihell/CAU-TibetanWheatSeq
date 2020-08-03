#!/usr/bin/env bash
set -euxo pipefail

declare -i numb_MATCH
numb_MATCH=0
#
for F in *.o; do
    if cat ${F} | grep -q "Finish"; then
    numb_MATCH+=1
    fi
done
echo "Parts done:"${numb_MATCH}