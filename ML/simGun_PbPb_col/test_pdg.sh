#!/bin/bash
if /root/alice/O2/prodtests/sim_challengeCust.sh -f simk -n 400 -j 4; then
    if root -b -q -l "SegmentationCkov2.C(1.75)"; then
        # Check if file size is greater than 1KB
        echo "run ok"
    else
        echo "ROOT commands failed"
    fi
else
    echo "sim_challengeCustom.sh failed"
fi

