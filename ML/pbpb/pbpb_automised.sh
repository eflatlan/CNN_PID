#!/bin/bash

set -e  # Stop on error

# Loop 50 times
for i in {0..2150}; do
    echo "Iteration: $i"

    # Run your sim_challenge.sh script
    /root/alice/O2/prodtests/sim_challenge.sh -s pbpb -n 1 -j 5 -e TGeant3|| { echo "sim_challenge.sh failed"; exit 1; }

    # Run ROOT framework and load macro, then call function with value 1.75 and exit
    root -b -q -l "SegmentationCkov.C(1.75)" || { echo "ROOT commands failed"; exit 1; }

    mv ParticleInfo2.h5 ParticlePbPb_Ckov${i}_.h5 || { echo "File rename failed"; exit 1; }
done

# taskwrapper hitstats.log root -q -b -l ${O2_ROOT}/share/macro/analyzeHits.C

