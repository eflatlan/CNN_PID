#!/bin/bash

set -e  # Stop on error

# Loop 50 times
for i in {0..100}; do
    echo "Iteration: $i"

    # Run your sim_challenge.sh script
    /root/alice/O2/prodtests/sim_challengeCustom.sh -f simpr -n 150 -j 5 || { echo "sim_challengeCustom.sh failed"; exit 1; }

    # Run ROOT framework and load macro, then call function with value 1.75 and exit
    root -b -q -l "SegmentationCkov.C(1.75)" || { echo "ROOT commands failed"; exit 1; }

    mv ParticleInfo2.h5 Particle2212_Ckov${i}_.h5 || { echo "File rename failed"; exit 1; }
done

# taskwrapper hitstats.log root -q -b -l ${O2_ROOT}/share/macro/analyzeHits.C

