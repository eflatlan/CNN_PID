#!/bin/bash

#set -e  # Stop on error

# Loop 50 times
for i in {0..10250}; do
    echo "Iteration: $i"

    # Run your sim_challenge.sh script
    /root/alice/O2/prodtests/sim_challenge_LEAD.sh -f simp -n 100 -j 5 || { echo "sim_challengeCustom.sh failed"; continue; }

    # Run ROOT framework and load macro, then call function with value 1.75 and exit
    root -b -q -l "SegmentationCkov.C(1.75)" || { echo "ROOT commands failed"; continue; }

    mv ParticleInfo2.h5 Particle211_Ckov${i}_.h5 || { echo "File rename failed"; continue; }
done

# taskwrapper hitstats.log root -q -b -l ${O2_ROOT}/share/macro/analyzeHits.C

