#!/bin/bash
set -e  # Stop on error

# Loop 50 times
for i in {1..49}; do
    echo "Iteration: $i"

    # Run your sim_challenge.sh script
    /root/alice/O2/prodtests/sim_challenge.sh -f simp -n 3000 || { echo "sim_challenge.sh failed"; exit 1; }

    # Run ROOT framework and load macro, then call function and exit
    root -b -q -l  Segmentation.C # << 'EOF' || { echo "ROOT commands failed"; exit 1; }
    #.q
#EOF

    # Rename the output file
    mv ParticleInfo2.h5 Particle211_Highconst_2cm_${i}_.h5 || { echo "File rename failed"; exit 1; }
done
#   taskwrapper hitstats.log root -q -b -l ${O2_ROOT}/share/macro/analyzeHits.C
