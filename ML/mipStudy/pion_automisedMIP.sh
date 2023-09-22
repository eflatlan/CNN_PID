#!/bin/bash
set -e  # Stop on error

# Loop 50 times
for i in {0..1}; do
    echo "Iteration: $i"

    # Run your sim_challenge.sh script
    /root/alice/O2/prodtests/sim_challengeMIPStudy.sh -f reco -n 4000 || { echo "sim_challenge.sh failed"; exit 1; }

    # Run ROOT framework and load macro, then call function and exit
    root -b -q -l  MIPstudy.C # << 'EOF' || { echo "ROOT commands failed"; exit 1; }
    #.q
#EOF

    # Rename the output file
    mv MLOUTPUT.root MLOUTPUTPion.root || { echo "File rename failed"; exit 1; }
    #mv ParticleInfo2.h5 ParticleY_211_Highconst_2cm_${i}_.h5 || { echo "File rename failed"; exit 1; }
done
#   taskwrapper hitstats.log root -q -b -l ${O2_ROOT}/share/macro/analyzeHits.C
