#!/bin/bash


for i in {0..100250}; do
    success=0
    while [[ $success -eq 0 ]]; do
        echo "Iteration: $i"


        if /root/alice/O2/prodtests/sim_challenge.sh -f sim -n 1 -j 5 -s pbpb; then

            if root -b -q -l "SegmentationCkov.C(1.75)"; then
                # Rename the file
                if mv ParticleInfo2.h5 ParticlePbPb_Ckov${i}_.h5; then
                    success=1
                else
                    echo "File rename failed"
                fi
            else
                echo "ROOT commands failed"
            fi
        else
            echo "sim_challengeCustom.sh failed"
        fi
    done
done
