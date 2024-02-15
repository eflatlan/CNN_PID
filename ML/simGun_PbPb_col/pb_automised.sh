#!/bin/bash

for i in {0..100250}; do
    success=0
    while [[ $success -eq 0 ]]; do
        echo "Iteration: $i"

        if /root/alice/O2/prodtests/sim_challenge.sh -f sim -n 1 -j 5 -s pbpb; then
            if root -b -q -l "SegmentationCkov.C(1.75)"; then
                # Check if file size is greater than 1KB
                fileSize=$(stat -c %s ParticleInfo2.h5)
                if [[ $fileSize -gt 1024 ]]; then
                    # Rename and store the file only if it's larger than 1KB
                    if mv ParticleInfo2.h5 ParticlePbPb_Ckov${i}_.h5; then
                        success=1
                    else
                        echo "File rename failed"
                    fi
                else
                    echo "File size is less than or equal to 1KB, not storing the file"   
                    # do not iter
                    #success=1 # Mark success to move on to the next iteration
                fi
            else
                echo "ROOT commands failed"
            fi
        else
            echo "sim_challengeCustom.sh failed"
        fi
    done
done

