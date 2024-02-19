#!/usr/bin/env bash
export LIBC_FATAL_STDERR_=1
o2-sim-digitizer-workflow -b --run --shm-segment-size 10000000000 --disable-trd-trapsim --interactionRate 50000 --configKeyValues HBFUtils.runNumber=300000 --early-forward-policy always --combine-devices;
RC=$?; echo "TASK-EXIT-CODE: ${RC}"; exit ${RC}
