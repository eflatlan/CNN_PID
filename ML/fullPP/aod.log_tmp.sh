#!/usr/bin/env bash
export LIBC_FATAL_STDERR_=1
o2-aod-producer-workflow -b --run --shm-segment-size 10000000000 --aod-writer-keep dangling --aod-writer-resfile AO2D --aod-writer-resmode UPDATE --aod-timeframe-id 1 --run-number 300000;
RC=$?; echo "TASK-EXIT-CODE: ${RC}"; exit ${RC}
