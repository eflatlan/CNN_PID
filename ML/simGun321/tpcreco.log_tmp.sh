#!/usr/bin/env bash
export LIBC_FATAL_STDERR_=1
o2-tpc-reco-workflow -b --run --shm-segment-size 10000000000 --input-type digits --output-type clusters,tracks,send-clusters-per-sector --configKeyValues GPU_rec.maxTrackQPtB5=20;
RC=$?; echo "TASK-EXIT-CODE: ${RC}"; exit ${RC}
