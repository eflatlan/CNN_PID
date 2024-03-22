#!/usr/bin/env bash
export LIBC_FATAL_STDERR_=1
o2-hmpid-matcher-workflow -b --run --shm-segment-size 10000000000;
RC=$?; echo "TASK-EXIT-CODE: ${RC}"; exit ${RC}
