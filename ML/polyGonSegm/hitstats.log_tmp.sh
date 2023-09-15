#!/usr/bin/env bash
export LIBC_FATAL_STDERR_=1
root -q -b -l /root/alice/sw/slc8_x86-64/O2/segment-local1/share/macro/analyzeHits.C;
RC=$?; echo "TASK-EXIT-CODE: ${RC}"; exit ${RC}
