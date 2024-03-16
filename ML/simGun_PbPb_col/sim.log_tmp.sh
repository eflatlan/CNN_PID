#!/usr/bin/env bash
export LIBC_FATAL_STDERR_=1
o2-sim -n4 --configKeyValues Diamond.width[2]=6. -g pythia8hi -e TGeant3 -j 4 --run 300000;
RC=$?; echo "TASK-EXIT-CODE: ${RC}"; exit ${RC}
