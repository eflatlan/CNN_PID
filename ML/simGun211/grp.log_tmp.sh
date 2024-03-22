#!/usr/bin/env bash
export LIBC_FATAL_STDERR_=1
o2-grp-simgrp-tool createGRPs --run 300000 --publishto GRP -o mcGRP;
RC=$?; echo "TASK-EXIT-CODE: ${RC}"; exit ${RC}
