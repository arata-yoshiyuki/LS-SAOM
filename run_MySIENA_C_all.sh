#!/bin/bash
#
#PJM --rsc-list "rscgrp=small"
#PJM --rsc-list "node=256"
#PJM --rsc-list "elapse=24:00:00"
#PJM --stg-transfiles all
#PJM --stgin "./MySIENA_C_all ./"
#PJM --stgin-dir "./data/all_sec ./"
#PJM -s
#
. /work/system/Env_base
export PARALLEL=8
export OMP_NUM_THREADS=8
mpiexec ./MySIENA_C_all
