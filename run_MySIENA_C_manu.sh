#!/bin/bash
#
#PJM --rsc-list "rscgrp=small"
#PJM --rsc-list "node=64"
#PJM --rsc-list "elapse=02:00:00"
#PJM --stg-transfiles all
#PJM --stgin "./MySIENA_C_manu ./"
#PJM --stgin-dir "./data/manu ./"
#PJM -s
#
. /work/system/Env_base
export PARALLEL=8
export OMP_NUM_THREADS=8
mpiexec ./MySIENA_C_manu
