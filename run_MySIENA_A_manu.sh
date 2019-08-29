#!/bin/bash
#
#PJM --rsc-list "rscgrp=small"
#PJM --rsc-list "node=64"
#PJM --rsc-list "elapse=02:00:00"
#PJM --stg-transfiles all
# --mpi "use-rankdir"
#PJM --stgin "./MySIENA_A_manu ./"
#PJM --stgin-dir "./data/manu ./"
# --stgout "rank=* ./sim_test.txt /data/hp160259/arata/sim_result/sim_test_0_%r.txt"
#PJM -s
#
. /work/system/Env_base
export PARALLEL=8
export OMP_NUM_THREADS=8
mpiexec ./MySIENA_A_manu
