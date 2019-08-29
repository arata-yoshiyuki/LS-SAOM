#!/bin/bash
#
#PJM --rsc-list "rscgrp=small"
#PJM --rsc-list "node=128"
#PJM --mpi "proc=1024"
#PJM --rsc-list "elapse=01:00:00"
#PJM --stg-transfiles all
#PJM --mpi "use-rankdir"
#PJM --stgin "./MySIENA_sim_C_all ./"
#PJM --stgin-dir "./data/all_sec ./"
#PJM --stgout "rank=* ./sim_test.txt /data/hp160259/arata/sim_result_large/sim_test_C_all_%r.txt"
#PJM -s
#
. /work/system/Env_base
mpiexec ./MySIENA_sim_C_all
