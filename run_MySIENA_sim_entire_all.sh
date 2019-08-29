#!/bin/bash
#
#PJM --rsc-list "rscgrp=small"
#PJM --rsc-list "node=8"
#PJM --mpi "proc=64"
#PJM --rsc-list "elapse=03:00:00"
#PJM --stg-transfiles all
#PJM --mpi "use-rankdir"
#PJM --stgin "./MySIENA_sim_entire_all ./"
#PJM --stgin-dir "./data/all_sec ./"
#PJM --stgout "rank=* ./sim_entire.txt /data/hp160259/arata/sim_result_entire_null/sim_entire_all_%r.txt"
#PJM -s
#
. /work/system/Env_base
mpiexec ./MySIENA_sim_entire_all
