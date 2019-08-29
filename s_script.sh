#!/bin/bash

mpiFCCpx -c -Xg -std=c++11 -Nstl=libc++ -Kfast,parallel,openmp,optmsg=2 -Nlst=t -o MySIENA.o ./script/MySIENA.cpp
mpiFCCpx -c -Xg -std=c++11 -Nstl=libc++ -Kfast,parallel,openmp,optmsg=2 -Nlst=t -o constant.o ./script/constant.cpp
# mpiFCCpx -Xg -std=c++11 -Nstl=libc++ -Kfast,parallel,openmp,optmsg=2 -Nlst=t -o MySIENA_A_manu constant.o MySIENA.o -SSL2
# mpiFCCpx -Xg -std=c++11 -Nstl=libc++ -Kfast,parallel,openmp,optmsg=2 -Nlst=t -o MySIENA_B_manu constant.o MySIENA.o -SSL2
# mpiFCCpx -Xg -std=c++11 -Nstl=libc++ -Kfast,parallel,openmp,optmsg=2 -Nlst=t -o MySIENA_C_manu constant.o MySIENA.o -SSL2
# mpiFCCpx -Xg -std=c++11 -Nstl=libc++ -Kfast,parallel,openmp,optmsg=2 -Nlst=t -o MySIENA_A_all constant.o MySIENA.o -SSL2
mpiFCCpx -Xg -std=c++11 -Nstl=libc++ -Kfast,parallel,openmp,optmsg=2 -Nlst=t -o MySIENA_B_all constant.o MySIENA.o -SSL2
# mpiFCCpx -Xg -std=c++11 -Nstl=libc++ -Kfast,parallel,openmp,optmsg=2 -Nlst=t -o MySIENA_C_all constant.o MySIENA.o -SSL2
# sed -i 's/\r//' run_MySIENA_A_manu.sh
# sed -i 's/\r//' run_MySIENA_B_manu.sh
# sed -i 's/\r//' run_MySIENA_C_manu.sh
# sed -i 's/\r//' run_MySIENA_A_all.sh
sed -i 's/\r//' run_MySIENA_B_all.sh
# sed -i 's/\r//' run_MySIENA_C_all.sh
# pjsub run_MySIENA_A_manu.sh
# pjsub run_MySIENA_B_manu.sh
# pjsub run_MySIENA_C_manu.sh
# pjsub run_MySIENA_A_all.sh
pjsub run_MySIENA_B_all.sh
# pjsub run_MySIENA_C_all.sh

# for simulations gof_1_flag == true
# mpiFCCpx -c -Xg -std=c++11 -Nstl=libc++ -Kfast,optmsg=2 -Nlst=t -o MySIENA.o ./script/MySIENA.cpp
# mpiFCCpx -c -Xg -std=c++11 -Nstl=libc++ -Kfast,optmsg=2 -Nlst=t -o constant.o ./script/constant.cpp
# mpiFCCpx -Xg -std=c++11 -Nstl=libc++ -Kfast,optmsg=2 -Nlst=t -o MySIENA_sim_C_all constant.o MySIENA.o -SSL2
# sed -i 's/\r//' run_MySIENA_sim_C_all.sh
# pjsub run_MySIENA_sim_C_all.sh

# for gof_2_flag == true
# mpiFCCpx -c -Xg -std=c++11 -Nstl=libc++ -Kfast,optmsg=2 -Nlst=t -o MySIENA.o ./script/MySIENA.cpp
# mpiFCCpx -c -Xg -std=c++11 -Nstl=libc++ -Kfast,optmsg=2 -Nlst=t -o constant.o ./script/constant.cpp
# mpiFCCpx -Xg -std=c++11 -Nstl=libc++ -Kfast,optmsg=2 -Nlst=t -o MySIENA_sim_entire_all constant.o MySIENA.o -SSL2
# sed -i 's/\r//' run_MySIENA_sim_entire_all.sh
# pjsub run_MySIENA_sim_entire_all.sh
