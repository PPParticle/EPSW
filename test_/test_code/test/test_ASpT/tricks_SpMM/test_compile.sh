#!/bin/bash
source ~/.bashrc
riscv-2.6.1 SpMM_EPSW_SP_trick_xt.cc -o SpMM_EPSW_SP_trick_xt -O2 -static -march=rv64imafdcv0p7_zfh_xtheadc -flax-vector-conversions -mabi=lp64d -mtune=c910
riscv-2.0 SpMM_EPSW_SP_trick_r.cc -o SpMM_EPSW_SP_trick_r -O2 -static -march=rv64imafdcvxtheadc -flax-vector-conversions -mabi=lp64d -mtune=c906 
scp SpMM_EPSW_SP_trick_r root@192.168.1.102:/opt/tricks
scp SpMM_EPSW_SP_trick_xt root@192.168.1.112:/mnt/qfc/tricks


