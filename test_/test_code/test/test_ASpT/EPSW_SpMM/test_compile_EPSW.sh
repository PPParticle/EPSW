#!/bin/bash
source ~/.bashrc
riscv-2.6.1 SpMM_EPSW_SP_xt_m2.cc -o SpMM_EPSW_SP_xt_m2 -O3 -static -march=rv64imafdcv0p7_zfh_xtheadc -flax-vector-conversions -mabi=lp64d -mtune=c910
riscv-2.6.1 SpMM_EPSW_SP_xt_m8.cc -o SpMM_EPSW_SP_xt_m8 -O3 -static -march=rv64imafdcv0p7_zfh_xtheadc -flax-vector-conversions -mabi=lp64d -mtune=c910
scp SpMM_EPSW_SP_xt_m8 SpMM_EPSW_SP_xt_m2 root@192.168.1.105:/mnt/qfc/spmm

riscv-2.0 SpMM_EPSW_SP_r_m2.cc -o SpMM_EPSW_SP_r_m2 -O3 -static -march=rv64imafdcvxtheadc -flax-vector-conversions -mabi=lp64d -mtune=c906 
riscv-2.0 SpMM_EPSW_SP_r_m8.cc -o SpMM_EPSW_SP_r_m8 -O3 -static -march=rv64imafdcvxtheadc -flax-vector-conversions -mabi=lp64d -mtune=c906 
riscv-2.0 SpMM_EPSW_SP_r_m8_0.cc -o SpMM_EPSW_SP_r_m8_0 -O3 -static -march=rv64imafdcvxtheadc -flax-vector-conversions -mabi=lp64d -mtune=c906 
riscv-2.0 SpMM_EPSW_SP_r_m8_1.cc -o SpMM_EPSW_SP_r_m8_1 -O3 -static -march=rv64imafdcvxtheadc -flax-vector-conversions -mabi=lp64d -mtune=c906 
scp SpMM_EPSW_SP_r_m8 SpMM_EPSW_SP_r_m2 root@192.168.1.103:/opt/EPSW/spmm


