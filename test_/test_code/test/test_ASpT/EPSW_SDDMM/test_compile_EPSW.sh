#!/bin/bash
source ~/.bashrc
riscv-2.6.1 SDDMM_EPSW_SP_xt_row_m2.cc -o SDDMM_EPSW_SP_xt_m2 -O2 -static -march=rv64imafdcv0p7_zfh_xtheadc -flax-vector-conversions -mabi=lp64d -mtune=c910
riscv-2.6.1 SDDMM_EPSW_SP_xt_row_m8.cc -o SDDMM_EPSW_SP_xt_m8 -O2 -static -march=rv64imafdcv0p7_zfh_xtheadc -flax-vector-conversions -mabi=lp64d -mtune=c910
scp SDDMM_EPSW_SP_xt_m8 SDDMM_EPSW_SP_xt_m2 root@192.168.1.112:/mnt/qfc/sddmm

riscv-2.0 SDDMM_EPSW_SP_r_row_m2.cc -o SDDMM_EPSW_SP_r_m2 -O2 -static -march=rv64imafdcvxtheadc -flax-vector-conversions -mabi=lp64d -mtune=c906 
riscv-2.0 SDDMM_EPSW_SP_r_row_m8.cc -o SDDMM_EPSW_SP_r_m8 -O2 -static -march=rv64imafdcvxtheadc -flax-vector-conversions -mabi=lp64d -mtune=c906 
scp SDDMM_EPSW_SP_r_m8 SDDMM_EPSW_SP_r_m2 root@192.168.1.102:/opt/sddmm


