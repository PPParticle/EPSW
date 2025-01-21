#!/bin/bash
source ~/.bashrc
riscv-2.6.1 SpMV_EPSW_SP_xt_m8.cc -o SpMV_EPSW_SP_xt_m8 -O2 -static -march=rv64imafdcv0p7_zfh_xtheadc -flax-vector-conversions -mabi=lp64d -mtune=c910
riscv-2.6.1 SpMV_EPSW_SP_xt_m2.cc -o SpMV_EPSW_SP_xt_m2 -O2 -static -march=rv64imafdcv0p7_zfh_xtheadc -flax-vector-conversions -mabi=lp64d -mtune=c910
riscv-2.6.1 SpMV_CSR_SP_xt_m8.cc -o SpMV_CSR_SP_xt_m8 -O2 -static -march=rv64imafdcv0p7_zfh_xtheadc -flax-vector-conversions -mabi=lp64d -mtune=c910
riscv-2.6.1 SpMV_CSR_SP_xt_m2.cc -o SpMV_CSR_SP_xt_m2 -O2 -static -march=rv64imafdcv0p7_zfh_xtheadc -flax-vector-conversions -mabi=lp64d -mtune=c910
riscv-2.6.1 SpMV_ASpT_SP_xt.cc -o SpMV_ASpT_SP_xt -O2 -static -march=rv64imafdcv0p7_zfh_xtheadc -flax-vector-conversions -mabi=lp64d -mtune=c910

riscv-2.0 SpMV_EPSW_SP_r_m8.cc -o SpMV_EPSW_SP_r_m8 -O2 -static -march=rv64imafdcvxtheadc -flax-vector-conversions -mabi=lp64d -mtune=c906 
riscv-2.0 SpMV_EPSW_SP_r_m2.cc -o SpMV_EPSW_SP_r_m2 -O2 -static -march=rv64imafdcvxtheadc -flax-vector-conversions -mabi=lp64d -mtune=c906 
riscv-2.0 SpMV_CSR_SP_r_m8.cc -o SpMV_CSR_SP_r_m8 -O2 -static -march=rv64imafdcvxtheadc -flax-vector-conversions -mabi=lp64d -mtune=c906 
riscv-2.0 SpMV_CSR_SP_r_m2.cc -o SpMV_CSR_SP_r_m2 -O2 -static -march=rv64imafdcvxtheadc -flax-vector-conversions -mabi=lp64d -mtune=c906 
riscv-2.0 SpMV_ASpT_SP_r.cc -o SpMV_ASpT_SP_r -O2 -static -march=rv64imafdcvxtheadc -flax-vector-conversions -mabi=lp64d -mtune=c906 