#!/bin/bash
source ~/.bashrc
riscv-2.6.1 SDDMM_EPSW_SP_xt_m2.cc -o SDDMM_EPSW_SP_xt_m2 -O2 -static -march=rv64imafdcv0p7_zfh_xtheadc -flax-vector-conversions -mabi=lp64d -mtune=c910
riscv-2.6.1 SDDMM_EPSW_SP_xt_m8.cc -o SDDMM_EPSW_SP_xt_m8 -O2 -static -march=rv64imafdcv0p7_zfh_xtheadc -flax-vector-conversions -mabi=lp64d -mtune=c910
riscv-2.6.1 SDDMM_CSR_SP_xt_m2.cc -o SDDMM_CSR_SP_xt_m2 -O2 -static -march=rv64imafdcv0p7_zfh_xtheadc -flax-vector-conversions -mabi=lp64d -mtune=c910
riscv-2.6.1 SDDMM_CSR_SP_xt_m8.cc -o SDDMM_CSR_SP_xt_m8 -O2 -static -march=rv64imafdcv0p7_zfh_xtheadc -flax-vector-conversions -mabi=lp64d -mtune=c910
riscv-2.6.1 SDDMM_ASpT_SP_xt.cc -o SDDMM_ASpT_SP_xt -O2 -static -march=rv64imafdcv0p7_zfh_xtheadc -flax-vector-conversions -mabi=lp64d -mtune=c910

riscv-2.0 SDDMM_EPSW_SP_r_m2.cc -o SDDMM_EPSW_SP_r_m2 -O2 -static -march=rv64imafdcvxtheadc -flax-vector-conversions -mabi=lp64d -mtune=c906 
riscv-2.0 SDDMM_EPSW_SP_r_m8.cc -o SDDMM_EPSW_SP_r_m8 -O2 -static -march=rv64imafdcvxtheadc -flax-vector-conversions -mabi=lp64d -mtune=c906 
riscv-2.0 SDDMM_CSR_SP_r_m2.cc -o SDDMM_CSR_SP_r_m2 -O2 -static -march=rv64imafdcvxtheadc -flax-vector-conversions -mabi=lp64d -mtune=c906 
riscv-2.0 SDDMM_CSR_SP_r_m8.cc -o SDDMM_CSR_SP_r_m8 -O2 -static -march=rv64imafdcvxtheadc -flax-vector-conversions -mabi=lp64d -mtune=c906 
riscv-2.0 SDDMM_ASpT_SP_r.cc -o SDDMM_ASpT_SP_r -O2 -static -march=rv64imafdcvxtheadc -flax-vector-conversions -mabi=lp64d -mtune=c906 

riscv-1.0 SDDMM_ASpT_SP_xt.cc -o SDDMM_ASpT_SP_xt -O3 -static -march=rv64imafdcvxtheadc -mabi=lp64d -mtune=c908 -Wconversion-null -flax-vector-conversions
riscv-1.0 SDDMM_CSR_SP_xt_m2.cc -o SDDMM_CSR_SP_xt_m2 -O3 -static -march=rv64imafdcvxtheadc -mabi=lp64d -mtune=c908 -Wconversion-null -flax-vector-conversions
riscv-1.0 SDDMM_CSR_SP_xt_m8.cc -o SDDMM_CSR_SP_xt_m8 -O3 -static -march=rv64imafdcvxtheadc -mabi=lp64d -mtune=c908 -Wconversion-null -flax-vector-conversions
riscv-1.0 SDDMM_EPSW_SP_xt_m2.cc -o SDDMM_EPSW_SP_xt_m2 -O3 -static -march=rv64imafdcvxtheadc -mabi=lp64d -mtune=c908 -Wconversion-null -flax-vector-conversions
riscv-1.0 SDDMM_EPSW_SP_xt_m8.cc -o SDDMM_EPSW_SP_xt_m8 -O3 -static -march=rv64imafdcvxtheadc -mabi=lp64d -mtune=c908 -Wconversion-null -flax-vector-conversions