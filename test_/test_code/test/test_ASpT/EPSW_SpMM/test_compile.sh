#!/bin/bash
source ~/.bashrc
riscv-2.6.1 SpMM_ASpT_SP_xt.cc -o SpMM_ASpT_SP_xt -O3 -static -march=rv64imafdcv0p7_zfh_xtheadc -flax-vector-conversions -mabi=lp64d -mtune=c910 -Wconversion-null
riscv-2.6.1 SpMM_CSR_SP_xt_m2.cc -o SpMM_CSR_SP_xt_m2 -O3 -static -march=rv64imafdcv0p7_zfh_xtheadc -flax-vector-conversions -mabi=lp64d -mtune=c910 -Wconversion-null
riscv-2.6.1 SpMM_CSR_SP_xt_m8.cc -o SpMM_CSR_SP_xt_m8 -O3 -static -march=rv64imafdcv0p7_zfh_xtheadc -flax-vector-conversions -mabi=lp64d -mtune=c910 -Wconversion-null
riscv-2.6.1 SpMM_EPSW_SP_xt_m2.cc -o SpMM_EPSW_SP_xt_m2 -O3 -static -march=rv64imafdcv0p7_zfh_xtheadc -flax-vector-conversions -mabi=lp64d -mtune=c910 -Wconversion-null
riscv-2.6.1 SpMM_EPSW_SP_xt_m8.cc -o SpMM_EPSW_SP_xt_m8 -O3 -static -march=rv64imafdcv0p7_zfh_xtheadc -flax-vector-conversions -mabi=lp64d -mtune=c910 -Wconversion-null

riscv-2.0 SpMM_ASpT_SP_r.cc -o SpMM_ASpT_SP_r -O3 -static -march=rv64imafdcvxtheadc -mabi=lp64d -mtune=c906 -Wconversion-null
riscv-2.0 SpMM_CSR_SP_r_m2.cc -o SpMM_CSR_SP_r_m2 -O3 -static -march=rv64imafdcvxtheadc -mabi=lp64d -mtune=c906 -Wconversion-null
riscv-2.0 SpMM_CSR_SP_r_m8.cc -o SpMM_CSR_SP_r_m8 -O3 -static -march=rv64imafdcvxtheadc -mabi=lp64d -mtune=c906 -Wconversion-null
riscv-2.0 SpMM_EPSW_SP_r_m2.cc -o SpMM_EPSW_SP_r_m2 -O3 -static -march=rv64imafdcvxtheadc -mabi=lp64d -mtune=c906 -Wconversion-null
riscv-2.0 SpMM_EPSW_SP_r_m8.cc -o SpMM_EPSW_SP_r_m8 -O3 -static -march=rv64imafdcvxtheadc -mabi=lp64d -mtune=c906 -Wconversion-null

riscv-1.0 SpMM_ASpT_SP_xt.cc -o SpMM_ASpT_SP_xt -O3 -static -march=rv64imafdcv_zihintpause_zfh_zba_zbb_zbc_zbs_xtheadc -mabi=lp64d -mtune=c908 -Wconversion-null -flax-vector-conversions
riscv-2.1 SpMM_CSR_SP_xt_m2.cc -o SpMM_CSR_SP_xt_m2 -O3 -static -march=rv64imafdcv0p7_zfh_xthead -mabi=lp64d -mtune=c908 -Wconversion-null -flax-vector-conversions
riscv-2.1 SpMM_CSR_SP_xt_m8.cc -o SpMM_CSR_SP_xt_m8 -O3 -static -march=rv64imafdcv0p7_zfh_xthead -mabi=lp64d -mtune=c908 -Wconversion-null -flax-vector-conversions
riscv-2.1 SpMM_EPSW_SP_xt_m2.cc -o SpMM_EPSW_SP_xt_m2 -O3 -static -march=rv64imafdcv0p7_zfh_xthead -mabi=lp64d -mtune=c908 -Wconversion-null -flax-vector-conversions
riscv-2.1 SpMM_EPSW_SP_xt_m8.cc -o SpMM_EPSW_SP_xt_m8 -O3 -static -march=rv64imafdcv0p7_zfh_xthead -mabi=lp64d -mtune=c908 -Wconversion-null -flax-vector-conversions