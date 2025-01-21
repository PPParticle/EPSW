#!/bin/bash
source ~/.bashrc
riscv-2.6.1 SpMSpM_EPSW_SP_xt.cc -o SpMSpM_EPSW_SP_xt -O3 -static -march=rv64imafdcv0p7_zfh_xtheadc -flax-vector-conversions -mabi=lp64d -mtune=c910
riscv-2.6.1 transform.cc -o transform_xt -O3 -static -march=rv64imafdcv0p7_zfh_xtheadc -mabi=lp64d -mtune=c910
riscv-2.0 SpMSpM_EPSW_SP_r.cc -o SpMSpM_EPSW_SP_r -O3 -static -march=rv64imafdcvxtheadc -flax-vector-conversions -mabi=lp64d -mtune=c906 
riscv-2.0 transform.cc -o transform_r -O3 -static -march=rv64imafdcvxtheadc -mabi=lp64d -mtune=c906 
riscv-2.6.1 SpMSpM_EPSW_SP_xt.cc -o SpMSpM_EPSW_SP_xt -O3 -static -march=rv64imafdcvxtheadc -flax-vector-conversions -mabi=lp64d -mtune=c908
riscv-2.6.1 transform.cc -o transform_xt -O3 -static -march=rv64imafdcvxtheadc -mabi=lp64d -mtune=c908

