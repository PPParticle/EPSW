#!/bin/bash
# source ~/.bashrc
riscv test_vf_xt_32f.cc -o test_vf_xt_32f -O2 -static -march=rv64imafdcv0p7_zfh_xtheadc -mabi=lp64d -mtune=c910  
riscv test_vi_xt.cc -o test_vi_xt -O2 -static -march=rv64imafdcv0p7_zfh_xtheadc -mabi=lp64d -mtune=c910  
riscv test_vl_xt.cc -o test_vl_xt -O2 -static -march=rv64imafdcv0p7_zfh_xtheadc -mabi=lp64d -mtune=c910  

nezha test_vf_r_32f.cc -o test_vf_r_32f -O2 -static -march=rv64imafdcvxtheadc -mabi=lp64d -mtune=c906 
nezha test_vi_r.cc -o test_vi_r -O2 -static -march=rv64imafdcvxtheadc -mabi=lp64d -mtune=c906 
nezha test_vl_r.cc -o test_vl_r -O2 -static -march=rv64imafdcvxtheadc -mabi=lp64d -mtune=c906 