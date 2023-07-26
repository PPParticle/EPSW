#!/bin/bash

echo "dataset, MKL_GFLOPs(K=32), MKL_GFLOPs(K=128), TACO_GFLOPs(K=32), TACO_GFLOPs(K=128), CSB_GFLOPs(K=32), CSB_GFLOPs(K=128), ASpT_GFLOPs(K=32), ASpT_GFLOPs(K=128), ASpT_diff_%(K=32+K=128)" >> SpMM_KNL_SP.out
echo "dataset, preprocessing_ratio" >> SpMM_KNL_SP_preprocessing.out

echo -n 144 >> SpMM_KNL_SP.out
echo -n "," >> SpMM_KNL_SP.out
echo -n 144 >> SpMM_KNL_SP_preprocessing.out
echo -n "," >> SpMM_KNL_SP_preprocessing.out

SpMM_ASpT_SP 144.mtx
echo >> SpMM_KNL_SP.out
echo >> SpMM_KNL_SP_preprocessing.out




