# QPET integration for SZ3/QoZ/HPEZ (pointwise and regional QoI)

## Introduction

This branch includes the codes for SZ3-QPET, QoZ-QPET, and HPEZ-QPET on pointwise and regional QoIs (single-snapshot compression).

For vector QoI preservation on cross-snapshot compression, check the codes at: [link] 

## Dependencies

Please Install the following dependencies before compiling HPEZ:

* cmake>=3.13
* gcc>=6.0
* SymEngine (https://github.com/symengine/symengine)
* GMP (https://gmplib.org/) (The dependency of SymEngine)
* Zstd >= 1.3.5 (https://facebook.github.io/zstd/). It is not mandatory to manually install as Zstandard v1.4.5 is included and will be used if libzstd can not be found by pkg-config.

## Installation

* mkdir build && cd build
* cmake -DCMAKE_INSTALL_PREFIX:PATH=[INSTALL_DIR] ..
* make
* make install

Then, you'll find all the executables in [INSTALL_DIR]/bin and header files in [INSTALL_DIR]/include. A Cmake version >= 3.13.0 is needed. 
Before you proceed to the following evaluations, you may add the installation path of HPEZ to your system path so that you can directly run the 'hpez' command in your machine for further evaluation. 
Otherwise, regard the 'hpez' command name in the following as '[INSTALL_DIR]/bin/hpez' (the path of your installed executable).

## HPEZ (base) Compression/Decompression Examples

This section is irrelevant to QoI-preserving compression but provides fundamental instructions for how to run several error-bounded lossy compressors in the SZ family. 

T和 'hpez' command integrates 5 different compression levels by the argument -q, corresponding to 3 compressors:

* hpez -q 0: SZ3.1 compression.
* hpez -q 1: QoZ 1.1 compression.
* hpez -q 2/3/4: different optimization levels of HPEZ compression (level 3 recommended, which is the default).

Notice: the integrated SZ3.1 and QoZ 1.1 in HPEZ (QoZ 2.0) have already leveraged the Fast-varying-first interpolation (proposed in our paper), therefore their compression ratios are sometimes higher than the original public released versions of SZ3.1 and QoZ 1.1.

## Test Dataset

Please download test datasets from: https://sdrbench.github.io/. 

