#!/bin/bash
export PATH=$PATH:$PWD/../lisa/
export PATH=$PATH:$PWD/../rna-seq/

# dependants are gflags and divsufsort
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/ijerkovic/gflags/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/ijerkovic/divsufsort/lib
