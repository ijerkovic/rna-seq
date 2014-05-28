#!/bin/bash
REL_PATH_TO_RNA_SEQ=/home/ijerkovic/rna-seq

PATH=$PATH:$REL_PATH_TO_RNA_SEQ/lisa
PATH=$PATH:$REL_PATH_TO_RNA_SEQ/rna-seq
#PATH=$PATH:$PWD/../lisa/
#PATH=$PATH:$PWD/../rna-seq/
export PATH

# dependants are gflags and divsufsort
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/ijerkovic/gflags/lib
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/ijerkovic/divsufsort/lib
export LB_LIBRARY_PATH

