#!/bin/sh

# Sickle cell AS
# Normal AA

# python ../Sickle-trait-RNAseq/TAKT/python/rna_rna_diss.py -t1  -t2  -c 10000 -s 15 -o 'Results'

# Questions
# are RNA delayed between traits of each parasite: AS 3D7 vs AA 3D7, AS FUP vs AA FUP.

ECHO AS 3D7 vs AA 3D7
python ../Sickle-trait-RNAseq/TAKT/python/rna_rna_diss.py -t1 AS_3D7_1stage_periodic_jtk05.txt -t2 AA_3D7_1stage_periodic_jtk05.txt -c 10000 -s 15 -o 'Results'
ECHO ------------------------------------

ECHO AS FUP vs AA FUP
python ../Sickle-trait-RNAseq/TAKT/python/rna_rna_diss.py -t1 AS_FUP_1stage_periodic_jtk05.txt -t2 AA_FUP_1stage_periodic_jtk05.txt -c 10000 -s 15 -o 'Results'
ECHO ------------------------------------


# are RNA delayed between parasites in AS: AS 3D7 vs AS FUP
ECHO AS 3D7 vs AS FUP
python ../Sickle-trait-RNAseq/TAKT/python/rna_rna_diss.py -t1 AS_3D7_1stage_periodic_jtk05.txt -t2 AS_FUP_1stage_periodic_jtk05.txt -c 10000 -s 15 -o 'Results'
ECHO ------------------------------------

# are RNA delayed between parasites in AA: AA 3D7 vs AA FUP
ECHO AA 3D7 vs AA FUP
python ../Sickle-trait-RNAseq/TAKT/python/rna_rna_diss.py -t1 AA_3D7_1stage_periodic_jtk05.txt -t2 AA_FUP_1stage_periodic_jtk05.txt -c 10000 -s 15 -o 'Results'
ECHO ------------------------------------
