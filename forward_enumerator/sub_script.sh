#!/bin/bash

#PBS -N forward_enumeration
#PBS -o /home/hwkim/simple_reaction_prediction/template_application/GitHub/out.out
#PBS -l nodes=cnode13:ppn=16
#PBS -l walltime=7:00:00:00
#PBS -e /home/hwkim/simple_reaction_prediction/template_application/GitHub/error.txt


##### Run ##### 
#NODE=$(uniq < $PBS_NODEFILE)

date

conda activate HWKim

cd /home/hwkim/simple_reaction_prediction/template_application/GitHub/forward_enumerator/

python main_function.py --config main_config.yaml

date

