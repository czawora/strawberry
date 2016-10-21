

import sys

exp_name = sys.argv[1]
lst_name = sys.argv[2]

bash = open("/cbcb/project-scratch/ZCL/wgcna/consensus/" + exp_name + "/indicator.sh", 'w')

bash.write("#PBS -q large\n#PBS -l mem=120GB,walltime=12:00:00,ncpus=10\n")

bash.write("/cbcb/lab/smount/programs/R-3.1.2/bin/Rscript /cbcb/lab/smount/ZCL/bioconductor_scripts/indicator_mat.R " + exp_name + " " + lst_name + " &> /cbcb/project-scratch/ZCL/wgcna/consensus/" + exp_name + "/indicator.log\n")

bash.close()

