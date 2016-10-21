

import sys

exp_name = sys.argv[1]
adjmat_list = sys.argv[2]

big_bash = open("/cbcb/project-scratch/ZCL/wgcna/consensus/" + exp_name + "/indicator.sh", 'w')

adjmat = open("/cbcb/project-scratch/ZCL/wgcna/consensus/" + exp_name + "/" + adjmat_list)


for l in adjmat:

	num = l.strip("\n").split(".")[0]

	big_bash.write("qsub indmat_bash/"+ num + ".sh\n")

	bash = open("/cbcb/project-scratch/ZCL/wgcna/consensus/" + exp_name + "/indmat_bash/" + num + ".sh", 'w')

	bash.write("#PBS -q throughput\n#PBS -l mem=36GB,walltime=18:00:00,ncpus=2\n")

	bash.write("/cbcb/lab/smount/programs/R-3.1.2/bin/Rscript /cbcb/lab/smount/ZCL/bioconductor_scripts/paral_indicator_mat.R " + exp_name + " " + num)

	bash.close()


big_bash.close()
adjmat.close()