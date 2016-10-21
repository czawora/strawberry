

import sys
import os

exp = sys.argv[1]

try:
	os.mkdir("/cbcb/project-scratch/ZCL/wgcna/consensus/" + exp + "/cluster_bash")
except:
	pass

clist = open("/cbcb/project-scratch/ZCL/wgcna/consensus/" + exp + "/cluster_list.txt")

out = open("/cbcb/project-scratch/ZCL/wgcna/consensus/" + exp + "/clust_to_mat.sh", 'w')

bcount = 0

for line in clist:

	bash = open("/cbcb/project-scratch/ZCL/wgcna/consensus/" + exp + "/cluster_bash/" + str(bcount) + ".sh", 'w')

	out.write("qsub " + str(bcount) + ".sh\n")

	bash.write("#PBS -q throughput\n#PBS -l mem=36GB,walltime=18:00:00\n")

	bash.write("/cbcb/lab/smount/programs/R-3.1.2/bin/Rscript /cbcb/lab/smount/ZCL/bioconductor_scripts/paral_con_mat.R " + exp + " " + line.split(".")[0] + "\n")

	bash.close()

	bcount += 1

out.close()
clist.close()