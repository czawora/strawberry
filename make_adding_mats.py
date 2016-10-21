

import sys
import os
import subprocess

exp_name = sys.argv[1]
rnd = sys.argv[2]
adjmat_name = sys.argv[3]
sort = sys.argv[4]

adj_list = open("/cbcb/project-scratch/ZCL/wgcna/consensus/" + exp_name + 
	"/" + adjmat_name )

big_bash = open("/cbcb/project-scratch/ZCL/wgcna/consensus/" + exp_name + 
	"/merge" + rnd + ".sh", 'w')

resulting_mats = open("/cbcb/project-scratch/ZCL/wgcna/consensus/" + exp_name + "/merge" + rnd + "_mats.txt", 'w')

count = 1
bash_count = 0

m1 = None
m2 = None

for l in adj_list:

	l = l.strip("\n")

	if count % 2 == 1:

		m1 = l

	if count % 2 == 0:

		m2 = l

		outname = rnd + "_" + str(bash_count)

		print(outname)

		resulting_mats.write(outname + ".csv\n")

		big_bash.write("qsub merge_bash/" + outname + ".sh\n")

		bash = open("/cbcb/project-scratch/ZCL/wgcna/consensus/" + exp_name + "/merge_bash/" + outname + ".sh", 'w')

		bash.write("#PBS -q throughput\n#PBS -l mem=36GB,walltime=12:00:00,ncpus=2\n")

		if sort == "sort":

			bash.write("/cbcb/lab/smount/programs/R-3.1.2/bin/Rscript /cbcb/lab/smount/ZCL/bioconductor_scripts/adding_mats.R " + exp_name + " " + m1 + " " + m2 + " " + outname + " sort\n")


		else:

			bash.write("/cbcb/lab/smount/programs/R-3.1.2/bin/Rscript /cbcb/lab/smount/ZCL/bioconductor_scripts/adding_mats.R " + exp_name + " " + m1 + " " + m2 + " " + outname + " _\n")

		#print(exp_name + " " + m1 + " " + m2 + " " + outname)

		bash.close()

		bash_count += 1

		m1 = None
		m2 = None

	count += 1

#odd number of mats to merge
#copy the last mat with rnd name
if m1 != None and m2 == None:

	print("copying extra matrix")
	print(m1)
	subprocess.call(["cp", "/cbcb/project-scratch/ZCL/wgcna/consensus/" + exp_name + "/adjmat/" + m1, "/cbcb/project-scratch/ZCL/wgcna/consensus/" + exp_name + "/adjmat/" + rnd + "_" + str(bash_count) + ".csv"])

	resulting_mats.write( rnd + "_" + str(bash_count) + ".csv\n")



big_bash.close()
resulting_mats.close()
adj_list.close()
