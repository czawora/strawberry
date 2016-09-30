


import sys

fname = sys.argv[1]

f = open(fname)
big_bash = open("big_bash.sh", 'w')

count = 0
for line in f:

	bash = open("bash" + str(count) + ".sh", 'w')
	big_bash.write("qsub bash" + str(count) + ".sh\n")

	line = line.strip("\n")
	merge = line.split("_")[-1]

	bash.write("#PBS -q throughput\n#PBS -l mem=36GB,ncpus=2\n")
	bash.write("/cbcb/lab/smount/programs/R-3.1.2/bin/Rscript /cbcb/lab/smount/ZCL/bioconductor_scripts/wgcna_consensus/wgcna_paral_clust_list.R " + line + " " + merge + "\n")

	count += 1
	bash.close()

big_bash.close()
f.close()

