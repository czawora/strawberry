
import sys

fname = sys.argv[1]
fdir = sys.argv[2]
outdir = sys.argv[3]

f = open(fname)
big_bash = open(outdir + "/big_bash.sh", 'w')

count = 0
for line in f:

	bash = open(outdir+ "/bash" + str(count) + ".sh", 'w')
	big_bash.write("qsub bash" + str(count) + ".sh\n")

	line = line.strip("\n")
	

	bash.write("#PBS -q throughput\n#PBS -l mem=36GB,ncpus=2\n")
	bash.write("/cbcb/lab/smount/programs/R-3.1.2/bin/Rscript /cbcb/lab/smount/ZCL/bioconductor_scripts/wgcna_consensus/list2mat.R " + fdir + "/"+ line + "\n")

	count += 1
	bash.close()

big_bash.close()
f.close()
