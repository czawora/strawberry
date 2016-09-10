
import csv

f = open("/cbcb/lab/smount/ZCL/bioconductor_scripts/params.csv" , 'r')

big_bash = open("/cbcb/lab/smount/ZCL/bioconductor_scripts/wgcna_bash/big_bash.sh", 'w')

reader = csv.reader(f)

count = 1
   	
for row in reader:

	bash = open("/cbcb/lab/smount/ZCL/bioconductor_scripts/wgcna_bash/" +str(count)+ "_bash.sh", 'w')


	power = row[1]
	module = row[2]
	tom = row[3]

	if tom == "1":

		bash.write("#PBS -q large\n#PBS -l mem=120GB,walltime=24:00:00\n")

	else:
		
		bash.write("#PBS -q throughput\n#PBS -l mem=36GB,walltime=12:00:00\n")

	bash.write("/cbcb/lab/smount/programs/R-3.1.2/bin/Rscript " + "/cbcb/lab/smount/ZCL/bioconductor_scripts/rpkm.R " + power + " " + module + " " + tom + " &> /cbcb/lab/smount/ZCL/bioconductor_scripts/wgcna_bash/output/"+str(count)+ "_bash.out\n")

	bash.close()

	big_bash.write("qsub "+  str(count)+ "_bash.sh\n")

	count += 1

big_bash.close()
