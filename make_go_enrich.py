import csv

big_bash = open("/cbcb/lab/smount/ZCL/bioconductor_scripts/go_enrich/lcm_big_bash.sh", 'w')

# annotations = ["hd_zh_sub_lrpkm_power_8_modulesize_90_TOM_0_merge_0", "hd_zh_sub_lrpkm_power_1_modulesize_40_TOM_0_merge_0", "hd_zh_sub_lrpkm_power_16_modulesize_120_TOM_0_merge_0"]
dir_file = open("/cbcb/project-scratch/ZCL/wgcna/exon_lcm.txt", 'r')
annotations = [l.strip("\n") for l in dir_file]
print(annotations)

count = 100000

for run in annotations:

	dataset = None
	if run[0] == "h":
		dataset = "hd"
	else:
		dataset = "lcm"

	norm = None
	splits = run.split("_")
	if "lrpkm" in splits:
		norm = "lrpkm"
	elif "rld" in splits:
		norm = "rld"
	elif "lcpm" in splits:
		norm = "lcpm"
	elif "log" in splits:
		norm = "log"

	clust_dir = "/cbcb/project-scratch/ZCL/wgcna/"

	##check if merged modules or not
	if run.endswith('0') == True:
		merged = 0
	elif run.endswith('1') == True:
		merged = 1

	print(str(merged))

	lines = []
	#read in cluster assignemnts
	if merged == 1:

		with open(clust_dir + run + "/merged_clusters.csv") as csvfile:
			reader = csv.reader(csvfile, delimiter=',')
			for row in reader:
				if row[2] != "group":
					lines.append(row[2])


	else:

		with open(clust_dir + run + "/original_clusters.csv") as csvfile:
			reader = csv.reader(csvfile, delimiter=',')
			for row in reader:
				if row[2] != "group":
					lines.append(row[2])

	for mod in set(lines):

		if mod != "grey":

			bash = open("/cbcb/lab/smount/ZCL/bioconductor_scripts/go_enrich/" + str(count) + "_bash.sh", 'w')
			big_bash.write("qsub " + str(count) + "_bash.sh\n")

			bash.write("#PBS -q high_throughput\n#PBS -l mem=8GB,walltime=6:00:00,ncpus=2\n")
			bash.write("/cbcb/lab/smount/programs/R-3.1.2/bin/Rscript /cbcb/lab/smount/ZCL/bioconductor_scripts/go_enrich/param_go_mapping.R " + mod + " " + run + " " + dataset + " " +norm + " " + str(merged) + "\n")

			bash.close()

			count += 1

big_bash.close()


		#run the enrichment for each cluster 
	#for mod in  :
