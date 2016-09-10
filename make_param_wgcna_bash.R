
powers = [1,2,4,8,12,16]
minModuleSize = [40, 60, 90, 120, 150, 180, 210]
TOM_on = [0]
files = [ ("hd_zh_sub_lcpm.csv", "hd_zh"), ("hd_zh_sub_log.csv", "hd_zh"), ("hd_zh_sub_lrpkm.csv", "hd_zh") , ("hd_zh_sub_rld.csv", "hd_zh") , ("lcm_sub_lcpm.csv", "lcm") , ("lcm_sub_log.csv" ,"lcm") , ("lcm_sub_lrpkm.csv", "lcm") , ("lcm_sub_rld.csv", "lcm")]

merge_eigengene = [0,1]

file_dir = "/cbcb/lab/smount/ZCL/bioconductor_scripts/counts/subsets/"
dir = "/cbcb/lab/smount/ZCL/bioconductor_scripts/wgcna/"

count = 1 

big_bash = open(dir+"big_bash.sh", 'w')

for p in powers:
	for m in minModuleSize:
		for t in TOM_on:
			for e in merge_eigengene:
				for tup in files:

					bash = open(dir + "bash"+str(count)+".sh", 'w')
					big_bash.write("qsub "+"bash"+str(count)+".sh\n")

					bash.write("#PBS -q throughput\n#PBS -l mem=36GB,walltime=12:00:00,ncpus=4\n")
					bash.write("/cbcb/lab/smount/programs/R-3.1.2/bin/Rscript /cbcb/lab/smount/ZCL/bioconductor_scripts/wgcna/param_wgcna.R" + " " + str(p) + " " + str(m) + " " + str(t) + " " + str(e) + " " + file_dir + tup[0] + " " + tup[1] + "\n")

					bash.close()
					count += 1

big_bash.close()



