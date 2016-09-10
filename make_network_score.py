
networks = open("/cbcb/project-scratch/ZCL/go_enrich_results/exon_hd.txt")

big_bash = open("/cbcb/lab/smount/ZCL/bioconductor_scripts/net_scores/hd_big_bash.sh", 'w')

count = 1

for net in networks:

	net = net.strip("\n")

	bash = open("/cbcb/lab/smount/ZCL/bioconductor_scripts/net_scores/" + str(count) + "_bash.sh" ,'w')

	bash.write("#PBS -q high_throughput\n#PBS -l mem=2GB,walltime=6:00:00\n")

	bash.write("/cbcb/lab/smount/programs/R-3.1.2/bin/Rscript /cbcb/lab/smount/ZCL/bioconductor_scripts/net_scores/param_network_score.R " +net +  "\n")

	big_bash.write("qsub " + str(count) + "_bash.sh\n")

	count += 1

	bash.close()

big_bash.close()
networks.close()



