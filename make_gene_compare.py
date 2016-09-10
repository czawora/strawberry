
to_compare = open("/cbcb/project-scratch/ZCL/comparison_mats/lcm_top50.txt", 'r')

to_compare = [l.strip("\n") for l in to_compare]

big_bash = open("/cbcb/project-scratch/ZCL/comparison_mats/gene_comp/big_bash.sh", 'w')
count = 1

seen_list = []

for n1 in to_compare:

	n1 = n1[1:-1]

	for n2 in to_compare:

		n2 = n2[1:-1]

		if n2 not in seen_list:

			bash = open("/cbcb/project-scratch/ZCL/comparison_mats/gene_comp/bash" + str(count) + ".sh", 'w')

			big_bash.write("qsub " + "bash" + str(count) + ".sh\n" )

			bash.write("#PBS -q high_throughput\n#PBS -l mem=8GB,walltime=8:00:00\n")
			bash.write("/cbcb/lab/smount/programs/R-3.1.2/bin/Rscript /cbcb/lab/smount/ZCL/bioconductor_scripts/gene_compare.R " +n1 + " " + n2 + "\n")

			bash.close()

			count += 1

	seen_list.append(n1)

big_bash.close()
