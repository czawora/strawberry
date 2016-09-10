
import datetime
import time
import os

ts = time.time()
st = datetime.datetime.fromtimestamp(ts).strftime('%Y_%m_%d__%H_%M_%S')

network_file = open("/cbcb/project-scratch/ZCL/comparison_mats/compare_list.txt", 'r')
out_dir = "/cbcb/project-scratch/ZCL/comparison_mats/" + st
os.makedirs(out_dir)

networks = [l for l in network_file]

big_bash = open(out_dir + "/big_bash.sh", 'w')

count = 1

seen_list = []

for net in networks:

	net = net.strip("\n")
	print(net)
	print(str(len(seen_list)))

	##check if merged modules or not
	if net.endswith('0') == True:
		n1merged = 0
	elif net.endswith('1') == True:
		n1merged = 1

	for net2 in networks:

		net2 = net2.strip("\n")

		if net2 in seen_list:
			continue
		
		else:


			##check if merged modules or not
			if net2.endswith('0') == True:
				n2merged = 0
			elif net2.endswith('1') == True:
				n2merged = 1

			bash = open(out_dir + "/" + str(count) + "_bash.sh", 'w')

			bash.write("#PBS -q high_throughput\n#PBS -l mem=8GB,walltime=6:00:00\n")
			bash.write("/cbcb/lab/smount/programs/R-3.1.2/bin/Rscript /cbcb/lab/smount/ZCL/bioconductor_scripts/network_pair_compare.R " + net + " " + str(n1merged) + " " + net2 + " " + str(n2merged) + " " + out_dir+"/results" + "\n")

			bash.close()

			big_bash.write("qsub " + str(count) + "_bash.sh\n")


			count += 1

	seen_list.append(net)

big_bash.close()
network_file.close()

os.makedirs(out_dir+"/results")
os.makedirs(out_dir+"/results/failed")
os.makedirs(out_dir+"/results/success")



