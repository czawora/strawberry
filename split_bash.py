

import sys
import os

exp = sys.argv[1]

fls = os.listdir("/cbcb/project-scratch/ZCL/wgcna/consensus/" + exp.strip("\n") + "/bash")

big_count = 0
total_count = 1

big_bash = open("/cbcb/project-scratch/ZCL/wgcna/consensus/" + exp.strip("\n") + "/big_bash" + str(big_count) + ".sh", 'w')

for f in fls:

	if total_count % 300 == 0:

		big_bash.close()
		big_count += 1
		big_bash = open("/cbcb/project-scratch/ZCL/wgcna/consensus/" + exp.strip("\n") + "/big_bash" + str(big_count) + ".sh", 'w')

	big_bash.write("qsub " + "bash/" + f.strip("\n") + "\n")
	total_count += 1


big_bash.close()
