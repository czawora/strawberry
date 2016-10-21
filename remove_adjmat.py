
import subprocess
import sys

exp_name = sys.argv[1]
file_name = sys.argv[2]

f = open("/cbcb/project-scratch/ZCL/wgcna/consensus/" + exp_name + "/" + file_name)

for line in f:

	subprocess.call(["rm", "/cbcb/project-scratch/ZCL/wgcna/consensus/" + exp_name + "/adjmat/" + line.strip("\n")])


f.close()

