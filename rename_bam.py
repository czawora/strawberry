

import sys
import os

if len(sys.argv) < 3:

	print "first argument must be path name to dir name list"
	print "second argument must be base dir of tophat output folders"
	exit(0)


dir_names = sys.argv[1]
base_dir = sys.argv[2]

names = open(dir_names)
bam_paths = open(base_dir + "bam_paths.txt", 'w')


for line in names:

	line = line.strip("\n")
	old = base_dir + line + "/accepted_hits.bam" 
	new = base_dir + line + "/" + line+ ".bam"
	print old, new

	try:
		os.rename(old, new)
	except:
		pass

	bam_paths.write(base_dir + line + "/" + line+ ".bam\n")

bam_paths.close()
names.close()
