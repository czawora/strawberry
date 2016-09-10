
#creates bash file for running in the same directory as the bowtie index

#first command line item is path to filename list
#second command line item is path for sam file

import sys
import os

if len(sys.argv) < 4:
	print "please provide all command line arguments"
	print "the first is the path to the filename list file"
	print "the second is the path to directory where you would the output SAM files to go"
	print "the last is path and name of the bash file for this script to produce"
	exit(0)

filenames = open(sys.argv[1], 'r')

count = 1

for line in filenames:

	bash = open(str(count) +"_"+ sys.argv[3], 'w')

	line = line.strip('\n')
	splits = line.split("/")

	#removing file extension
	dot_split = splits[len(splits) -1].split(".")
	#want to glue back all but last split string aka dont want the file extension
	extensionless_filename = "".join(dot_split[:len(dot_split)-1])

#/cbcb/sw/RedHat-7-x86_64/common/local/bowtie2/2.2.4/bin/bowtie2
	bash.write("/cbcb/lab/smount/programs/bowtie2-2.2.6/bowtie2 -x /cbcb/lab/smount/ZCL/Fragaria_vesca_v2.0" 
		+ " -U " + "/cbcb/lab/smount/ZCL/data"+line 
		+ " -S " + "/cbcb/lab/smount/ZCL/"+sys.argv[2] + extensionless_filename + ".sam"#add SAM file extension
		+ "\n")

	count += 1

	bash.close()

