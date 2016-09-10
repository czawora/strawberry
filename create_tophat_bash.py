'''/cbcb/sw/RedHat-7-x86_64/common/local/tophat/2.1.0/bin/tophat2 -o /cbcb/project-scratch/ZCL/LCM_mapped_reads_tophat/ --max-intron-length 50000 /cbcb/lab/smount/ZCL/Fragaria_vesca_v2.0 /cbcb/project-scratch/ZCL/LCM-raw_reads_5bp_trim/Anther_6-7_A_5bp_trim.fq

'''

import sys

if len(sys.argv) < 3:

	print "first argument must be the input directory for the read files"
	print "there must be a filenames.txt in the at directory"
	print "second parameter must be the output directory for tophat"

	exit(0)

index_file = "/cbcb/lab/smount/ZCL/Fragaria_vesca_v2.0"
input_directory = sys.argv[1]
output_directory = sys.argv[2]

filename_list = input_directory + "filenames.txt"

count = 1

files = open(filename_list, 'r')
bash_filenames = open("bash_filenames.txt", 'w')

for line in files:

	bash = open("bash_" + str(count)+ ".sh", 'w')
	bash_filenames.write("bash_" + str(count)+ ".sh\n")

	bash.write("export PATH=/cbcb/lab/smount/programs/bowtie2-2.2.6:$PATH\n")
	bash.write("/cbcb/sw/RedHat-7-x86_64/common/local/tophat/2.1.0/bin/tophat2 -o "
		+ output_directory + line.strip(".fastq\n") 
		+ " --max-intron-length 50000"
		+ " /cbcb/lab/smount/ZCL/Fragaria_vesca_v2.0 "
		+ input_directory + line)

	bash.close()

	count += 1

bash_filenames.close()
files.close()

