

import sys

if len(sys.argv) < 4:

	print "first parameter must be path to filename list"
	print "second parameter must be where files in filename_list are located eg /some/dir/xxx/"
	print "third parameter must be path of directory fastqc output should be directed to"
	

	exit(0)

filename_list = sys.argv[1]
output_directory = sys.argv[3]
filename_dir = sys.argv[2]

count = 1

files = open(filename_list, 'r')
bash_filenames = open("bash_filenames.txt", 'w')

for line in files:

	bash = open("bash_" + str(count)+ ".sh", 'w')
	bash_filenames.write("bash_" + str(count)+ ".sh\n")

	bash.write("/cbcb/lab/smount/programs/FastQC/fastqc --outdir=" + output_directory + " " + filename_dir + line)

	bash.close()

	count += 1

bash_filenames.close()
files.close()
