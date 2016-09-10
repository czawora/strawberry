

import sys

if len(sys.argv) < 3:

	print "first argument must be the path and name of filenames list"
	print "second argument must be the destination directory for new folders being made"

	exit(0)


filenames = sys.argv[1]
output_directory = sys.argv[2]

files = open(filenames)

bash = open(output_directory+"make_directories.sh", 'w')

for line in files:

	line = line.strip(".fq\n")
	line = line + "\n"

	bash.write("mkdir " +  output_directory + line)

bash.close()
files.close()