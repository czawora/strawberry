
import sys

if len(sys.argv) < 5:

	print "fourth parameter must be name of output file"
	print "third parameter must be path to filename list"
	print "second parameter must be -l parameter"
	print "first parameter must be q sub -q parameter"

	exit(0)

q_param = sys.argv[1]
l_param = sys.argv[2]
shell_filename_list = sys.argv[3]
output_file = sys.argv[4]

shells = open(shell_filename_list,'r')

bash = open(output_file, 'w')

for line in shells:

	bash.write("qsub -q " + q_param + " -l " +  l_param + " " + line)

bash.close()