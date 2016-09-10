
import sys

if len(sys.argv) < 5:

	print "first parameter must be path to filename list"
	print "second parameter must be where files in filename_list are located eg /some/dir/xxx/"
	print "third parameter must be path of directory trimomatic output should be directed to"
	print "fourth parameter must be the number of base pairs to HEADCROP"
	

	exit(0)

filename_list = sys.argv[1]
output_directory = sys.argv[3]
filename_dir = sys.argv[2]
crop_num = sys.argv[4]

count = 1

files = open(filename_list, 'r')
bash_filenames = open("bash_filenames.txt", 'w')

for line in files:

	output_filename_split = line.strip("\n").split(".")
	output_filename_split[len(output_filename_split) - 2] = output_filename_split[len(output_filename_split) - 2] + "_" + crop_num + "bp_trim"
	output_filename = ".".join(output_filename_split)

	bash = open("bash_" + str(count)+ ".sh", 'w')
	bash_filenames.write("bash_" + str(count)+ ".sh\n")
	bash.write("java -jar /cbcb/lab/smount/programs/trimmomatic-0.36.jar SE " 
		+ filename_dir+ line.strip("\n") 
		+ " " +  output_directory + output_filename  
		+ " HEADCROP:" + crop_num + "\n")

	bash.close()

	count += 1

bash_filenames.close()
files.close()
