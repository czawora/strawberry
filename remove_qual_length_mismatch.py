

import sys

if len(sys.argv) == 1:

	print "argument must be file to edit"
	exit(0)

inputfile = sys.argv[1]

f = open(inputfile)
output = open("edited_" + inputfile , 'w')

now_seq = False
now_qual = False
current_seg_length = None
current_qual_length = None
count = 0
tmp = []

for line in f:

	if line[0] == '@':#new block

		tmp = []
		tmp.append(line)
		now_seq = True

	elif now_seq == True:

		now_seq = False
		current_seg_length = len(line)
		tmp.append(line)

	elif line[0] == '+':

		tmp.append(line)
		now_qual = True


	elif now_qual == True:

		now_qual = False
		current_qual_length = len(line)

		tmp.append(line)

		if current_qual_length == current_seg_length:

			for item in tmp:

				output.write(item)

		else:

			count += 1

			print "removed block: "
			print tmp

print "blocks removed: " + str(count)

f.close()
output.close()

