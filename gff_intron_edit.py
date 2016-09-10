

def find_introns(gene_range):

	pseudo_cds_list = []

	prev_idx = None
	pseudo_cds = []

	for idx in gene_range:

		if prev_idx != None:
			if idx - prev_idx > 1:
				pseudo_cds_list.append((pseudo_cds[0], pseudo_cds[len(pseudo_cds)-1]))
				pseudo_cds = []

		pseudo_cds.append(idx)
		prev_idx = idx

	if pseudo_cds != []:
		pseudo_cds_list.append((pseudo_cds[0], pseudo_cds[len(pseudo_cds)-1]))

	return pseudo_cds_list

gff = open("/Users/Chris/Documents/strawberry/ZCL/Fragaria_vesca_v2.0.a1.transcripts.gff3", 'r')
introns = open("/Users/Chris/Documents/strawberry/ZCL/Fragaria_vesca_v2.0.a1.transcripts_cds_as_introns.gtf", 'w')
# RNA_count_f = open("/Users/Chris/Documents/strawberry/scripts/RNA_count.txt", 'w')
lines = [l for l in gff]

saved_gene_line = None
saved_mRNA_line = None
gene_range = None
read_new_data = False
# RNA_count = 0
line_count = 0
none_count = 0
gene_count = 0
present_count = 0

for line in lines:

	line_count += 1
	print(str(line_count) + "/" + str(len(lines)))

	splits = line.split("\t")

	#check is gene 

	if splits[2] == "gene":

		gene_count += 1
		#deal with previosuly saved gene and write it out with introns

		if read_new_data == True and saved_gene_line != None:

			found = find_introns(gene_range)
			if found == []:
				none_count += 1
			else:
				present_count += 1

			
			# save_splits = saved_gene_line.split("\t")
			# gene_id_name = save_splits[8].split(";")[0].split("=")[1]
			# transcript_id = saved_mRNA_line.split("\t")[8].split(";")[0].split("=")[1]

			# for new_cds in find_introns(gene_range):
			# 	introns.write(save_splits[0] + "\t" + save_splits[1] + "\t" + "CDS" + "\t" + str(new_cds[0]) + "\t" + str(new_cds[1]) + "\t" + save_splits[5] + "\t" + save_splits[6] + "\t" + "0" + "\t" + "transcript_id \"" + transcript_id + "\"; " + "gene_id \"" + gene_id_name + "\"; " + "gene_name \"" + gene_id_name + "\";\n")

		start = int(splits[3])
		stop = int(splits[4])

		gene_range = range(start, stop + 1)

		saved_gene_line = line
		RNA_count = 0

	#check if mRNA

	if splits[2] == "mRNA":

		# RNA_count +=1
		# if RNA_count > 1:
		# 	RNA_count_f.write(line)

		saved_mRNA_line = line

	#check if CDS

	if splits[2] == "CDS":

		cds_start = int(splits[3])
		cds_stop = int(splits[4])

		for num in range(cds_start, cds_stop + 1):
			gene_range.remove(num)

		read_new_data = True

##write out the introns for the last gene

# save_splits = saved_gene_line.split("\t")
# gene_id_name = save_splits[8].split(";")[0].split("=")[1]
# transcript_id = saved_mRNA_line.split("\t")[8].split(";")[0].split("=")[1]

# for new_cds in find_introns(gene_range):
# 	introns.write(save_splits[0] + "\t" + save_splits[1] + "\t" + "CDS" + "\t" + str(new_cds[0]) + "\t" + str(new_cds[1]) + "\t" + save_splits[5] + "\t" + save_splits[6] + "\t" + "0" + "\t" + "transcript_id \"" + transcript_id + "\"; " + "gene_id \"" + gene_id_name + "\"; " + "gene_name \"" + gene_id_name + "\";\n")

print("total genes: " + str(gene_count))
print("gene with introns: " + str(present_count))
print("genes without introns: " + str(none_count))

print "complete"
#read each line
#look for the CDS gaps and add Intron markers there
