
gtf = open("/Users/Chris/Documents/strawberry/ZCL/Fragaria_vesca_v2.0.a1.transcripts.gtf", 'r')
introns = open("/Users/Chris/Documents/strawberry/ZCL/Fragaria_vesca_v2.0.a1.transcripts_with_introns.gtf", 'w')


for line in gtf:

	splits = line.split("\t")

	#check gene id

	extras = splits[8]
	gene_id = extras.split(";")[1].split(" ")[2]

	print(gene_id)

