

files  = open("xml_files.txt")
combined  = open("half1_xml.xml", 'w')

combined.write("<?xml version=\"1.0\"?>\n")
combined.write("<BlastXML2\n")
combined.write("xmlns=\"http://www.ncbi.nlm.nih.gov\"\n")
combined.write("xmlns:xs=\"http://www.w3.org/2001/XMLSchema-instance\"\n")
combined.write("xs:schemaLocation=\"http://www.ncbi.nlm.nih.gov http://www.ncbi.nlm.nih.gov/data_specs/schema_alt/NCBI_BlastOutput2.xsd\"\n")
combined.write(">\n")

count = 1
for f in files:

	count += 1

	if count == 175:
		combined.write("</BlastXML2>\n")
		combined.close()
		combined = open("half2_xml.xml", 'w')
		combined.write("<?xml version=\"1.0\"?>\n")
		combined.write("<BlastXML2\n")
		combined.write("xmlns=\"http://www.ncbi.nlm.nih.gov\"\n")
		combined.write("xmlns:xs=\"http://www.w3.org/2001/XMLSchema-instance\"\n")
		combined.write("xs:schemaLocation=\"http://www.ncbi.nlm.nih.gov http://www.ncbi.nlm.nih.gov/data_specs/schema_alt/NCBI_BlastOutput2.xsd\"\n")
		combined.write(">\n")

	current = open(f.strip("\n"))

	for line in current:

		if line != "<?xml version=\"1.0\"?>\n" and "<BlastXML2" not in line and "xmlns=\"http://www.ncbi.nlm.nih.gov\"" not in line and "xmlns:xs=\"http://www.w3.org/2001/XMLSchema-instance\"" not in line and "xs:schemaLocation=\"http://www.ncbi.nlm.nih.gov http://www.ncbi.nlm.nih.gov/data_specs/schema_alt/NCBI_BlastOutput2.xsd\"" not in line and line != ">\n" and line != "</BlastXML2>\n":
			combined.write(line)

	current.close()

combined.write("</BlastXML2>\n")
combined.close()
files.close()