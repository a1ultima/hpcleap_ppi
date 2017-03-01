
proteins = []
fi = open("../data/protein_to_pdb.tsv","r")
while True: 
	line = fi.readline()
	if "Protein" in line:
		continue
	if line=="":
		break
	line_split = line.split("\t")
	protein = line_split[0]
fi.close()


uniprot_to_vb = {}
fi = open("../data/UniProt_to_VB_ids.tsv","r")
while True: 
	line = fi.readline()
	if line=="":
		break
	if "Uniprot" in line:
		continue
	uniprot, vb = line.split("\t")
	uniprot_to_vb[uniprot] = vb
fi.close()




# fi = open("../data/7165.protein.links.detailed.v10.txt","r")

# while True: 

# 	line = fi.readline()

# 	if line=="":
# 		break

# 	line_split = line.split("\t")

# fi.close()



