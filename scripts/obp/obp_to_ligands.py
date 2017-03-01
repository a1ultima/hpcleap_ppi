
from itertools import combinations
import pickle
from splinter import Browser
from splinter import exceptions
import pdb as pdbi
import os.path
import time


if not os.path.exists("../../data/obp_input/"):
	os.makedirs("../../data/obp_input/")

#
# get uniprot ids of obps 
# 
# obtain a list of Anopheles gambiae OBP protein ids (uniProt and AGAP style) via http://www.bo-protscience.fr/mobpdb/?page_id=82
#
obps_uniprot 	= []
obps_geneid 	= []
obps_name 		= []

uniprot_to_geneid = {}

fi = open("../../data/obp_input/anopheles-gambiae_obp_ids.tsv","r")

fi.readline() # headers

while True: 

	line = fi.readline()

	if line=="":
		break

	line_split = line.split("\t")

	uniprot = line_split[2]
	geneid 	= line_split[0]
	name 	= line_split[1]

	uniprot_to_geneid[uniprot] = geneid

	obps_uniprot.append(uniprot)
	obps_geneid.append(geneid)
	obps_name.append(name)
	
fi.close()


#
# get uniprot ids -to- pdb structures
#
uniprot_to_pdb = {}

fi = open("../../data/protein_to_pdb.tsv","r")

fi.readline() # headers

while True: 

	line = fi.readline()

	if line=="":
		break

	line_split 	= line.split("\t")
	uniprot 	= line_split[0].rstrip()

	if uniprot in obps_uniprot:
		pdbs = [i.rstrip() for i in line_split[1:]]
		uniprot_to_pdb[uniprot] = pdbs

fi.close()

#
# query ipfam for pdb structures to obtain their ligands (pdb-to-ligand)
#
obp_pdbs = []

for uniprot in uniprot_to_pdb:
	obp_pdbs += uniprot_to_pdb[uniprot]

#
# open browser
#
br = Browser()

url = 'http://www.ipfam.org/search/keyword'

br.visit(url)

#  
# Search pdb structures vs. interactions
#
pdb_structures_query = " ".join(obp_pdbs)


# make a search qeury with all the pdb structures
br.find_by_css("#keywords")[0].fill(pdb_structures_query)

br.find_by_css("input.button").click()

# all structure interactions
br.find_by_css(".lozenge > ul:nth-child(2) > li:nth-child(3) > input:nth-child(1)").click()

# click "show all"
try:
	br.find_by_css("input.button:nth-child(3)").click()
except AttributeError:
	pass

# show 100 entries
br.find_by_css("#pdb_matches_table_length > label:nth-child(1) > select:nth-child(1)").first.select("-1")

#
# Collect all structure's and their interactions links, store into pdb_to_url
#
count = 0

pdb_to_url = []

while True:

	count += 1

	try: 
		pdb_id = br.find_by_css("#pdb_matches_table > tbody:nth-child(2) > tr:nth-child("+str(count)+") > td:nth-child(1) > a:nth-child(1)").first.text
		pdb_url = br.find_by_css("#pdb_matches_table > tbody:nth-child(2) > tr:nth-child("+str(count)+") > td:nth-child(1) > a:nth-child(1)").first['href']
		pdb_to_url.append((pdb_id,pdb_url))
	except exceptions.ElementDoesNotExist: 
		break



#
# obtain interactions per pdb
#
print "determining ligands binding to each Odorant Binding Protein..."

pdb_to_interactions = {}

interaction_to_url = {}

pdb_to_ligand = {} # e.g. {<pdb>:{ 'url': ... , 'ligand': ... } , <pdb2>: ... }

for pdb, url in pdb_to_url:

	print "pdb structure: "+pdb

	br.visit(url)

	# Copies the text next to the radio button: "Ligands (4)", in the "Annotations:" section
	interaction_status = br.find_by_css("div.lozenge:nth-child(1) > dl:nth-child(1) > dd:nth-child(2) > p:nth-child(2) > label:nth-child(2)").first.text

	# Collect ligand URLs for each pdb (later we need to visit them to grab the ligand names)
	n_ligands = int(interaction_status.replace("Ligands (","").replace(")",""))
	
	if n_ligands > 0:

		# add the url to the pdb_to_ligand dict, for later visiting in order to find the name of the ligand
		ligands = br.find_link_by_partial_href("/ligand/")

		pdb_to_ligand[pdb] = {}

		for ligand in ligands:

			ligand_url = ligand['href']

			ligand_id = ligand_url.replace("http://www.ipfam.org/ligand/","")

			pdb_to_ligand[pdb][ligand_id] = {"url":ligand_url,"name":"this should have been replaced! go back to @ligandName to debug"}			

	else:
		print "\t\t0 interactions found"

		pdb_to_ligand[pdb] = {}

#
# Visit each ligand URL to obtain their names
#
print "Visiting Ligand URLS to extract their names..."

for pdb in pdb_to_ligand.keys():

	print pdb

	for ligand in pdb_to_ligand[pdb].keys():

		print "\t"+ligand

		br.visit(pdb_to_ligand[pdb][ligand]['url'])

		pdb_to_ligand[pdb][ligand]['name'] = br.find_by_css(".buffer_top").first.text

		try:
			ligand_external_url_PDB_RCSB = br.find_by_css("div.left > section:nth-child(2) > div:nth-child(2) > div:nth-child(1) > p:nth-child(2) > a:nth-child(1)").first['href']
			ligand_external_url_PDBe = br.find_by_css("div.left > section:nth-child(2) > div:nth-child(2) > div:nth-child(2) > p:nth-child(2) > a:nth-child(1)").first['href']

			pdb_to_ligand[pdb][ligand]['ext_url_pdb_rcsb'] = ligand_external_url_PDB_RCSB
			pdb_to_ligand[pdb][ligand]['ext_url_pdbe'] = ligand_external_url_PDBe
		except:
			pass

#
# Assign Ligand names to UniProt Protein IDs 
#
print "Assigning Ligand names to UniProt Protein IDs of Odorant Binding Proteins..."

protein_to_ligands = {}

for protein in uniprot_to_pdb:

	pdbs = uniprot_to_pdb[protein]

	for pdb in pdbs:

		ligands = pdb_to_ligand[pdb].keys()

		for ligand in ligands:

			ligand_name = pdb_to_ligand[pdb][ligand]['name']

			try: 
				protein_to_ligands[protein].append(ligand_name)
			except KeyError:	
				protein_to_ligands[protein] = [ligand_name]


#
# Wwrite to file 
#
if not os.path.exists("../../data/obp_output"):
 	os.makedirs("../../data/obp_output")

fo = open("../../data/obp_output/iPfam_protein_to_ligands.txt","w")

headers = "UniProt\tGeneId\t"+"\t".join(["Ligand "+str(i) for i in range(0,11)])+"\n"



fo.write(headers)

for protein in protein_to_ligands.keys():

	ligands = protein_to_ligands[protein]

	geneid = uniprot_to_geneid[protein]

	fo.write(protein+"\t"+geneid+"\t"+"\t".join(ligands)+"\n")

fo.close()

print "ANDY REMEMBER!!!: Please check how many of all our OBP proteins have ligands available in iPfam, it may be that iPfam only has a small fraction of them... "

