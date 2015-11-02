
from itertools import combinations
import pickle
from splinter import Browser
from splinter import exceptions
import pdb as pdbi
import os.path


#
# Functions
#
def check_for_matching_interaction(p1,p2):
	"""

	Checks if two pdb structures are found to contain pfam domains interacting according to ipfam.

	Note: this relies on a bad assumption. Ideally, there needs to be a distinction between interactors belonging to the native protein, vs. external interactors.

		e.g. check_for_matching_interactions('1PN9','2IL3')

	"""

	p1_interactions = [sorted(i) for i in pdb_to_interactions[p1]]
	p2_interactions = [sorted(i) for i in pdb_to_interactions[p2]]

	for interaction in p1_interactions:
		if interaction in p2_interactions:
			return True

	return False

def check_proteins_have_interacting_pdbs(p1,p2):

	"""

	Checks if two proteins are likely to interact, given that pdb structures between the two proteins are deemed to interact according to "check_for_matching_interaction()". 

		e.g. check_proteins_have_interacting_pdbs('Q17027','Q9BIH3')

	"""

	p1_pdbs = protein_to_pdb[p1]
	p2_pdbs = protein_to_pdb[p2]

	for pdb1,pdb2 in pdb_to_pdb_interactions:

		if (pdb1 in p1_pdbs) and (pdb2 in p2_pdbs):
			return True
		if (pdb1 in p2_pdbs) and (pdb2 in p1_pdbs):
			return True

	return False

def get_protein_to_pdb_dict(species,uniprot_search_file):

	""" 

	Generates a dictionary mapping Uniprot Ids of proteins, to their PDB structures given a results file downloaded from a UniProt search file**.

	Usage:
		protein_to_pdb_dict = get_protein_to_pdb_dict("Anopheles gambiae","./data/uniProt_structures_anopheles_gambiae.tsv")
	Args:
		species (string), e.g. "Anopheles gambiae"
		uniprot_search_file**, tab-delimited output from searching "Anopheles gambiae" in http://www.uniprot.org/
	Returns:
		protein_to_pdb (dict), e.g. {'Q005N3': ['3PZF', '4RO9', '4ROA', '4ROB', '4RSQ'],'Q17027': ['2E3E', '2E3F', '2E3G', '2NY8', '2NY9', '2NZ3']}


	@uniprot_search_file**: 

	Go to http://www.uniprot.org/, and search for a species name, e.g. "Anopheles gambiae". 

	In the results, click on the pencil icon to edit search fields, and include the following fields,
	with the same order. Then click "download" (format: Text):

		Entry
		Organism
		Cross-ref (PDB)

	"""

	fi = open(uniprot_search_file,"r")

	protein_to_pdb = {}

	while True:

		line = fi.readline().replace("\n","")

		if line == "":
			break

		if "Organism" in line:
			headers = line.split("\t")
			for i,header in enumerate(headers):

				if header.rstrip()=="Cross-reference (PDB)":
					pdb_index = i 
				if header.rstrip()=="Organism":
					organism_index = i
				if header.rstrip()=="Entry":
					entry_index = i

		if not species in line:
			continue

		sline = line.split("\t")

		protein  			= sline[entry_index]
		organism 			= sline[organism_index]
		pdb 				= sline[pdb_index]

		if (organism.lower() == species.lower()) and (";" in pdb):

			#print "pdb structure found: "+pdb+" ("+organism+")"

			pdbs = pdb.split(";")

			if "" in pdbs:
				pdbs.remove("")

			if protein_to_pdb.has_key(protein):
				print "woah, something's wrong... multiple rows exist for this protein"

			protein_to_pdb[protein] = pdbs

	fi.close()

	return protein_to_pdb

def merge_dicts(dict_list):
	"""
	Merge two dicts, essential for obtaining: protein_to_pdb

	Args:
		dict_list (list of dicts), e.g. [{A:[1,2],B:[1]}, {A:[1],B:[1]}]
		^^						 , e.g. merge_dicts([protein_to_pdb_anopheles,protein_to_pdb_plasmodium])
	Returns:
		merged_dict

	"""
	import copy

	merged_dict = {}

	for i,dictX in enumerate(dict_list):

		for key in dictX:

			if merged_dict.has_key(key):

				print "merged_dict already has a key from a previous dict: "+str(key)+"...merging by renaming the key to <key>_<i>"

				key = key+"_"+str(i)

			merged_dict[key] = copy.copy(dictX[key])

	return merged_dict

def prepare_interaction_detection_inputs(protein_to_pdb_anopheles,protein_to_pdb_plasmodium):

	""" Takes dictionaries that map proteins to their PDB structures from: species1 and species2, then prepares inputs to the interaction prediction engine: @todo

	Args:
		protein_to_pdb_anopheles (dict), e.g. get_protein_to_pdb_dict("Anopheles gambiae","./data/uniProt_structures_anopheles_gambiae.tsv")
		protein_to_pdb_anopheles (dict), e.g. get_protein_to_pdb_dict("Plasmodium falciparum","./data/uniProt_structures_plasmodium_falciparum.tsv")
	Returns:
		pdb_structures (list), e.g. <a list of annotated structures from PDB, as PDB IDs>
		pdb_structures_query** (str), e.g. <the same as pdb_structures but as a string with each ID delimited by a space>

	@pdb_structures_query**: this query is used to search iPfam for interacting pfam sub-domains.

	"""

	protein_to_pdb = merge_dicts([protein_to_pdb_anopheles,protein_to_pdb_plasmodium])

	pdb_structures = []

	for protein in protein_to_pdb.keys():
		pdbs = protein_to_pdb[protein]
		pdb_structures += pdbs 

	pdb_structures = list(set(pdb_structures))

	pdb_structures_query = " ".join(pdb_structures)

	return pdb_structures, pdb_structures_query, protein_to_pdb

def query_iPfam( pdb_structures_query ):

	#
	# open browser
	#

	br = Browser()

	url = 'http://www.ipfam.org/search/keyword'

	br.visit(url)

	#  
	# Search pdb structures vs. interactions
	#

	# make a search qeury with all the pdb structures
	br.find_by_css("#keywords")[0].fill(pdb_structures_query)

	br.find_by_css("input.button").click()

	# all structure interactions
	br.find_by_css(".lozenge > ul:nth-child(2) > li:nth-child(3) > input:nth-child(1)").click()

	# all ligand interactions
	# ...

	# click "show all"
	br.find_by_css("input.button:nth-child(3)").click()

	# show 100 entries
	br.find_by_css("#pdb_matches_table_length > label:nth-child(1) > select:nth-child(1)").first.select("-1")

	# grab all structure's and their interactions links
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
	print "obtaining interactions for each pdb structure..."

	pdb_to_interactions = {}

	interaction_to_url = {}

	for pdb, url in pdb_to_url:

		print "pdb structure: "+pdb

		br.visit(url)

		interaction_status = br.find_by_css("div.lozenge:nth-child(1) > dl:nth-child(3) > dd:nth-child(2) > p:nth-child(1) > label:nth-child(2)").first.text

		n_family_interactions = int(interaction_status.replace("Family (","").replace(")",""))

		if n_family_interactions > 0:
			print "\t\t"+str(n_family_interactions)+" interactions found"

			br.find_by_value("fam_int").first.click() # click family interactions

			family_interactions = br.find_link_by_partial_href("/fam_int/") # @todo: test if this is a correct matcher

			for interaction in family_interactions:

				interaction_url = interaction['href']
				a, b = interaction_url.split("/fam_int/")
				a_pfam_id = a.split("/family/")[1]
				b_pfam_id = b.split("/sequence")[0]

				interaction_neat = (a_pfam_id,b_pfam_id)

				print "\t\t\titeraction: "+interaction_neat[0]+"-to-"+interaction_neat[1]+" url: "+interaction['href'] # e.g. RVP-to-RVP

				interaction_to_url[interaction_neat] = interaction['href']

				if pdb_to_interactions.has_key(pdb):
					pdb_to_interactions[pdb].append(interaction_neat)
				else:
					pdb_to_interactions[pdb] = [interaction_neat]
		else:
			print "\t\t"+str(n_family_interactions)+" interactions found"

			pdb_to_interactions[pdb] = []

	# # save interactions data
	# pickle.dump( pdb_to_interactions, open( "./data/pdb_to_interactions.p", "wb" ) )
	# pickle.dump( interaction_to_url,  open( "./data/interaction_to_url.p", "wb" ) )
	#
	# determine which pdb protein structures interact
	# 	Note: problem, we do not know which of the interacting pfams belong to the native protein
	#

	return pdb_to_interactions, interaction_to_url

def predict_pdb_to_pdb_interactions(pdb_structures):

	""" Determine which pdb structures might interact, by the interacting sub-domains shared between pdb structures

	Args:
		@todo 
	Return:
		@todo 
	...



	"""

	# pdb_structures 			= pdb_structures_query.split(" ")
	pdb_pairs 				= combinations(pdb_structures,2)
	pdb_to_pdb_interactions = []

	for p1,p2 in pdb_pairs:
		try:
			if check_for_matching_interaction(p1,p2):
				pdb_to_pdb_interactions.append((p1,p2))
		except KeyError:
			continue

	return pdb_to_pdb_interactions

def predict_protein_to_protein_interactions(protein_to_pdb,pdb_to_pdb_interactions):

	""" Determine which proteins might interact, based on the observation that two proteins share pfam domains shown to interact according to ipfam.

	Args:
		protein_to_pdb (dict), e.g. see: prepare_interaction_detection_inputs()


	"""
	#fo = open("./data/protein_to_protein_interactions.tsv",'w')

	protein_to_protein_interactions = []

	proteins = protein_to_pdb.keys()

	protein_pairs = combinations(proteins,2)

	#fo.write("Protein pairs that have PDB structures shown to previously interact with one another in iPfam:")

	#fo.write("Protein A\tProtein B\n")

	for p1,p2 in protein_pairs:
		if check_proteins_have_interacting_pdbs(p1,p2):
			protein_to_protein_interactions.append((p1,p2))
			#fo.write(p1+"\t"+p2+"\n")

	return protein_to_protein_interactions

def write_outputs():

	""" Write outputs

	"""

	#
	# Write protein interactors to files
	# 

	fo = open("../data/protein_to_protein_interactions.tsv",'w')

	fo.write("Protein pairs that have PDB structures shown to previously interact with one another in iPfam:")

	fo.write("Protein A\tProtein B\n")

	for p1,p2 in protein_to_protein_interactions:
		fo.write(p1+"\t"+p2+"\n")

	fo.close()

	#
	# Write other data to files
	#
	print("writing data to files...")

	fo = open("../data/pdb_to_interactions.tsv","w")

	fo.write("PDB structures and their interactions, shown by pairwise pfam ids:\nInteractor A\tInteractor B\n")

	fo.write("PDB Structure\tPfam interactor A\tPfam interactor B\n")

	for pdb in pdb_to_interactions.keys():

		for interaction in pdb_to_interactions[pdb]:
				a,b = interaction 
				fo.write(pdb+"\t"+a+"\t"+b+"\n")

	fo.close()

	#
	# URLs to interactions
	#
	fo = open("../data/interaction_to_url.tsv","w")

	fo.write("Pfam interactions observed in ./pdb_to_interactions.tsv, and their URLs:\n")

	for interaction in interaction_to_url.keys():

		a,b = interaction

		url = interaction_to_url[interaction]

		fo.write(a+"\t"+b+"\t"+url+"\n")

	fo.close()

########
# Main #
########

#
# Parse in UniProt search results
#
protein_to_pdb_anopheles 	= get_protein_to_pdb_dict("Anopheles gambiae","../data/uniProt_structures_anopheles_gambiae.tsv")
protein_to_pdb_plasmodium 	= get_protein_to_pdb_dict("Plasmodium falciparum","../data/uniProt_structures_plasmodium_falciparum.tsv")

#
# Prepare data
#
pdb_structures, pdb_structures_query, protein_to_pdb = prepare_interaction_detection_inputs(protein_to_pdb_anopheles,protein_to_pdb_plasmodium)

#
# Search iPfam for interactions between PDB structures
#
if os.path.exists("../data/cache/pdb_to_interactions.p"):
	print "iPfam results cached, loading..."
	pdb_to_interactions = pickle.load(open('../data/cache/pdb_to_interactions.p',"r"))
	interaction_to_url 	= pickle.load(open('../data/cache/interaction_to_url.p',"r"))
else:
	print "Generating iPfam results..."
	# determine interacting pFam subdomains from PDB structures
	pdb_to_interactions, interaction_to_url = query_iPfam( pdb_structures_query )
	# save interactions data
	pickle.dump( pdb_to_interactions, open( "../data/cache/pdb_to_interactions.p", "wb" ) )
	pickle.dump( interaction_to_url,  open( "../data/cache/interaction_to_url.p", "wb" ) )

#
# Determine which PDB structures are found to interact in iPfam
#
print "determining which pdb structures might interact..."
pdb_to_pdb_interactions = predict_pdb_to_pdb_interactions(pdb_structures)

#
# Determine which Proteins are found to interact <= those that contain PDB interactions
#
print "determining which proteins might interact..."
protein_to_protein_interactions = predict_protein_to_protein_interactions(protein_to_pdb,pdb_to_pdb_interactions)

print "writing data to file..."
write_outputs()


