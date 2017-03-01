# -*- coding: utf-8 -*-

import pdb

import difflib
import numpy as np
import csv

from matplotlib import pyplot as plt
import numpy as np
from matplotlib_venn import venn3, venn3_circles, venn2

#############
# Functions #
#############

def compare(a,b):

	""" Compares two strings for equality, returning a score: 0-to-1

	"""

	return difflib.SequenceMatcher(None, a, b).ratio()

def remove_non_mutual_ligands(mobp_sorted_ki_ligands, qobp_sorted_ki_ligands_wrong_spelling,intersection_mobp):
	""" Remove ligands from mobp and qobp ligand-vs-obp binding assays that are not mutual to both experiments
	"""

	quSpelling_to_mobpSpelling= \
		{'(-)-citronellal': 'citronellal',\
		 '(e)-2-hexenal': 'e2-hexenal',\
		 '6-methyl-5-hepten-2-one': '6-methylhept-5-en-2-one',\
		 'geranylacetone': 'geranyl acetone',\
		 "n,n-diethyl-m-toluamide (deet)": "deet"}

	#
	# mOBP: ignore ligands not present in the mutually-tested ligands list
	#
	mobp4_sorted_ki_ligands_mutual = []

	for ligand in mobp_sorted_ki_ligands:
		if (ligand in intersection_mobp):
			mobp4_sorted_ki_ligands_mutual.append(ligand)

	#
	# Qaio et al: ignore ligands not present in the mutually-tested ligands list
	#

	qobp4_sorted_int_ligands_wrong_spelling = [i.lower() for i in qobp_sorted_ki_ligands_wrong_spelling]

	qobp4_sorted_int_ligands = []

	# spellings 
	for ligand in qobp4_sorted_int_ligands_wrong_spelling:
		try:
			ligand_mobp = quSpelling_to_mobpSpelling[ligand]
		except KeyError:
			ligand_mobp = ligand
		qobp4_sorted_int_ligands.append(ligand_mobp)

	qobp4_sorted_int_ligands_mutual = []

	for ligand in qobp4_sorted_int_ligands:
		if (ligand in intersection_mobp):
			qobp4_sorted_int_ligands_mutual.append(ligand)
	#
	# mOBP: ignore ligands with Int values not available in the qOBP list
	#
	mobp4_sorted_ki_ligands_mutual_final = []

	for ligand in mobp4_sorted_ki_ligands_mutual:
		if (ligand in qobp4_sorted_int_ligands_mutual):
			mobp4_sorted_ki_ligands_mutual_final.append(ligand)

	return mobp4_sorted_ki_ligands_mutual_final, qobp4_sorted_int_ligands_mutual

def read_mobp_data(mobp_path):
	"""
		Returns lists of ligands in order of which is most likely to bind to which OBP

		mobp_path = "../../data/obp_input/mOBP4.txt"

	"""

	#
	# Read .csv file representing the punett
	#
	reader      = csv.reader(open(mobp_path,"rb"),delimiter="\t")
	x           = list(reader)
	mobp4      	= np.array(x)

	# 'FEB',   	7
	# 'Ki', 	8
	# 'pKi', 	9
	# 'LE', 	10
	# 'FQ', 	11
	# 'SILE'	12

	mobp4 = mobp4[1:]

	mobp4_sorted_ki 	= mobp4[mobp4[:,8].argsort()] 	# sorted by x-th column
	mobp4_sorted_pki 	= mobp4[mobp4[:,9].argsort()[::-1]] 	# sorted by x-th column
	mobp4_sorted_le 	= mobp4[mobp4[:,10].argsort()[::-1]] 	# sorted by x-th column
	mobp4_sorted_fq 	= mobp4[mobp4[:,11].argsort()[::-1]] 	# sorted by x-th column
	mobp4_sorted_sile 	= mobp4[mobp4[:,12].argsort()[::-1]] 	# sorted by x-th column

	mobp4_sorted_ki_ligands 	= list(mobp4_sorted_ki[:,1])
	mobp4_sorted_pki_ligands 	= list(mobp4_sorted_pki[:,1])
	mobp4_sorted_le_ligands 	= list(mobp4_sorted_le[:,1])
	mobp4_sorted_fq_ligands 	= list(mobp4_sorted_fq[:,1])
	mobp4_sorted_sile_ligands 	= list(mobp4_sorted_sile[:,1])

	mobp4_sorted_ki_ligands = [i.lower() for i in mobp4_sorted_ki_ligands]
	mobp4_sorted_pki_ligands = [i.lower() for i in mobp4_sorted_pki_ligands]
	mobp4_sorted_le_ligands = [i.lower() for i in mobp4_sorted_le_ligands]
	mobp4_sorted_fq_ligands = [i.lower() for i in mobp4_sorted_fq_ligands]
	mobp4_sorted_sile_ligands = [i.lower() for i in mobp4_sorted_sile_ligands]

	return mobp4_sorted_ki_ligands, mobp4_sorted_pki_ligands, mobp4_sorted_le_ligands, mobp4_sorted_fq_ligands, mobp4_sorted_sile_ligands

def read_qobp_data(qobp_path):

	""" Read qobp ligands data

	Args:
		qobp_path = "/home/qiime/Desktop/hpcleap_ppi/data/obp_input/qOBP4.txt"

	Usage: 
		qobp_ligands = read_qobp_data("/home/qiime/Desktop/hpcleap_ppi/data/obp_input/qOBP4.txt")

	"""

	fi = open(qobp_path,"r")

	qobp_ligands = []

	headers = fi.readline()

	while True: 

		line = fi.readline()

		if line=="":
			break

		line_split = line.split("\t")

		ligand = line_split[0].rstrip()

		qobp_ligands.append(ligand)

	fi.close()

	#pdb.set_trace()

	return qobp_ligands

def compare_ligand_sequence(mobp_input,qobp_input,intersection_mobp,m_sort_by="ki"):

	""" Generate sequences of ligands mutually shared between two ligand-vs-obp binding assays, in order of which ligands are most likely to bind the obp

	Args: 
		mobp_input = "../../data/obp_input/mOBP4.txt"

		qobp_input = "/home/qiime/Desktop/hpcleap_ppi/data/obp_input/qOBP4.txt"

		intersection_mobp = ['nonanal', 'indole', '1-dodecanol', '7-octenoic acid', 'octanoic acid', 'heptanoic acid', 'nonanoic acid', 'decanal', 'hexanal', '1-octen-3-ol', 'icaridin', 'benzaldehyde', 'pentanoic acid', 'octanal', '6-methylhept-5-en-2-one', 'citronellal', 'geranyl acetone', 'e2-hexenal', 'deet']

		m_sort_by = "ki"  # or: "pki", "le", "fq", "sile" 

	Usage:

		mobp4_ligands, qob4_ligands = compare_ligand_sequence("../../data/obp_input/mOBP4.txt","/home/qiime/Desktop/hpcleap_ppi/data/obp_input/qOBP4.txt",intersection_mobp)

	"""

	#
	# Read the qobp ligands data
	#
	qobp4_sorted_ki_ligands_wrong_spelling = read_qobp_data(qobp_input)

	#
	# Read the mobp ligands data
	#
	mobp4_ki, mobp4_pki, mobp4_le, mobp4_fq, mobp4_sile = read_mobp_data(mobp_input)

	#
	# Remove non-mutual ligands
	#
	if m_sort_by.lower() == "ki":
		mobp4_ligands, qob4_ligands = remove_non_mutual_ligands(mobp4_ki,qobp4_sorted_ki_ligands_wrong_spelling,intersection_mobp)
	elif m_sort_by.lower() == "pki":
		mobp4_ligands, qob4_ligands = remove_non_mutual_ligands(mobp4_pki,qobp4_sorted_ki_ligands_wrong_spelling,intersection_mobp)
	elif m_sort_by.lower() == "le":
		mobp4_ligands, qob4_ligands = remove_non_mutual_ligands(mobp4_le,qobp4_sorted_ki_ligands_wrong_spelling,intersection_mobp)
	elif m_sort_by.lower() == "fq":
		mobp4_ligands, qob4_ligands = remove_non_mutual_ligands(mobp4_fq,qobp4_sorted_ki_ligands_wrong_spelling,intersection_mobp)
	elif m_sort_by.lower() == "sile":
		mobp4_ligands, qob4_ligands = remove_non_mutual_ligands(mobp4_sile,qobp4_sorted_ki_ligands_wrong_spelling,intersection_mobp)
	else:
		print "incorrect argument for 'm_sort_by' in compare_ligand_sequence()... defaulting to 'ki'..."
		mobp4_ligands, qob4_ligands = remove_non_mutual_ligands(mobp4_ki,qobp4_sorted_ki_ligands_wrong_spelling,intersection_mobp)
	# Out[112]: 
	# ['hexanal',
	#  '6-methylhept-5-en-2-one',
	#  '7-octenoic acid',
	#  'e2-hexenal',
	#  'geranyl acetone',
	#  '1-octen-3-ol',
	#  'benzaldehyde',
	#  'octanoic acid',
	#  'indole',
	#  'citronellal',
	#  'decanal',
	#  '1-dodecanol',
	#  'nonanal',
	#  'octanal']

	# In [113]: mobp4_sorted_ki_ligands_mutual_final
	# Out[113]: 
	# ['e2-hexenal',
	#  'geranyl acetone',
	#  'citronellal',
	#  '1-dodecanol',
	#  'indole',
	#  '6-methylhept-5-en-2-one',
	#  '1-octen-3-ol',
	#  'decanal',
	#  '7-octenoic acid',
	#  'benzaldehyde',
	#  'nonanal',
	#  'octanal',
	#  'octanoic acid',
	#  'hexanal']

	return mobp4_ligands, qob4_ligands, qobp4_sorted_ki_ligands_wrong_spelling

def ladder_plots(obp_id,obp_to_comparison,plot_path,):

	""" Plots ligand preference sequences for mOBP and Qaio et al. with ladder plots to show obp_to_comparison


	Args:
		obp_id =  "OBP4"

		plot_path = "../../data/obp_output/"+obp_id+"_ladder.png"

		obp_to_comparison = obp_to_comparison["OBP"+obp_id] = {"mOBP":mobpX_ligands, "qOBP":qobpX_ligands}

	"""

	#
	# Get coordiantes for plotting on ladder for each ligand of each OBP
	#

	mobp4_ligands = obp_to_comparison[obp_id]['mOBP']
	qobp4_ligands = obp_to_comparison[obp_id]['qOBP']


	ligand_to_coordinate = {}

	for i,ligand in enumerate(mobp4_ligands):

		ligand_to_coordinate[ligand] = {'mOBP':[0,i],'qOBP':[]}

	for i,ligand in enumerate(qobp4_ligands):

		ligand_to_coordinate[ligand]['qOBP'] = [1,i]


	#
	# Plot ladders
	#

	plt.figure()

	for ligand in ligand_to_coordinate.keys():

		x = ligand_to_coordinate[ligand]['mOBP']

		y = ligand_to_coordinate[ligand]['qOBP']

		t = np.array([x,y]).T

		x_vals = list(t[0])

		y_vals = list(t[1])

		plt.plot(x_vals,y_vals,label=ligand)
		plt.text(x_vals[1],y_vals[1],ligand)
		plt.text(x_vals[0],y_vals[0],ligand)

	#plt.plot(x_vals,y_vals,label=ligand)

	#plt.legend(loc="upper right")

	plt.title(obp_id)

	plt.savefig(plot_path)

	plt.close()

def generate_ligands_intersection():

	""" Generates a list of ligands tested in both: Qaio et al. and mOBP

	"""

	#
	# Determining the intersection of Qaio et at. and mOBP
	#
	qu=['(Z)-3-Hexenol','1-Octanol','1-Octen-3-ol','1-Nonanol','1-Decanol','Linalool','Menthol','1-Dodecanol','3,7-Dimethyloctanol','a-Pentylcinnamyl alcohol','Farnesol','Retinol','Aldehydes and ketones','Hexanal','(E)-2-Hexenal','Octanal','(E)-2-Octenal','Nonanal','(E)-2-Nonenal','Decanal','(-)-citronellal','Benzaldehyde','a-Pentyl cinnamaldehyde','6-Methyl-5-hepten-2-one','Jasmone','Geranylacetone','Retinal','Carboxylic acids','Pentanoic acid','Heptanoic acid','Octanoic acid','7-Octenoic acid','Nonanoic acid','Undecanoic acid','Benzoates','Ethyl benzoate','Butyl benzoate','Hexyl benzoate','Octyl benzoate','3,7-Dimethyloctyl benzoate','3-Hexyl benzoate','4-Methylpentyl benzoate','p-Isopropylphenyl benzoate','Butyl p-tert-butylbenzoate','Phenyl benzoate','p-Tolyl benzoate','2-Phenylethyl benzoate','Benzyl benzoate','p-tert-Butylphenyl benzoate','Butyl p-nitrobenzoate','Isobutyl p-nitrobenzoate','Aromatic compounds','m-Cresol','p-Cresol','p-tert-Butylbenzophenone','Phenylbenzylhydrazine','Indole','Methyl cinnamate','Butyl cinnamate','o-Hydroxybenzaldehyde','p-tert-butyl benzaldehyde','N-p-isopropylphenyl-p-','hydroxybenzimine','4-Hydroxy-4’-isopropylazobenzene','2-Pyrrolyl-p-methyl-azobenzene','N,N-diethyl-m-toluamide (DEET)','Others','(E)-b-Farnesene','3,7-Dimethyloctyl acetate','Icaridin']

	mobp=['permethrin','5-alpha-androst-16-en-3alpha-ol','5-alpha-androst-16-one','(1S,4R)-(-)-fenchone','octadecanoic acid','2-ethylphenol','methyl salicylate','acetophenone','hexadecanoic acid','2-acetylpyridine','tetradecanoic acid','octanal','3-methyl-2-cyclohexen-1-ol','methyl benzoate','2-acetylthiophene','2-methylphenol','PMD','benzaldehyde','dodecanoic acid','linalool oxide','geranyl acetone','2-acetylthiazole','tridecanoic acid','3-methylphenol','cyclohexanone','E2-hexenal','decanoic acid','(+)-carvone','2-ethyl toluene','Nepetalactone','Icaridin','cis-3-hexen-1-ol','carvacrol','heptanal','2-propylphenol','(-)-carvone','3-methyl-2-hexenoic acid','eucamalol','4-methylcyclohexanol','Z2-hexenol','geranyl acetate','nonanoic acid','phenol','L(+) lactic acid','2,4,5-trimethyl thiazole','isobutyl acetate','heptane','DEPA','hexanoic acid','DEET','octanoic acid','limonene','isovaleric acid','hexanal','1-hexanol','6-methylhept-5-en-2-one','(+)-fenchone','2-heptanone','cuminyl alcohol','(1R)-(-)-fenchone','gamma-decalactone','Geraniol','indole','octadecanoic acid','delta-decalactone','1-hepten-3-ol','E2-hexenol','ethyl hexanoate','3-methylindole','2-ethoxythiazole','4,5-dimethyl thiazole','heptanoic acid','1-hexen-3-ol','2-phenoxy ethanol','2-oxohexanoic acid','(±)-β-Citronellol','3-methyl-1-butanol','Citronellal','4-ethylphenol','Eucalyptol','2-oxobutanoic acid','2-oxopentanoic acid','phenethyl acetate','benzyl acetate','isobutyric acid','propyl acetate','methyl-2-methyl benzoate','4-methylphenol','1-pentanol','isoamyl acetate','ethyl propanoate','4-methylthiazole','ethyl butyrate','1-dodecanol','pentanoic acid','1-chlorododecane','2-iso-butyl-thiazole','butanoic acid','cadaverine','acetic acid','2-oxopropanoic acid','IR3535','2,3-butanedione','2-butanone','ethyl acetate','2-nonanone','propanoic acid','cadaverine','7-octenoic acid','2-ethyl-1-hexanol','1-butanol','methyl propanoate','decanal','1-octen-3-ol','amyl acetate','3-octanone','methyl caprylate','putrescine','nonanal','ammonia','thiazole','ethyl formate','acetone','methanol','dimethylsulfide','ethanol']

	qu2=[i.lower() for i in qu] 		# lowercase'd ligands list for Qaio et al.  @lower
	mobp2=[i.lower() for i in mobp]		# ^^ for mOBP

	intersection = []
	counts = []

	qu2_to_mobp2 = {}

	count = 0

	# STRICT
	for i in qu2:
		for j in mobp2:

			count += 1

			if compare(i,j) > 0.8:
				#print i + ":::" + j

				if ("anoic acid" in i) and ("anoic acid" in j):
					#print "\t"+i+":::"+j
					if i == j:
						intersection.append(i)
						continue
					continue

				if i + ":::" + j=="1-decanol:::1-dodecanol":
					continue

				if i == "geranylacetone":
					if j == "geranyl acetate":
						continue
				if i == "ethyl benzoate":
					continue

				if i == "octanoic acid":
					continue

				if i == "7-octenoic acid":
					if i == j:
						intersection.append(i)
						continue
					continue
				intersection.append(i)

				try:
					qu2_to_mobp2[i].append(j)
				except KeyError:
					qu2_to_mobp2[i] = [j]

	intersection = list(set(qu2) & set(mobp2))

	# then accounting for ligands differently spelt

	intersection = intersection + ["6-methyl-5-hepten-2-one","(-)-citronellal","geranylacetone","(e)-2-hexenal","n,n-diethyl-m-toluamide (deet)"]

	quSpelling_to_mobpSpelling= \
	{'(-)-citronellal': 'citronellal',\
	 '(e)-2-hexenal': 'e2-hexenal',\
	 '6-methyl-5-hepten-2-one': '6-methylhept-5-en-2-one',\
	 'geranylacetone': 'geranyl acetone',\
	 "n,n-diethyl-m-toluamide (deet)": "deet"}

	intersection_mobp = []

	for ligand in intersection:
		try:
			ligand_mobp = quSpelling_to_mobpSpelling[ligand]
		except KeyError:
			ligand_mobp = ligand
		intersection_mobp.append(ligand_mobp)

	# ['nonanal',
	#  'indole',
	#  '1-dodecanol',
	#  '7-octenoic acid',
	#  'octanoic acid',
	#  'heptanoic acid',
	#  'nonanoic acid',
	#  'decanal',
	#  'hexanal',
	#  '1-octen-3-ol',
	#  'icaridin',
	#  'benzaldehyde',
	#  'pentanoic acid',
	#  'octanal',
	#  '6-methylhept-5-en-2-one',
	#  'citronellal',
	#  'geranyl acetone',
	#  'e2-hexenal',
	#  'deet']

	# Maybe add linanool-vs-linanool oxide? (qu and mobp respectively)

	return intersection_mobp, quSpelling_to_mobpSpelling 

def plot_venn_ligands(output_path, qu2, mobp2, intersection_mobp):

	""" Plots venn diagram of mOBP vs. Qaio et al. ligands tested in their binding-assays

	Args:

		output_path = "../../data/obp_output/venn_mobp_vs_qaio.png"

		qu2 = see: @lower

		mobp2 = see: @lower

		intersection = see: 

	Usage: 

		plot_venn_ligands("../../data/obp_output/venn_mobp_vs_qaio.png", qu2, mobp2, intersection)

	"""

	plt.figure(figsize=(4,4))
	#v = venn3(subsets=(1, 1, 1, 1, 1, 1, 1), set_labels = ('A', 'B', 'C'))

	v = venn2(subsets = (len(qu2), len(mobp2), len(intersection_mobp)))

	v.get_label_by_id('A').set_text('Exp')
	v.get_label_by_id('B').set_text('MD')


	plt.savefig(output_path)

	plt.close()

########
# Main #
########

intersection_mobp, quSpelling_to_mobpSpelling = generate_ligands_intersection()

#########################
# Plotting venn daigram #
#########################

plot_venn_ligands("../../data/obp_output/venn_mobp_vs_qaio.png", qu2, mobp2, intersection_mobp)

#######################
# Plotting lader plot # @latest
#######################

#
# Generate ligand preference sequences
#

data_path = "/home/qiime/Desktop/hpcleap_ppi/data/obp_input/"

mobps = ['mOBP3.txt','mOBP4.txt','mOBP12.txt','mOBP19.txt']
qobps = ['qOBP3.txt','qOBP4.txt','qOBP12.txt','qOBP19.txt']
#qobps = ['qOBP3_int.txt','qOBP4_int.txt','qOBP12_int.txt','qOBP19_int.txt']

mobps_paths = [data_path+i for i in mobps]
qobps_paths = [data_path+i for i in qobps]
mq_paths 	= zip(mobps_paths,qobps_paths)

obp_to_comparison = {}

for m_path, q_path in mq_paths:

	obp_id = m_path[51:53].replace(".","")

	mobpX_ligands, qobpX_ligands, qu_raw = compare_ligand_sequence(m_path,q_path,intersection_mobp,m_sort_by="sile")

	obp_to_comparison["OBP"+obp_id] = {"mOBP":mobpX_ligands, "qOBP":qobpX_ligands}

#
# Generate ligand preference sequences
#

for obp_id in obp_to_comparison.keys():

	ladder_plots(obp_id,obp_to_comparison,"../../data/obp_output/"+obp_id+"_ladder.png")



