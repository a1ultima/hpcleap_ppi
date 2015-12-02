
from operator import itemgetter

data = []

fi = open("../data/string_ppis/7165.protein.links.detailed.v10.txt","r")

while True: 

	line = fi.readline()

	if line=="":
		break

	if "protein1" in line:
		headers = line 
		continue

	line_split = line.split(" ")

	# protein1		= line_split[0]
	# protein2  		= line_split[1]
	# neighborhood 	= line_split[2]
	# fusion 			= line_split[3]
	# cooccurence 	= line_split[4]
	# coexpression 	= line_split[5]
	# experimental 	= line_split[6]
	# database 		= line_split[7]
	# textmining 		= line_split[8]
	# combined_score 	= line_split[9]

	data.append(line_split)

fi.close()


data = sorted(data, key=itemgetter(9))




