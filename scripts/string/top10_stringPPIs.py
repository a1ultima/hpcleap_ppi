
#
# Imports
#

from operator import itemgetter
import numpy as np
import pdb

#
# Functions
#

def string_parser( fi_name ):
    """ returns a list containing the data-rows, header line excluded, of STRING PPIs outputted by: ./STRING_agam_PPIs_download_sort.sh, e.g. "../../data/string_ppis/7165.protein.links.detailed.v10_sorted.txt", hereby: @@BASHDOWNLOADER

    ARGS:
        fi_name (string):   pathname to the output file of @BASHDOWNLOADER, e.g. "../../data/string_ppis/7165.protein.links.detailed.v10_sorted.txt"

    USAGE:
        data = string_parser( "../../data/string_ppis/7165.protein.links.detailed.v10_sorted.txt" )

    RETURNS:
        data (list):    list of all @data_rows (excl. header line) of the output file of @BASHDOWNLOADER, hereby: @@data

    """

    fi = open( fi_name, "r") # @@fi_name (formerly @@_fi), e.g. "../../data/string_ppis/7165.protein.links.detailed.v10_sorted.txt"

    while True: 

        line = fi.readline()

        if line=="":
            break

        if "protein1" in line:
            #headers = line # @ANDY:headers may be useful in future // splices away the header line, thus #data_rows = #file_rows - 1 from @fi
            continue

        line_split = line.rstrip().split(" ") # e.g. ['7165.AGAP028012-PA','7165.AGAP009539-PA','594','0','0','994','978','805','865','999']

        #
        # Column headers: (i.e. header line/row)
        #
        # protein1      = line_split[0]
        # protein2      = line_split[1]
        # neighborhood  = line_split[2]
        # fusion        = line_split[3]
        # cooccurence   = line_split[4]
        # coexpression  = line_split[5]
        # experimental  = line_split[6]
        # database      = line_split[7]
        # textmining    = line_split[8]
        # combined_score= float(line_split[9]) # @andy:converting combined_socre field/column of the data file to floats, prior to formatting into np.array() type. For smoother calculations downstream.
        
        #line_split.pop() # throw away the string format "@combined_score" value... 
        #line_split.append(combined_score) # ...to replace it with the float() version, @DONE:test: histogram looks the same // if the output to data is consistent with the previous version just apending line_split without flaot-converted "combined_score" field

        data.append(line_split)

    fi.close()

    return data 

def gen_human_readable_ppi_summary(ntotal_ppis):
    """ returns the number of rows parsed in by @data <-- @fi (../../data/string_ppis/7165.protein.links.detailed.v10_sorted.txt), as a human readable string of three consecutive digits delimited by commas 

    ARGS: 
        ntotal_ppis (int): number of rows of @data, i.e. the number of STRING PPIs available (@fi)

    RETURNS:
        str @todo

    @@@TESTING: 
        checked number of @@lines in @BSAHDOWNLOADER file was == number of data_rows + 1 (1 being the header-line of the file (i.e. not a PPI data row)) hereby @@data_rows

            wc -l < 7165.protein.links.detailed.v10_sorted.txt (returns 1958813) == 1 + ntotal_ppis (i.e. @lines) (returns 1958812))
    """

    #
    # Sub-functions
    #
    def chunks(l, n):
        """ Yield successive n-sized chunks from list l. """
        for i in xrange(0, len(l), n):
            yield l[i:i+n]

    #
    # Generate human readable #total ppis found (data rows, excl. headers)
    #
    ntotoal_ppis_list    = list(str(ntotal_ppis)) # e.g. ['1', '9', '5', '8', '8', '1', '2']
    ntotal_ppis_reversed = list(reversed(ntotoal_ppis_list)) # e.g. ['2', '1', '8', '8', '5', '9', '1']
    ntotal_ppis_chunked  = list(chunks(ntotal_ppis_reversed,3)) # e.g. [['2', '1', '8'], ['8', '5', '9'], ['1']]
    ntotal_ppis_readable_list = list(reversed([list(reversed(i)) for i in ntotal_ppis_chunked])) # e.g. [['1'], ['9', '5', '8'], ['8', '1', '2']]
    ntotal_ppis_readable = ",".join(["".join(i) for i in ntotal_ppis_readable_list]) # e.g. '1,958,812' <-- ['1', '958', '812']
    return ntotal_ppis_readable



#
# Instantiation
#

data = []

prev_sorted = True

scoring_threshold_type = "percentile" #, hereby @REL, OR, e.g. scoring_threshold_type =  "absolute" @ABS

# allternative 1. {{ 

#ppi_acceptance_threshold = 9000     # filter PPIs, keeping those classed as "@likely PPIs", by means of an absolute user-specified @scoring_threshold, and @scoring_threshold_type (="absolute"), hereby @@ABS

# }} 1 alternative 2. {{
    
# }}

scoring_threshold = 99.5 # filter PPIs, keeping those classed as "@likely PPIs", by means of a statistical user-specified @scoring_threshold (cutoff). User-specified: N-th percentile cutoff, corresponding to the user-defined scoring_colum, such that PPIs with scores below the percentile_cutoff value are classed as "non-interactions (non-PPIs)", or "unlikely interactions (unlikely-PPIs)", and PPIs whose scores are above percentile_cutoff are classed as, "likley interactions (likely-PPIs)"

#
# Parse STRING PPI data
#

print "Reading score-sorted PPI (STRING) input file, please be patient..."

data = string_parser( "../../data/string_ppis/7165.protein.links.detailed.v10_sorted.txt" )

#
# Sort by column: "combined score", if not previously sorted by ./STRING_agam_PPIs_download_sort.sh @bashdownloader
#
if prev_sorted == False:
    data = sorted(data, key=itemgetter(9)) # @@data

#
# Completion summary for @string_parser
#

ntotal_ppis_in = len(data) # e.g. 1958812
ntotal_ppis_readable_out = gen_human_readable_ppi_summary(ntotal_ppis_in)

print "\tComplete! ... Total umber of STRING PPIs parsed: "+ntotal_ppis_readable_out

#
# Return top N-scoring PPIs, specified by: ppi_acceptance_threshold
#

print "Determining most likely PPI candidates..."

data_arr = np.array(data)

combined_score_col = data_arr[1:,-1].astype(np.float)

if scoring_threshold_type == "percentile":

    print "\t user-specified scoring_threshold_type: "+scoring_threshold_type
    print "\t\tremoving PPIs scoring < "+str(scoring_threshold)+"-th percentile..."
    # p = np.percentile(combined_score_col, 50) # e.g. np.percentile(a, 50), will return 50th percentile, e.g median, of a numpy array of floats
    p_cutoff = np.percentile(combined_score_col, scoring_threshold) # e.g. np.percentile(a, scoring_threshold=50), will return 50th percentile, e.g median, of a numpy array of floats
    # @TODO:test:is using the np.array faster than using the plain list()? 
    # throw away all PPIs with a score < @p_cutoff

    likely_ppis = [(data_row[0],data_row[1]) for data_row in data if (float(data_row[9])>p_cutoff)] # e.g. 0, 1 and 9: 1st and 2nd columns, and 10th columns corresponding to "protein1", "protein2", and "combined_score" fields of @BASHDOWNLOADER output file; e.g. data_row = ['7165.AGAP028012-PA','7165.AGAP009539-PA','594','0','0','994','978','805','865','999']

    n_likely_ppis = len(likely_ppis)

    n_likely_ppis_readable = gen_human_readable_ppi_summary(n_likely_ppis)

    print "\t\t\tNo. of most likely PPIs detected: "+n_likely_ppis_readable+" ("+str((float(len(likely_ppis))/ntotal_ppis_in)*100)+"%"+" of total input PPIs!)"

elif scoring_threshold_type == "absolute":

    print "\t user-specified scoring_threshold_type: "+scoring_threshold_type
    print "\t\tremoving PPIs scoring < "+str(scoring_threshold)+"..."

    likely_ppis = [(data_row[0],data_row[1],data_row[-1]) for data_row in data if (float(data_row[9])>scoring_threshold)] # e.g. 0, 1 and 9: 1st and 2nd columns, and 10th columns corresponding to "protein1", "protein2", and "combined_score" fields of @BASHDOWNLOADER output file; e.g. data_row = ['7165.AGAP028012-PA','7165.AGAP009539-PA','594','0','0','994','978','805','865','999']



# #
# # Plot histogram of scores
# #

# import matplotlib.pyplot as plt

# # gaussian_numbers = np.random.randn(1000) # @toy-example, gaussian_numers --> combined_score_col
# plt.hist(combined_score_col)
# plt.title("PPI Scores Histogram")
# plt.xlabel("Value")
# plt.ylabel("Frequency")
# plt.show()


