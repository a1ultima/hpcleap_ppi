
#
# Imports
#

from operator import itemgetter
import numpy as np
import matplotlib.pyplot as plt
import sys 
import networkx as nx

import time

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
            headers = line.rstrip().split(" ") # @ANDY:headers may be useful in future // splices away the header line, thus #data_rows = #file_rows - 1 from @fi
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
        # combined_score= float(line_split[9]) # @andy:converting combined_score field/column of the data file to floats, prior to formatting into np.array() type. For smoother calculations downstream.
        
        #line_split.pop() # throw away the string format "@combined_score" value... 
        #line_split.append(combined_score) # ...to replace it with the float() version, @DONE:test: histogram looks the same // if the output to data is consistent with the previous version just apending line_split without flaot-converted "combined_score" field

        data.append(line_split)

    fi.close()

    return data, headers 

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


def plot_histogram( fo_name_suffix, likely_ppis, combined_score_col_all ):

    """ Plots and saves histogram comparing the score distributions between "All PPIs" vs. "Likely PPIs", saving to filepath specified by @hist_fo_path, a modification of @fo_name_suffix

        ARGS:

            fo_name_suffix (str):   @BASHDOWNLOAD outputfile name, recycled to name the histogram plot, for consistency. e.g. ['../../data/string_ppis/7165.protein.links.detailed.v10_sorted',[]], hereby: @@fo_name_suffix
            likely_ppis (list):     list of rows, each row consisting of three column elements (protein1, protein2, combined_score), see: @TODO
            combined_score_col_all (numpy_array): vector containing jus the combined_scores of each PPI, comprising the unfiltered "All PPIs" class, see: @TODO

        RETURNS:
            None

        USAGE:
            plot_histogram( fo_name_suffix, likely_ppis, combined_score_col_all )


    """

    # extending name of histogram plot, to make it clear it is not the output .txt file
    hist_fo_name_suffix     = fo_name_suffix
    hist_fo_name_suffix[1]  = "_histogram_likely_vs_all_PPIs.png"
    hist_fo_path            = "".join(hist_fo_name_suffix)

    # histogram configuration
    bins = np.linspace(100, 1000, 25) # start bin, end bin, total number of bins

    #
    # Frequency distribution of scores, compared between "Likely PPIs" (after thresholding) vs. "All PPIs" class sets, hereby:@@FREQ 
    #
    likely_ppis_arr           = np.array(likely_ppis)
    combined_score_col_likely = likely_ppis_arr[1:,2].astype(np.float) 

    #
    # Comparative histograms, pertaining to: @FREQ comparisons of scores between "likely PPIs" vs. "All PPIs"
    #
    plt.hist(combined_score_col_all, bins, alpha=0.5, label='All PPIs')
    plt.hist(combined_score_col_likely, bins, alpha=0.5, label='Likely PPIs')

    print "\tWriting histogram (.png) file to: \n\t\t"+hist_fo_path

    plt.legend(loc='upper right')
    plt.title("PPI Scores Histogram")
    plt.xlabel("STRING PPI Score")
    plt.ylabel("Frequency")
    #plt.show()
    plt.savefig(hist_fo_path)
    plt.close()

    #return combined_score_col_likely


def write_likely_PPIs_to_file( fo_name_suffix, likely_ppis ):

    """ Writes minimal formatted, PPI data after filtering via scoring_threshold, to file.

    ARGS:
        fo_name_suffix (str):   @BASHDOWNLOAD outputfile name, recycled to name the histogram plot, for consistency. e.g. ['../../data/string_ppis/7165.protein.links.detailed.v10_sorted',[]], hereby: @@fo_name_suffix
        likely_ppis (list):     list of rows, each row consisting of three column elements (protein1, protein2, combined_score), see: @TODO

    RETURNS:
        fo_path (str):  modified output file name.

    USAGE:
        write_likely_PPIs_to_file( fo_name_suffix, likely_ppis )        

    """

    # extend filename consistenly named to @BASHLOADER's output file name, with extension indicating "Likely PPIs" (after filering) data
    fo_name_suffix[1] = "_likely-PPIs_cutoff_"+str(scoring_threshold)+"_"+scoring_threshold_type
    fo_path = "".join(fo_name_suffix)

    print "Writing 'Likely PPIs' (filtered) data file to...\n\t"+fo_path+"\n"

    # write to file
    with open(fo_path,'w') as fo:
        fo.write("protein1\tprotein2\tcombined_score\n") # column headers
        for ppi in likely_ppis:
            fo.write("\t".join(ppi)+"\n") # data_rows

    return fo_path


def draw_graph( fo_name_suffix, likely_ppis ): 

    """ Plots a network representation of PPIs given a list, in the format of: see: "likley_ppis"
        
        ARGS:
            fo_name_suffix (str):   @BASHDOWNLOAD outputfile name, recycled to name the histogram plot, for consistency. e.g. ['../../data/string_ppis/7165.protein.links.detailed.v10_sorted',[]], hereby: @@fo_name_suffix
            likely_ppis (list):     list of rows, each row consisting of three column elements (protein1, protein2, combined_score), see: @TODO

        RETURNS:
            None

        USAGE:
            draw_graph( fo_name_suffix, likely_ppis )            

    """

    # extend filename consistenly named to @BASHLOADER's output file name, with extension indicating "Likely PPIs" (after filering) data.
    graph_fo_name_suffix = fo_name_suffix
    graph_fo_name_suffix[1] = "_likely-PPIs_network_"+str(scoring_threshold)+"_"+scoring_threshold_type+".png"
    graph_fo_path = "".join(graph_fo_name_suffix)

    print "Writing graph/network plot file to:\n\t"+graph_fo_path+"\n"

    # NetworkX graph object
    ppi_graph = nx.Graph()

    for ppi in likely_ppis:

        p1 = ppi[0] # e.g. '7165.AGAP011298-PA'
        p2 = ppi[1] # e.g. '7165.AGAP010252-PA'
        #p1p2_score = ppi[2] # e.g. '999'

        ppi_graph.add_node(p1)
        ppi_graph.add_node(p1)
        ppi_graph.add_edge(p1,p2)

    nx.draw(ppi_graph,pos=nx.spring_layout(ppi_graph))
    #plt.show()
    plt.title("Spring-board graph representation of threshdold-filtered STRING PPIs")
    plt.savefig(graph_fo_path)



#
# @TESTING:timing efficiency alternatives
#
start_time = time.time()  # @time


#
# Instantiation, configuration
#

data = []

prev_sorted = True

bashDownloader_fo_path = "../../data/string_ppis/7165.protein.links.detailed.v10_sorted.txt"

scoring_threshold_type = "percentile" #, hereby @REL, OR, e.g. scoring_threshold_type = "absolute" @@ABS

# allternative 1. {{  Absolute-score-filter: removes PPIs scoring less than the specified absolute: scoring_threshold, @ABS

# scoring_threshold = 900 # 608 median score, for e.g. for @BASHDOWNLOADER output file: "../../data/string_ppis/7165.protein.links.detailed.v10_sorted.txt", note that 7165 is the tax id for: Anopheles gambiae (as of 2016) # filter PPIs, keeping those classed as "@likely PPIs", by means of an absolute user-specified @scoring_threshold, and @scoring_threshold_type (="absolute"), hereby @@ABS

# }} 1 alternative 2. {{ Percentile-score-fileter: removes PPIs scoring less than the specific relative (percentile): scoring_threshold (0.00-to-100.00), where 50.00 is the median

scoring_threshold = 99.5 # e.g. with 99.5, 8722 PPIs are kept (out of 1,958,812); i.e. 99.5-th percentile is used to keep only PPIs scoring higher than the 99.5th score percentile, for "combined_scores" column  
#scoring_threshold = 99.5 # filter PPIs, keeping those classed as "@likely PPIs", by means of a statistical user-specified @scoring_threshold (cutoff). User-specified: N-th percentile cutoff, corresponding to the user-defined scoring_colum, such that PPIs with scores below the percentile_cutoff value are classed as "non-interactions (non-PPIs)", or "unlikely interactions (unlikely-PPIs)", and PPIs whose scores are above percentile_cutoff are classed as, "likley interactions (likely-PPIs)"
#@TODO:change, to allow colouring of edges of the STRING network viz. script-based on differnt evidence-type

# }} 2 alternative.


#
# Parse STRING PPI data
#

print "Reading score-sorted PPI (STRING) input file, please be patient..."

data, headers = string_parser( bashDownloader_fo_path )

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

print "\tTotal number of STRING PPIs parsed (before filtering): "+ntotal_ppis_readable_out

#
# Return top N-scoring PPIs, specified by: ppi_acceptance_threshold
#

# @TODO:turn into function {{ ...

print "Determining most 'Likely PPI' candidates (score threshold filtering)..."

data_arr = np.array(data)

combined_score_col_all = data_arr[1:,-1].astype(np.float)

if scoring_threshold_type == "percentile":

    print "\t user-specified scoring_threshold_type: "+scoring_threshold_type
    print "\t\tremoving PPIs scoring < "+str(scoring_threshold)+"-th percentile..."
    # p = np.percentile(combined_score_col_all, 50) # e.g. np.percentile(a, 50), will return 50th percentile, e.g median, of a numpy array of floats
    p_cutoff = np.percentile(combined_score_col_all, scoring_threshold) # e.g. np.percentile(a, scoring_threshold=50), will return 50th percentile, e.g median, of a numpy array of floats
    # @TODO:test:is using the np.array faster than using the plain list()? 
    # throw away all PPIs with a score < @p_cutoff

    # @TODO:code-redunancy    
    # alternative 1. {{ Traverse non-numpy vanilla python list (Slower?) @TODO test <--

    likely_ppis = [[data_row[0],data_row[1],data_row[-1]] for data_row in data if (float(data_row[9])>p_cutoff)] # e.g. 0, 1 and 9: 1st and 2nd columns, and 10th columns corresponding to "protein1", "protein2", and "combined_score" fields of @BASHDOWNLOADER output file; e.g. data_row = ['7165.AGAP028012-PA','7165.AGAP009539-PA','594','0','0','994','978','805','865','999']

    # @TODO:complete file writing and threshold-filtering steps in the same iteration (instead of separateley later)

    # }} 1 alternative 2. {{ Traverse numpy array (Faster?) @TODO: test <---

    #likely_ppis = [[data_row[0],data_row[1],data_row[-1]] for data_row in data_arr if (float(data_row[9])>p_cutoff)] # e.g. 0, 1 and 9: 1st and 2nd columns, and 10th columns corresponding to "protein1", "protein2", and "combined_score" fields of @BASHDOWNLOADER output file; e.g. data_row = ['7165.AGAP028012-PA','7165.AGAP009539-PA','594','0','0','994','978','805','865','999']

    # }} 2. alternative

elif scoring_threshold_type == "absolute":

    print "\t user-specified scoring_threshold_type: "+scoring_threshold_type
    print "\t\tremoving PPIs scoring < "+str(scoring_threshold)+"..."

    # alternative 1. {{ Traverse non-numpy vanilla python list (Slower?) @TODO test <--

    likely_ppis = [[data_row[0],data_row[1],data_row[-1]] for data_row in data if (float(data_row[9])>scoring_threshold)] # e.g. 0, 1 and 9: 1st and 2nd columns, and 10th columns corresponding to "protein1", "protein2", and "combined_score" fields of @BASHDOWNLOADER output file; e.g. data_row = ['7165.AGAP028012-PA','7165.AGAP009539-PA','594','0','0','994','978','805','865','999']

    # }} 1 alternative 2. {{ Traverse numpy array (Faster?) @TODO: test <---

    # likely_ppis = [[data_row[0],data_row[1],data_row[-1]] for data_row in data_arr if (float(data_row[9])>scoring_threshold)] # e.g. 0, 1 and 9: 1st and 2nd columns, and 10th columns corresponding to "protein1", "protein2", and "combined_score" fields of @BASHDOWNLOADER output file; e.g. data_row = ['7165.AGAP028012-PA','7165.AGAP009539-PA','594','0','0','994','978','805','865','999']

    # }} 2. alternative

# @DONE:code redundancy // factored out
n_likely_ppis = len(likely_ppis)
n_likely_ppis_readable = gen_human_readable_ppi_summary(n_likely_ppis)
print "\t\t\tNo. of most 'Likely PPIs' detected: "+n_likely_ppis_readable+" ("+str((float(len(likely_ppis))/ntotal_ppis_in)*100)+"%"+" of total input PPIs)"

#  }} ... @TODO:turn into function 

### 
### Consistent handle for downstream output file names, to @BASHLOADER's outputfile names
### 
fo_name_suffix = bashDownloader_fo_path.split(".txt")

###
### Plot histogram of score distributions (likely vs. all) + network visualisation of likely PPIs
###
print "Comparing score distributions between 'likely PPIs' (filtered) vs. 'All PPIs' (raw input) class sets..."
plot_histogram( fo_name_suffix, likely_ppis, combined_score_col_all )
print "Generating netowrk visualisation, of all 'likely PPIs' interactions..."
draw_graph( fo_name_suffix, likely_ppis )  

### 
### Writing output to file
### 
fo_name_likely_PPIs = write_likely_PPIs_to_file( fo_name_suffix, likely_ppis )        

###
### End
###
print "All jobs complete!"

#
# @TESTING: time efficiency alternatives
#
print "Elapsed time: "+str(time.time() - start_time)   # @time