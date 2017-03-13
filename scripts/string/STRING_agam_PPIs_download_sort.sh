
#
# Download PPIs from STRING, pertaining to Tax id: 7165 (Anopheles gambiae), then sorts by total score outputting to file:  

# USAGE: bash ./(must be in same dir as ./STRING*.sh)


#
# Create data directory if missing...
#

STRINGDIR="../../data/string_ppis";

echo "$STRINGDIR";

if [ ! -d "$STRINGDIR" ]; then
	echo "Oh dear, missing: '../../data/string_ppis', creating new...";
	mkdir "$STRINGDIR";
fi

#
# CD data dir, and download STRING file...
#


echo "Entering the STRING data directory, to start downloading...";

cd "$STRINGDIR";


# @TODO: uncomment after testing ^

#wget http://string-db.org/newstring_download/protein.links.detailed.v10/7165.protein.links.detailed.v10.txt.gz;

echo "Sorting PPIs according to combined score rank...";

sort -t' ' -k 10,10 -nr ./7165.protein.links.detailed.v10.txt > ./7165.protein.links.detailed.v10_sorted.txt;
