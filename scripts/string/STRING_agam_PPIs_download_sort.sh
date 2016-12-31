
# Download PPIs from STRING, pertaining to Tax id: 7165 (Anopheles gambiae), then sorts by total score outputting to file:  

cd ../data/string_ppis

wget http://string-db.org/newstring_download/protein.links.detailed.v10/7165.protein.links.detailed.v10.txt.gz

sort -t' ' -k 10,10 -nr ../data/string_ppis/7165.protein.links.detailed.v10.txt > ../data/string_ppis/7165.protein.links.detailed.v10_sorted.txt