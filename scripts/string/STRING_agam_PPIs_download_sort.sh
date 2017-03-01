
cd ~/Desktop/

wget http://string-db.org/newstring_download/protein.links.detailed.v10/7165.protein.links.detailed.v10.txt.gz

sort -t' ' -k 10,10 -nr ./7165.protein.links.detailed.v10.txt > 7165.protein.links.detailed.v10_sorted.txt