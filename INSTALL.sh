#!/usr/bin/bash

pip install biopython
pip install PyYAML
wget "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz"
tar xzf new_taxdump.tar.gz
cat fullnamelineage.dmp|sed "s/\s\+\|\s\+//g"|sed 's/\|//2g'|awk 'BEGIN{FS="|"}{print $1,$3,$4,$2}'|sed "s/\s/\t/g"|sort -k1,1 > fullnamelineage_taxid_sorted.dmp
