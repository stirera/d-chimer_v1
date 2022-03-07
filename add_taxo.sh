#!/usr/bin/bash
sharedbase=$(basename $1 .tsv)
tax_lineages=$2
blasttag=$3

sortedtsv=${sharedbase}.sorted.tsv



#ls -l $tax_lineages

sort -k2,2 $1|sed "s/\s\+/\t/g" > $sortedtsv
# one have to create a taxid sorted "fullname lineage" file form https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
#join -t$'\t' -1 2 -2 1 $sortedtsv ~/work_directory/db/fullnamelineage_taxid_sorted.dmp|cut -f1- > ${sharedbase}.taxo

join -t$'\t' -1 2 -2 1 $sortedtsv $tax_lineages|cut -f1- > ${sharedbase}.taxo

a=$(wc -l $1|cut -d " " -f1)
b=$(wc -l ${sharedbase}.taxo|cut -d " " -f1)


#the input file contained 69250 k_cu_fp.1.bh.filtered.sorted.tsv lines; and the output 69247 k_cu_fp.1.bh.filtered.sorted.taxo lines

echo "the input file ($1) contained $a lines; and the output file (${sharedbase}.taxo) $b lines"

echo "searching for unfound lines in the taxonomic database, please check, ${sharedbase}.notfound_taxo.tsv"

join -t$'\t' -v 1 -1 2 -2 1 $sortedtsv $tax_lineages|cut -f1- > ${sharedbase}.notfound_taxo.tsv
