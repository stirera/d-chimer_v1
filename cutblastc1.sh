#!/usr/bin/bash
outname=$(basename $1 .csv)
cut -f1 $1 |sort -u > ${outname}.list