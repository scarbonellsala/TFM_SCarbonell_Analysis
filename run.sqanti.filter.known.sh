#!/bin/bash -ex

# get list of isoforms from SQANTI evalation that are considered
# "known" relative to gencode annotation
# then use 
#   gffread --ids known.trans.tsv --gtf flair.collapse.isoforms.gtf > known_isoforms.gtf
# to get a GFF of know isoforms


# flair.collapse.isoforms_classification.txt columns
#  1	isoform
#  6	structural_category

awk -v 'FS=\t' -v 'OFS=\t'  'BEGIN{print "isoform"}
      $6~/full-splice_match|incomplete-splice_match|novel_in_catalog|novel_not_in_catalog/{print $1}' \
          sqanti-all/flair.collapse.isoforms_classification.txt > known.trans.tsv


