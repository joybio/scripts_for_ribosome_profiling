#!/bin/bash

ls *.unaligned.sam | while read id;do(awk '{print ">"$1"\n"$10}' $id > $(basename $id '.unaligned.sam').unaligned.fasta);done



