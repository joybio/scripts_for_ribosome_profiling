#!/bin/bash

ls *_R1.fastq.gz | while read id;do(mv $id $(basename $id '_R1.fastq.gz').R1.fastq.gz);done
ls *_R2.fastq.gz | while read id;do(mv $id $(basename $id '_R2.fastq.gz').R2.fastq.gz);done
ls *R1.fastq.gz | while read id;do(echo $(basename $id '.R1.fastq.gz') >> index.csv);done
