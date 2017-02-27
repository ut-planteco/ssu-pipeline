# ssu-pipeline

Series of scripts to analyse 454 and Illumina sequences in SSU amplicon against [Maarj_AM_](http://maarjam.botany.ut.ee) database. This is pipeline that has been developed and used in paper "TODO" (Vasar et al xxxx). Both 454 and Illumina pipelines are optimized on analyzing SSU sequences against MaarjAM database using BLAST based OTU picking approach. For 454, we assume that sequences are multiplexed into one file and for Illumina sequences are demultiplexed into separate fastq or fastq.tar.gz (packed) pairs for each sample as it is general approach. 

Prerequisite software:

* [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) 
* [USEARCH](http://www.drive5.com/usearch/)
* [FLASh](https://ccb.jhu.edu/software/FLASH/)
* [Python 2.7](https://www.python.org/)

## 1.1. Clean raw 454 sequences

Clean 454 sequences by preparing barcode and primer list for demultiplexing. Example file layout provided as follow, where columns are seperated with tabular and for each column sample, barcode and primer are shown:

```
Sample1    ACGTACGT    TCCGTCCGTTTGGTTA
Sample2    ACGTATCT    TCCGTCCGTTTGGTTA
.......    ........    ................
```

Sample | Barcode | Primer
--- | --- | ---
Sample1 | ACGTACGT | TCCGTCCGTTTGGTTA
Sample2 | ACGTATCT | TCCGTCCGTTTGGTTA

Run python from command line as following:
```
python pipeline_clean_454.py -f input.fasta -qf input.qual -b input.barcode -bs 1 -bb 2 -bp 3 -q 25 -trimq 20 -trimw 50 -ml 170 -tl 520
```

Command help:
```
python pipeline_clean_454.py
    -f FASTA_FILE [-qf QUALITY_FILE] -b
    BARCODE_FILE -bs SAMPLE_COLUMN -bb BARCODE_COLUMN
    [-bp PRIMER_COLUMN] [-q AVERAGE_QUALITY]
    [-trimq TRIM_QUALITY] [-trimw TRIM_WINDOW]
    [-ml MINIMUM_LENGTH] [-tl TRUNCATE_LENGTH]
```

## 1.2. Clean raw Illumina sequences

Clean Illumina sequences by defining the folder where paired reads are located and provide forward and reverse for both primers with average quality. Make sure that demultiplexed file names coming from Illumina MiSeq platform are correct. Script will gather files named as `SAMPLE_R1_001.fastq` or `SAMPLE_R1_001.fastq.tar.gz`. Script will interleave correct forward and reverse reads together that can be easily used by FLASh software to pair them. Because FLASh output FASTQ, we need to convert it to FASTA to make it understandable for BLAST. In order to skip intermediate files, pipe each step into one command as following:
```
python pipeline_clean_illumina.py [-h] -folder FOLDER -fprimer TCCTTCCAA [-rprimer TCCTTCCAA] | flash | python pipeline_fastq_fasta.py -l 400 > cleaned.fasta
```
## 2. Remove chimeric sequences

Once files are cleaned, we need to remove chimeric sequences that are introduced using PCR. We use USEARCH in reference database mode against Maarj_AM_ database. 
usearch 

## 3. Identify reads against reference database

Once we have removed chimeric reads, we can start identifying sequences using BLAST+ software and Maarj_AM_ database. 
blastn | python pipeline_parseblast.py

## 4. Summarize BLAST results

Finally, we can summarize BLAST result using parsed output. Providing FASTA file will output also nohit selection that can be used for further BLAST against additional databases.

python pipeline_summarize.py

## 5. Final results

Two files are generated from previous step. One is named as `*.nohits.fasta`, where sequences that did not get significant hit against reference database are written out and `*.tsv`, where results in pivot table form are written and hits are sorted in descending order.

___

Use following citation when using our python scripts:
```
Vasar M, Öpik M, ... (xxxx) Cleaning 454 and MiSeq Illumina, Paper 5: 336-352
```

License: [CC-BY](https://creativecommons.org/licenses/by/3.0/)