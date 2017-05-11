# ssu-pipeline

Series of scripts to analyse 454 and Illumina sequences in SSU amplicon against [Maarj*AM*](http://maarjam.botany.ut.ee) database. This is pipeline that has been developed and used in paper *"Comparison of 454 and Illumina sequencing methods to study arbuscular mycorrhizal fungal community diversity"* (Vasar et al. xxxx). Both 454 and Illumina pipelines are optimized on analyzing SSU sequences against MaarjAM database using BLAST based OTU picking approach. For 454, we assume that sequences are multiplexed into one file and for Illumina sequences are demultiplexed into separate fastq or fastq.tar.gz (packed) pairs for each sample as it is general approach. 

## Prerequisite software

* [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) 
* [USEARCH](http://www.drive5.com/usearch/)
* [FLASh](https://ccb.jhu.edu/software/FLASH/)
* [Python 2.7](https://www.python.org/)

## Install software and clone project

First download ssu-pipeline

```
git clone https://github.com/ut-planteco/ssu-pipeline
```

Update repository and install mandatory software
```
sudo apt-get update
sudo apt-get install ncbi-blast+
wget https://sourceforge.net/projects/flashpage/files/FLASH-1.2.11.tar.gz/download
```

You need to download and install [USEARCH](http://www.drive5.com/usearch/download.html) from their homepage as they provide download link and license per e-mail.

Example datasets for 454 and Illumina are located in `example_data.tar.gz` and unpacking will create folders `454` and `illumina`. All the example commands are using these files. Unpack dataset with following command

```
tar -xzvf example_data.tar.gz
```

Inside `maarjam` folder is located Maarj*AM* database (status October 2016) with FASTA and BLAST+ formatted file formats that can be directly used to identify sequences. Use BLAST+ formatted files as it will allow to use multiple cores compared to only using FASTA file. You can format FASTA file into BLAST+ format as following

```
makeblastdb -in reference.fasta -dbtype nucl -title CustomDB -out reference
```

## 1.1. Clean raw 454 sequences

Clean 454 sequences by preparing barcode and primer list for demultiplexing. Example file layout provided as following, where columns are seperated with tabular and for each column sample, barcode and primer:

```
JS1RS5     AGTGAGTG    TTGGAGGGCAAGTCTGGTGCC
JS2CM2     AGTGTCTG    TTGGAGGGCAAGTCTGGTGCC
.......    ........    .....................
```

Sample | Barcode | Primer
--- | --- | ---
JS1RS5 | AGTGAGTG | TTGGAGGGCAAGTCTGGTGCC
JS2CM2 | AGTGTCTG | TTGGAGGGCAAGTCTGGTGCC

We need to define for the script our fasta and quality file locations, sample list and how the tab delimited sample list file is formated: in which column sample, barcode and primer are listed. We define average quality `-q` to be at least 25 (0 to 40), minimum length `-ml` after barcode and primer removal should be at least 170nt and longer `-tl` than 520nt sequences are trimmed to shorter to remove reverse primer AML2. SSU amplicon length between primers is approximately 520nt. Run python from command line as following:
```
python pipeline_clean_454.py -f 454/example.fasta -qf 454/example.qual -b 454/example.barcode -bs 1 -bb 2 -bp 3 -q 25 -ml 170 -tl 520
```

Command help

```
python pipeline_clean_454.py
    -f FASTA_FILE [-qf QUALITY_FILE] -b
    BARCODE_FILE -bs SAMPLE_COLUMN -bb BARCODE_COLUMN
    [-bp PRIMER_COLUMN] [-q AVERAGE_QUALITY]
    [-trimq TRIM_QUALITY] [-trimw TRIM_WINDOW]
    [-ml MINIMUM_LENGTH] [-tl TRUNCATE_LENGTH]
arguments:
  -f FASTA_FILE        FASTA file, where sequences are stored
  -qf QUALITY_FILE     QUALITY file, where sequence qualities are stored for
                       each nucleotide. This file headers should match FASTA
                       file headers.
  -b BARCODE_FILE      BARCODE file, where sample ID, barcodes and primers are
                       stored in tabular file format
  -bs SAMPLE_COLUMN    sample column in BARCODE file (first column is 1)
  -bb BARCODE_COLUMN   barcode column in BARCODE file
  -bp PRIMER_COLUMN    primer column in BARCODE file
  -q AVERAGE_QUALITY   lower limit of average quality of the sequence to be
                       filtered out (recommended = 25)
  -trimq TRIM_QUALITY  lower limit of average quality that is trimmed away
                       when it drops below the threshold (recommended = 20)
  -trimw TRIM_WINDOW   window size to calculate average quality for trimming
                       sequence end (recommended = 50)
  -ml MINIMUM_LENGTH   minimum allowed length of the sequences after trimming
                       (170)
  -tl TRUNCATE_LENGTH  truncate sequences longer than provided length to
                       remove reverse primer
```

Output is written into file `454/example.cleaned.fasta`, move it into main folder for simplicity.

```
mv 454/example.cleaned.fasta 454.cleaned.fasta
```

## 1.2. Clean raw Illumina sequences

Clean Illumina sequences by defining the folder `-folder` where paired reads are located and provide forward and reverse for both primers with average quality. Make sure that demultiplexed file names coming from Illumina MiSeq platform are correct. Script will gather files named as `SAMPLE_R1_001.fastq` or `SAMPLE_R1_001.fastq.tar.gz`. Script will interleave correct forward and reverse reads together that can be easily used by FLASh software to pair them. Because FLASh output FASTQ, we need to convert it to FASTA to make it understandable for BLAST. We also define Illumina Nextera adapters first 10 nucleotides to remove sequences containing part of the adapter for forward `-fadapter` and reverse `-radapter` reads. As the example data is using tagmentation based Illumina, we do not need to define forward `-fprimer` and reverse `-rprimer` primers. Finally we define that average quality `-quality` for both reads needs to be at least 30 (0-40). In order to skip intermediate files, pipe each step into one command as following:
```
python pipeline_clean_illumina.py -folder illumina/ -fprimer "" -rprimer "" -fadapter CTGTCTCTTA -radapter CTGTCTCTTA -quality 30 | ~/applications/FLASH/flash -m 10 -M 300 --interleaved-input - -c | python pipeline_fastq_fasta.py > illumina.cleaned.fasta
```

Nextera adapters R1 GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG and R2 TGTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG needs to be reverse complement. We only need to match first ten bases to find adapters. These 10 bases have been checked against Maarj*AM* database and no interference using short adapter sequence is found to catch false positives. 

Command help

```
python pipeline_clean_illumina.py 
    -folder FOLDER [-fprimer SEQUENCE]
    [-rprimer SEQUENCE] [-fadapter SEQUENCE]
    [-radapter SEQUENCE] [-quality QUALITY]
    [-phred QUALITY]
arguments:
  -folder FOLDER      define FOLDER where FASTQ or FASTQ.tar.gz files are
                      stored
  -fprimer SEQUENCE   define forward read primer
  -rprimer SEQUENCE   define reverse read primer
  -fadapter SEQUENCE  define adapter for forward read
  -radapter SEQUENCE  define adapter for reverse read
  -quality QUALITY    average quality of sequence to be accepted
  -phred QUALITY      FASTQ file phred quality score (33)
```

## 1.3 Correct strand of the sequences for Illumina

Tagmentation based Illumina produces sequences that are not in the same direction, but USEARCH software needs to have all the reads in same direction as the reference database, we need to change them into correct strand. All the sequences should start from NS31 primer and end with AML2 primer. To achieve this, we use Maarj*AM* database with our cleaned sequences and run BLAST+ software to identify strand of the sequences. Sequences identified as +/- by the BLAST+ needs to be reverse complemented. 

```
blastn -query illumina.cleaned.fasta -evalue 1e-50 -max_target_seqs 1 -num_threads 4 -db maarjam/maarjam -outfmt 5 | python pipeline_parse_blast.py > illumina.strand.blast
```

Now run python script that reads BLAST results and fasta input to change direction of the sequences

```
python pipeline_correct_direction.py -f illumina.cleaned.fasta -b illumina.strand.blast > illumina.correct.fasta
```

## 2. Remove chimeric sequences

Once files are cleaned, we need to remove chimeric sequences that are introduced using PCR. We use USEARCH in reference database mode against Maarj*AM* database. Make sure to use correct input file for 454 and Illumina, `454.cleaned.fasta` and `illumina.correct.fasta` respectivelly.

```
usearch -uchime_ref 454.cleaned.fasta -db maarjam/maarjam.fasta -nonchimeras 454.cf.fasta -strand plus
usearch -uchime_ref illumina.correct.fasta -db maarjam/maarjam.fasta -nonchimeras illumina.cf.fasta -strand plus
```

## 3. Identify reads against reference database

Once we have removed chimeric reads, we can start identifying sequences using BLAST+ software and Maarj*AM* database. 

```
blastn -query 454.cf.fasta -evalue 1e-50 -max_target_seqs 1 -num_threads 4 -db maarjam/maarjam -outfmt 5 | python pipeline_parse_blast.py > 454.cf.blast
blastn -query illumina.cf.fasta -evalue 1e-50 -max_target_seqs 1 -num_threads 4 -db maarjam/maarjam -outfmt 5 | python pipeline_parse_blast.py > illumina.cf.blast
```

## 4. Summarize BLAST results

Finally, we can summarize BLAST result using parsed output. Providing FASTA file will output also nohit selection that can be used for further BLAST against additional databases. We use parameters `-vs` and `-ve` do define reference database variable region location. Because we use Maarj*AM* database in this example, all the referene sequences start after NS31 primer and variable region on the amplicon is located from 70nt to 300nt after the NS31 primer. We also define hit identity `-i` to be at least 97% and alignment length `-l` for the hit at least 95% to be counted as a hit.

```
python pipeline_summarize_blast.py -f 454.cf.fasta -b 454.cf.blast -i 97 -l 95 -t 0 -vs 70 -ve 300
python pipeline_summarize_blast.py -f illumina.cf.fasta -b illumina.cf.blast -i 97 -l 95 -t 0 -vs 70 -ve 300
```

Command help

```
python pipeline_summarize_blast.py
    -b BLAST_FILE [-f FASTA_FILE] -i
    IDENTITY[0-100] -l ALIGNMENT[0-100]
    [-vs VARIABLE_START] [-ve VARIABLE_END] -t
    BLAST_TYPE[0-2]
arguments:
  -b BLAST_FILE        BLAST tabulated output that was generated with
                       pipeline_parseblast.py
  -f FASTA_FILE        FASTA file to be used to output list of bad hits that
                       did not match thresholds
  -i IDENTITY[0-100]   hit identity in percentage to be accepted as a hit,
                       recommended 97
  -l ALIGNMENT[0-100]  hit aliginment length in percentage to be accepted a
                       hit, recommended 95
  -vs VARIABLE_START   reference sequence variable region start
  -ve VARIABLE_END     reference sequence variable region end
  -t BLAST_TYPE[0-2]   defines which section of the BLAST to be used to
                       summarize results. 0 - suitable for MaarjAM, only last
                       portion of hit description is used, 1 - all hit
                       description is used, 2 - hit identificator is used
```

## 5. Final results

Two files are generated from previous step. One is named as `*.nohits.fasta`, where sequences that did not get significant hit against reference database are written out and `*.tsv`, where results are written as pivot table with samples and hits sorted in descending order.

## 6. BLAST against INSDC (Optional)

To identify nohits we can use INSDC database to understand what else the sequences are containing. We first need to download database from NCBI FTP server and conduct BLAST on the downloaded database. BLAST can be run with all the INSDC data partitions together (large memory usage) or if memory usage is limited, by separately.

## 6.1. Download INSDC database from NCBI FTP server. 

Please check number of NT sequences in FTP to download them all, change `41` from below line accordingly. As the database contains only GenBank accessions, GenBank ID and short description, we need to download taxonomy information (`gi_taxid_nucl.dmp.gz`, `taxdump.tar.gz`) to build taxonomy tree for each hit. 
```
for i in {00..41}; do echo "Downloading NT.$i"; wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.$i.tar.gz; tar xzvf nt.$i.tar.gz; done
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
gunzip gi_taxid_nucl.dmp.gz
gunzip taxdump.tar.gz

```

## 6.2 Run BLAST against INSDC

To run BLAST with all the INSDC data partitions together, we can simply define database parameter in BLAST as `nt` or if we want to run them separately, we need to define all the partions one by one or use for loop.

```
blastn -query 454.cf.blast.i97.a95.nohits.fasta -evalue 1e-50 -max_target_seqs 1 -num_threads 4 -db nt -outfmt 5 | python pipeline_parse_blast.py > 454.nohits.blast
blastn -query illumina.cf.blast.i97.a95.nohits.fasta -evalue 1e-50 -max_target_seqs 1 -num_threads 4 -db nt -outfmt 5 | python pipeline_parse_blast.py > illumina.nohits.blast
```

or run separately and combine BLAST results together by selecting best hut for each sequence based on BLAST score (change number 41 accordinly to downloaded partitions)

```
for i in {00..41}; do blastn -query 454.cf.blast.i97.a95.nohits.fasta -evalue 1e-50 -max_target_seqs 1 -num_threads 4 -db nt.$i -outfmt 5 | python pipeline_parse_blast.py > 454.nohits.$i.blast; done
less 454.nohits.*.blast | python pipeline_merge_blasts.py > 454.nohits.blast
for i in {00..41}; do blastn -query illumina.cf.blast.i97.a95.nohits.fasta -evalue 1e-50 -max_target_seqs 1 -num_threads 4 -db nt.$i -outfmt 5 | python pipeline_parse_blast.py > illumina.nohits.$i.blast; done
less illumina.nohits.*.blast | python pipeline_merge_blasts.py > illumina.nohits.blast
```

## 6.3 Summarize BLAST results

We use relaxed parameters to filter potential hits by reducing identity threshold to be at least 90% and length at least 90%. We also need to define files, where taxonomy information is stored for `nt` database. Warning: as the `*.dmp` files are relatively large and below script is not optimized, it can use large ammount of memory and will take time to process file `gi_taxid_nucl.dmp`.

```
python pipeline_summarize_gbblast.py -b 454.nohits.blast -i 90 -l 90 -ti gi_taxid_nucl.dmp -tt names.dmp -tn nodes.dmp
python pipeline_summarize_gbblast.py -b illumina.nohits.blast -i 90 -l 90 -ti gi_taxid_nucl.dmp -tt names.dmp -tn nodes.dmp
```

Command help

```
python pipeline_summarize_gbblast.py [-h] -b BLAST_FILE -ti ID_FILE -tt
                                     TAXONOMY_FILE -tn NODE_FILE -i
                                     IDENTITY[0-100] -l ALIGNMENT[0-100]
arguments:
  -h, --help           show this help message and exit
  -b BLAST_FILE        BLAST tabulated output that was generated with
                       pipeline_parseblast.py
  -ti ID_FILE          Taxonomy file, where for each GenBank ID node ID is
                       specified
  -tt TAXONOMY_FILE    Taxonomy file, where for each node ID scientific name
                       is provided
  -tn NODE_FILE        Taxonomy file, where full tree node connections of node
                       IDs are provided to build full taxonomy tree
  -i IDENTITY[0-100]   hit identity in percentage to be accepted as a hit,
                       recommended 90
  -l ALIGNMENT[0-100]  hit aliginment length in percentage to be accepted a
                       hit, recommended 90
```

___

Use following citation when using our python scripts:
```
Vasar M, Andreson R, Davison J, Jairus T, Moora M, Remm M, Young JPW, Zobel M, Öpik M (XXXX) Increased sequencing depth does not increase captured diversity of arbuscular mycorrhizal fungi.
```

License: [CC-BY](https://creativecommons.org/licenses/by/3.0/)
