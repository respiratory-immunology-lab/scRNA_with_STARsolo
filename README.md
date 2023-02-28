# scRNA_with_STARsolo
Processing DB Rhapsody scRNA data with STARsolo

# Setup

## Download STAR

Download and compile the latest version of STAR (https://github.com/alexdobin/STAR).

You can find the STAR manual here https://raw.githubusercontent.com/alexdobin/STAR/master/doc/STARmanual.pdf.

## Download Genome

Download species genome from GENCODE (recommended - https://www.gencodegenes.org/) or ENSEMBL (latest release - http://ftp.ensembl.org/pub/release-104/)

Download both GTF file and fasta (dna) file for your species. The most comprehensive annotations for yours species is strongly recommended.


Example using Mouse Release M27 (GRCm39)
```
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M27/gencode.vM27.chr_patch_hapl_scaff.annotation.gtf.gz
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M27/GRCm39.genome.fa.gz
```

### If samples have been tagged




## Download BD Rhapsody Whitespace files

The structure of BD Rhapsody is explained here https://teichlab.github.io/scg_lib_structs/methods_html/BD_Rhapsody.html

STARsolo requires all 3 of BD Rhapsodys' whitespace files.

```
# I chose to download my whitespace files to the db directory
cd {path}/db/
wget https://teichlab.github.io/scg_lib_structs/data/BD_CLS1.txt
wget https://teichlab.github.io/scg_lib_structs/data/BD_CLS2.txt
wget https://teichlab.github.io/scg_lib_structs/data/BD_CLS3.txt
```

# Generating genome indexes files

From STAR manual.

The basic options to generate genome indices are as follows:

```
--runThreadN NumberOfThreads
--runMode genomeGenerate
--genomeDir /path/to/genomeDir
--genomeFastaFiles /path/to/genome/fasta1 /path/to/genome/fasta2 ...
--sjdbGTFfile /path/to/annotations.gtf
--sjdbOverhang ReadLength-1
```

Create genomeDir before running.

Example bash script 
```
#!/bin/bash

dir="{path}/db/"

/path/to/STAR-2.7.10a/source/STAR \
	--runThreadN 20 \
	--genomeDir "${dir}STARgenomeIndex" \
	--runMode genomeGenerate \
	--genomeFastaFiles "${dir}GRCm39.genome.fa" \
	--sjdbGTFfile "${dir}gencode.vM27.chr_patch_hapl_scaff.annotation.gtf" \
	--sjdbOverhang 99
  
```

# Mapping reads to the genome 

## BD Rhapsody

Example bash script.

```
#!/bin/bash

genomeDir="path/to/DB"
seqDir="/path/to/sequencing/files"

/path/to/STAR-2.7.10a/source/STAR \
        --outFileNamePrefix "${seqDir}/sample1" \
        --genomeDir "${genomeDir}/STARgenomeIndex" \
        --runThreadN 20 \
        --readFilesIn "${seqDir}/raw_fastq/sample1_R2.fastq.gz" "${seqDir}/raw_fastq/sample1_R1.fastq.gz" \
        --readFilesCommand zcat \
        --sjdbGTFfile "${genomeDir}/db/gencode.vM27.chr_patch_hapl_scaff.annotation.gtf" \
        --soloType CB_UMI_Complex \
        --soloCBmatchWLtype 1MM \
        --soloCBwhitelist "${genomeDir}/BD_CLS1.txt" "${genomeDir}/BD_CLS2.txt" "${genomeDir}/BD_CLS3.txt" \
        --soloUMIlen 8 \
        --soloCBposition 0_0_0_8 0_21_0_29 0_43_0_51 \
        --soloUMIposition 0_52_0_59 \
        --soloCellFilter EmptyDrops_CR 3000 0.99 10 45000 90000 500 0.01 20000 0.01 10000
```

## BD Rhapsody With Sample tags

```
#!/bin/bash

genomeDir="path/to/DB"
seqDir="/path/to/sequencing/files"

/path/to/STAR-2.7.10a/source/STAR \
        --outFileNamePrefix "${seqDir}/sample1_w_tags" \
        --genomeDir "${genomeDir}/STARgenomeIndex_sampleTag" \
        --runThreadN 19 \
        --readFilesIn "${seqDir}/raw_fastq/sample1_w_tags_R2.fastq.gz" "${seqDir}/raw_fastq/s/sample1_w_tags_R1.fastq.gz" \
        --readFilesCommand zcat \ ## if files not gzip remove this line
        --sjdbGTFfile "${genomeDir}/db_sampleTag_BD_Mus/gencode.vM29.annotation.gtf" \
        --soloType CB_UMI_Complex \
        --soloCBmatchWLtype 1MM \
        --soloCBwhitelist "${genomeDir}/BD_CLS1.txt" "${genomeDir}/BD_CLS2.txt" "${genomeDir}/BD_CLS3.txt" \
        --soloUMIlen 8 \
        --soloCBposition 0_0_0_8 0_21_0_29 0_43_0_51 \
        --soloUMIposition 0_52_0_59 \
        --alignIntronMax 1 \
        --alignMatesGapMax 1 \
        --winAnchorMultimapNmax 1 \
        --chimMainSegmentMultNmax 1 \
        --outFilterMultimapNmax 1 \
        --seedMultimapNmax 1 \
        --outFilterScoreMinOverLread 0 \
        --outFilterMatchNminOverLread 0

```



