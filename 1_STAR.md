# Generating Solo files
by Sarah Salisbury

## Getting Started

### 1. Make a folder for your project

Generate the following folder hierarchy:
```
├──data <- this is where we'll put our raw samples
├──refmt <- this is where we'll put the relevant reference genome files and do our indexing
│   ├──star_genome_index <- this is where STAR will put the indexed genome
├──outmt <- this is where we'll do our mapping of our samples against the reference genome
│   ├──sample1 <- this is where we'll put our mapping output files for sample 1
│   ├──sample2 <- this is where we'll put our mapping output files for sample 2
│   ├──sample3 <- this is where we'll put our mapping output files for sample 3
```

### 2. Get the Data

First, unzip the samples in data directory:
```
cd data
tar zxvf first-sample.tgz
tar zxvf second-sample.tgz
tar zxvf third-sample.tar.gz
```
> [!NOTE]
> There are many R1, R2, I1, I2 files for each sample because each sample was run across multiple lanes during sequencing!

## Downloading Reference Genomes

### 1. Download from NCBI
Download ```.gff``` and ```.fna``` files for this shark species from: https://www.ncbi.nlm.nih.gov/genome/10771.

Place them in the refmt folder.

> [!NOTE]
> A mitochondrial genome is included in the reference genome.

### 2. Unzip the files

Unzip the fasta and gff files.
```
cd refmt
gunzip GCF_902713615.1_sScyCan1.1_genomic.fna.gz
gunzip GCF_902713615.1_sScyCan1.1_genomic.gff.gz
```

## Fixing the Mitochondrial Genome

Ok so the mitochondrial genome, despite being included in the annotation file has some issues.

First, there is no "exon" reported in the ```.gff``` file for protein coding mitochondrial genes (e.g., COI, ND1, etc.)

Second, for those non-coding genes (e.g., tRNAs, rRNAs) there is no gene ID reported.

Let's fix this.

### 1. Getting exons for protein coding mtDNA genes

Only CDS and gene is listed for each protein coding mtDNA gene, there is no RNA or exon. The latter in particular is a problem since STAR maps reads to exons!

You can see that by doing the following (note that the code for the mitochondrial genome in this species is "NC_001950.1":
```
grep 'NC_001950.1' GCF_902713615.1_sScyCan1.1_genomic.gff | head
```
Output:
```
##sequence-region NC_001950.1 1 16697
NC_001950.1     RefSeq  region  1       16697   .       +       .       ID=NC_001950.1:1..16697;Dbxref=taxon:7830;Is_circular=true;Name=MT;dev-stage=adult;gbkey=Src;genome=mitochondrion;isolate=single individual;mol_type=genomic DNA;tissue-type=muscle
NC_001950.1     RefSeq  gene    1       975     .       +       .       ID=gene-ND1;Dbxref=GeneID:808294;Name=ND1;gbkey=Gene;gene=ND1;gene_biotype=protein_coding;gene_synonym=NADH 1
NC_001950.1     RefSeq  CDS     1       975     .       +       0       ID=cds-NP_007614.1;Parent=gene-ND1;Dbxref=UniProtKB/Swiss-Prot:O21408,Genbank:NP_007614.1,GeneID:808294;Name=NP_007614.1;gbkey=CDS;gene=ND1;product=NADH dehydrogenase subunit 1;protein_id=NP_007614.1;transl_table=2
NC_001950.1     RefSeq  tRNA    979     1048    .       +       .       ID=rna-NC_001950.1:979..1048;gbkey=tRNA;product=tRNA-Ile
NC_001950.1     RefSeq  exon    979     1048    .       +       .       ID=exon-NC_001950.1:979..1048-1;Parent=rna-NC_001950.1:979..1048;gbkey=tRNA;product=tRNA-Ile
NC_001950.1     RefSeq  tRNA    1050    1122    .       -       .       ID=rna-NC_001950.1:1050..1122;gbkey=tRNA;product=tRNA-Gln
NC_001950.1     RefSeq  exon    1050    1122    .       -       .       ID=exon-NC_001950.1:1050..1122-1;Parent=rna-NC_001950.1:1050..1122;gbkey=tRNA;product=tRNA-Gln
NC_001950.1     RefSeq  tRNA    1123    1192    .       +       .       ID=rna-NC_001950.1:1123..1192;gbkey=tRNA;product=tRNA-Met
NC_001950.1     RefSeq  exon    1123    1192    .       +       .       ID=exon-NC_001950.1:1123..1192-1;Parent=rna-NC_001950.1:1123..1192;gbkey=tRNA;product=tRNA-Met
```

Ok so what we're going to try is changing the CDS rows to exon rows. We're going to try to be as minimally invasive as possible by first changing all instances of "CDS" to "exon" and "cds" to "exon".

You can check to make sure that "CDS" and "cds" aren't in any gene names or something that we wouldn't want to change with the following lines of code:
```
grep 'NC_001950.1' GCF_902713615.1_sScyCan1.1_genomic.gff | grep 'cds'
grep 'NC_001950.1' GCF_902713615.1_sScyCan1.1_genomic.gff | grep 'CDS'
```

Ok now that we're confident that we can change "CDS" and "cds" without messing up our file, let's substitute "CDS" and "cds" for "exon" using sed and a global (g) substitution, after first selecting for only those lines with the string 'NC_001950.1'.
```
sed '/NC_001950.1/s/CDS/exon/g' GCF_902713615.1_sScyCan1.1_genomic.gff > GCF_902713615.1_sScyCan1.1_genomic_fixedCDS.gff
sed '/NC_001950.1/s/cds/exon/g' GCF_902713615.1_sScyCan1.1_genomic_fixedCDS.gff > GCF_902713615.1_sScyCan1.1_genomic_fixedCDS_fixedcds.gff
```

Ok let's look at just those lines in the mtDNA genome to make sure this substitution was ok:
```
grep 'NC_001950.1' GCF_902713615.1_sScyCan1.1_genomic_fixedCDS_fixedcds.gff
```

### 2. Getting gene ID for non-coding RNAs in mtDNA genome

Ok so we keep getting an error during our indexing in the ```Log.out``` file like the following:
```
WARNING: while processing pGe.sjdbGTFfile=GCF_902713615.1_sScyCan1.1_genomic.gtf: no gene_id for line:
NC_001950.1     RefSeq  exon    8889    8957 
```

I believe this is the reason why in our features.tsv file in our Solo.out directory we get counts for "MissingGeneID".

I think this comes down to the non-coding rRNAs and tRNAs being labeled as "ID=rna-something". So what I'm going to try is to replace the "rna" with "gene" here, both for the exon line (which refers to the Parent id as "Parent=rna-something") and for the non-coding line itself (labeled tRNA or rRNA).
Have a look at this structure here:
```
grep 'NC_001950.1' GCF_902713615.1_sScyCan1.1_genomic_fixedCDS_fixedcds.gff
```
E.g. you get something like this for each tRNA/rRNA. You have a line for the "exon" and a line for the "tRNA".
```
NC_001950.1     RefSeq  tRNA    1050    1122    .       -       .       ID=rna-NC_001950.1:1050..1122;gbkey=tRNA;product=tRNA-Gln
NC_001950.1     RefSeq  exon    1050    1122    .       -       .       ID=exon-NC_001950.1:1050..1122-1;Parent=rna-NC_001950.1:1050..1122;gbkey=tRNA;product=tRNA-Gln
```
Note that for the "tRNA" line: "ID=rna", and the third column is "tRNA" not "gene". For the "exon" line, the Parent="rna code of tRNA/rRNA".

So, for those lines in our gff with "NC_001950.1", substitute (s) "rna-" with "gene-" (this will fix the ID for the tRNA/rRNA and fix the Parent= for the exon!).
```
sed '/NC_001950.1/s/rna-/gene-/g' GCF_902713615.1_sScyCan1.1_genomic_fixedCDS_fixedcds.gff > GCF_902713615.1_sScyCan1.1_genomic_fixedCDS_fixedcds_fixedncRNA.gff
```

I'm also going to replace "tRNA" and "rRNA" in the third column with "gene". To do this I take only those lines with "NC_001950.1", then I substitute (s), a ```tab``` (\t) before "tRNA" (because tRNA appears elsewhere in some lines and I don't want to change those) with a ```tab``` followed by "gene", globally (g). Then I do the same for "rRNA".
```
sed '/NC_001950.1/s/\ttRNA/\tgene/g' GCF_902713615.1_sScyCan1.1_genomic_fixedCDS_fixedcds_fixedncRNA.gff > GCF_902713615.1_sScyCan1.1_genomic_fixedCDS_fixedcds_fixedncRNA_fixedtRNAlabel.gff
sed '/NC_001950.1/s/\trRNA/\tgene/g' GCF_902713615.1_sScyCan1.1_genomic_fixedCDS_fixedcds_fixedncRNA_fixedtRNAlabel.gff > GCF_902713615.1_sScyCan1.1_genomic_fixedCDS_fixedcds_fixedncRNA_fixedtRNAlabel_fixedrRNAlabel.gff
```

### 3. Checking the altered .gff file
Have a look at the mtDNA genes to make sure you're happy:
```
grep 'NC_001950.1' GCF_902713615.1_sScyCan1.1_genomic_fixedCDS_fixedcds_fixedncRNA_fixedtRNAlabel_fixedrRNAlabel.gff
```


E.g. you get something like this for each protein-coding gene. You have a line for the "exon" and a line for the "gene". Note that CDS and cds have been changed to "exon".
```
NC_001950.1     RefSeq  gene    1       975     .       +       .       ID=gene-ND1;Dbxref=GeneID:808294;Name=ND1;gbkey=Gene;gene=ND1;gene_biotype=protein_coding;gene_synonym=NADH 1
NC_001950.1     RefSeq  exon    1       975     .       +       0       ID=exon-NP_007614.1;Parent=gene-ND1;Dbxref=UniProtKB/Swiss-Prot:O21408,Genbank:NP_007614.1,GeneID:808294;Name=NP_007614.1;gbkey=exon;gene=ND1;product=NADH dehydrogenase subunit 1;protein_id=NP_007614.1;transl_table=2
```

E.g. you get something like this for each tRNA/rRNA. You have a line for the "exon" and a line for the "tRNA". Note that tRNA/rRNA in the third column has been changed to "gene", and that "rna-" in the exon "Parent=" spot and in the gene "ID=" spot have been changed to "gene-".
```
NC_001950.1     RefSeq  gene    1050    1122    .       -       .       ID=gene-NC_001950.1:1050..1122;gbkey=tRNA;product=tRNA-Gln
NC_001950.1     RefSeq  exon    1050    1122    .       -       .       ID=exon-NC_001950.1:1050..1122-1;Parent=gene-NC_001950.1:1050..1122;gbkey=tRNA;product=tRNA-Gln
```

Have a look at the first 50 lines of the altered gff to make sure CDS and -rna weren't changed in other chromosomes:
```
head -n 50 GCF_902713615.1_sScyCan1.1_genomic_fixedCDS_fixedcds_fixedncRNA_fixedtRNAlabel_fixedrRNAlabel.gff
```
Looks good.

Compare length of original .gff with new .gff
```
wc -l GCF_902713615.1_sScyCan1.1_genomic.gff
1392858 GCF_902713615.1_sScyCan1.1_genomic.gff

wc -l GCF_902713615.1_sScyCan1.1_genomic_fixedCDS_fixedcds_fixedncRNA_fixedtRNAlabel_fixedrRNAlabel.gff
1392858 GCF_902713615.1_sScyCan1.1_genomic_fixedCDS_fixedcds_fixedncRNA_fixedtRNAlabel_fixedrRNAlabel.gff
```

### 4. Getting mtDNA gene IDs for subsequent Seurat analysis

Ok so you'll want to know the IDs for each of the mtDNA genes for downstream Seurat analysis. We can get those using the following code. This greps the mtDNA code "NC_001950.1" in our gff file, then greps only those lines that have "ID=gene" (excluding exons), then we select for everything after "ID=", and then remove everything after the first ";" in each line, then sort and get only the unique values.
```
grep "NC_001950.1" GCF_902713615.1_sScyCan1.1_genomic_fixedCDS_fixedcds_fixedncRNA_fixedtRNAlabel_fixedrRNAlabel.gff | grep 'ID=gene' | sed 's/.*ID=\(.*\)/\1/' | sed 's/;.*//'| sort | uniq
```
Output:
```
gene-ATP6
gene-ATP8
gene-COX1
gene-COX2
gene-COX3
gene-CYTB
gene-NC_001950.1:1050..1122
gene-NC_001950.1:1123..1192
gene-NC_001950.1:11445..11514
gene-NC_001950.1:12661..12732
gene-NC_001950.1:12733..12801
gene-NC_001950.1:13852..13920
gene-NC_001950.1:13921..14877
gene-NC_001950.1:14878..14949
gene-NC_001950.1:14950..16622
gene-NC_001950.1:16623..16697
gene-NC_001950.1:2239..2307
gene-NC_001950.1:2309..2377
gene-NC_001950.1:2378..2450
gene-NC_001950.1:2488..2554
gene-NC_001950.1:2556..2625
gene-NC_001950.1:4181..4251
gene-NC_001950.1:4256..4325
gene-NC_001950.1:5025..5098
gene-NC_001950.1:6729..6798
gene-NC_001950.1:7148..7217
gene-NC_001950.1:8889..8957
gene-NC_001950.1:8958..9024
gene-NC_001950.1:9025..9096
gene-NC_001950.1:979..1048
gene-ND1
gene-ND2
gene-ND3
gene-ND4
gene-ND4L
gene-ND5
gene-ND6
```

Great! Now we can move forward!

## Read alignment and expression matrices using STAR

### 1. Convert GFF to GTF file

You will need to transform the gff genome annotation file to gtf format for STAR, you can do this using gffread.

Stay in the refmt directory.

Now run gffread:
```
gffread GCF_902713615.1_sScyCan1.1_genomic_fixedCDS_fixedcds_fixedncRNA_fixedtRNAlabel_fixedrRNAlabel.gff -T -o GCF_902713615.1_sScyCan1.1_genomic_fixedCDS_fixedcds_fixedncRNA_fixedtRNAlabel_fixedrRNAlabel.gtf
```
NOTES:
- The -T option specifies that you want it in GTF2 version, rather than default GFF3 version
- The -o option specifies you want to name the output file something specific (as specified after the flag).

Then the genome and its annotation can be indexed with STAR.

### 2. Get Cell Barcode (CB) Whitelist

So all of our nuclei/cells are barcoded with a unique cell barcode (CB). These are the barcodes that are unique to each bead, and because each nuclei is paired with a different bead during the 10X protocol, each nuclei has a unique barcode. We can get a list of each of these potential barcodes here: https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/3M-february-2018.txt.gz. During the mapping stage we will remove sequences without a cell barcode that is on this list.

So place the file in the refmt directory and then unzip it:
```
gunzip 3M-february-2018.txt.gz
```

### 3. Make output folder for indexing step
So we need to make a folder within the refmt directory that we will then output our STAR outputs, this is specified in the indexing code below by ```--genomeDir star_genome_index```. Please note that this directory must be empty before you run STAR.

If you haven't already, make star_genome_index within the refmt directory:
```
mkdir star_genome_index
```

### 4. Indexing

So before we look at our samples, we need to prepare our reference genome using STAR. This indexing step will make a STAR annotated genome that we can then map our samples to in order to generate our matrices of nuclei/cells and features.

```
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir star_genome_index --genomeFastaFiles GCF_902713615.1_sScyCan1.1_genomic.fna --sjdbGTFfile GCF_902713615.1_sScyCan1.1_genomic_fixedCDS_fixedcds_fixedncRNA_fixedtRNAlabel_fixedrRNAlabel.gtf
```

NOTES:

```--limitGenomeGenerateRAM``` indicates total amount of RAM STAR can use

```--runThreadN``` option tells you how many threads to have

```--runMode genomeGenerate``` option directs STAR to run genome indices generation job

```--genomeDir``` specifies path to the directory (henceforth called ”genome directory” where the genome indices are stored. This empty directory has to be created (with mkdir) before STAR run and needs to have writing permissions. The file system needs to have at least 100GB of disk space available for a typical mammalian genome. It is recommended to remove all files from the genome directory before running the genome generation step. This directory path will have to be supplied at the mapping step to identify the reference genome.

```--genomeFastaFiles``` specifies your fasta file

```--sjdbGTFfile``` /path/to/annotations.gtf

### 5. Mapping
Now we're ready to map our samples against the indexed genome we made in the previous step.

You will want to run this step once per sample in the directory for each given sample.

E.g. for sample 1 mapping you'll want to run this in the following directory:
```
cd outmt/sample1
```

Run mapping independently for each sample.
E.g. for sample 1:
```
STAR --genomeDir ../../refmt/star_genome_index --readFilesIn ../../data/raw-first-sample/SN340_preliminary_R2_001.fastq.gz,../../data/raw-first-sample/SN340_S1_L001_R2_001.fastq.gz,../../data/raw-first-sample/SN340_S1_L002_R2_001.fastq.gz,../../data/raw-first-sample/SN340_S1_L003_R2_001.fastq.gz,../../data/raw-first-sample/SN340_S1_L004_R2_001.fastq.gz ../../data/raw-first-sample/SN340_preliminary_R1_001.fastq.gz,../../data/raw-first-sample/SN340_S1_L001_R1_001.fastq.gz,../../data/raw-first-sample/SN340_S1_L002_R1_001.fastq.gz,../../data/raw-first-sample/SN340_S1_L003_R1_001.fastq.gz,../../data/raw-first-sample/SN340_S1_L004_R1_001.fastq.gz --soloMultiMappers Unique EM --soloType CB_UMI_Simple --soloUMIlen 12 --soloCBwhitelist ../../ref/3M-february-2018.txt --soloFeatures GeneFull --clipAdapterType CellRanger4 --outFilterScoreMin 20 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR  --runThreadN 4 --outMultimapperOrder Random --outFilterMultimapScoreRange 1 --outFilterMultimapNmax 10 --outFilterMismatchNmax 10  --readFilesCommand zcat --outSAMtype BAM Unsorted
```

NOTES:

```--genomeDir``` is the location where your indexed genome is

```--readFilesIn``` is your raw library reads (read2 then read1)

```--soloMultiMappers Unique EM``` is the counting method for reads mapping to multiple genes. So "Unique" is the default, where you count only reads that map to unique genes. "EM" specifies that UMIs which map to multiple genes are "counted" distributed across all the genes they map to using the Expectation Maximization algorithm, so basically you're splitting up your UMI "count" over many possible genes (so your final count table might have non-integer values).

```--soloType CB_UMI_Simple``` indicates that you have 10X Chromium libraries (so star should look for one UMI and one cell barcode)

```--soloUMIlen``` the length of the UMI, the default here is 10 but the new chemistry of 10X uses 12bp UMI, so this must be specified!

```--soloCBwhitelist``` this is our list of cell barcodes (CB) we expect (this is the file we got from Chromium earlier)

```--soloFeatures GeneFull``` so this is specifying which genomic features we want associated with our UMI and cell barcodes. "GeneFull" indicates that we want full gene (pre-mRNA) info, including if reads are in genes' exons and introns.

```--clipAdapterType CellRanger4``` so this just specifies that we want star to cut out our adapters (5p and 3p) as is also done in CellRanger4, Utilizes Opal package by Martin ˇSoˇsi ́c: https://github.com/Martinsos/opal.

```--outFilterScoreMin``` integer: alignment will be output only if its score is higher than or equal to this value, default = 0

```--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts``` multiple matches in whitelist with 1 mismatched base allowed, posterior probability calculation is used choose one of the matches AND pseudocounts of 1 are added to all whitelist barcodes AND  multimatching to whitelist is allowed for CBs with N-bases. This option matches best with CellRanger >=3.0.0

```--soloUMIfiltering MultiGeneUMI_CR``` remove UMIs with N and homopolymers AND remove lower-count UMIs that map to more than one gene, matching CellRanger > 3.0.0

```--soloUMIdedup 1MM_CR``` CellRanger2-4 algorithm for 1MM UMI collapsing

```--runThreadN``` defines the number of threads to be used for genome generation, it has to be set to the number of available cores on the server node

```--outMultimapperOrder Random``` order of multimapping alignments in the output files random order of alignments for each multi-mapper. Read mates (pairs) are always adjacent, all alignment for each read stay together. This option will become default in the future releases.

```--outFilterMultimapScoreRange 1``` integer: the score range below the maximum score for multimapping alignments, default is 1

```--outFilterMultimapNmax 10``` int: maximum number of loci the read is allowed to map to. Alignments (all of them) will be output only if the read maps to no more loci than this value. Otherwise no alignments will be output, and the read will be counted as ”mapped to too many loci” in the Log.final.out. Default is 10

```--outFilterMismatchNmax 10``` integer: alignment will be output only if it has no more mismatches than this value. Default is 10

```--readFilesCommand zcat``` command line to execute for each of the input file. This command should generate FASTA or FASTQ text and send it to stdout, zcat will uncompress .gz files

```--outSAMtype BAM Unsorted``` output BAM without sorting

SANITY CHECKS: IS YOUR ```--genomeDir``` THE SAME AS THE ONE YOU JUST MADE IN THE INDEX STEP?

IS ```--readFilesIn``` ASSOCIATED WITH THE CORRECT FILES AND IS READ2 BEFORE READ1?
