# DEEP
Introduction
============

DEEP: 

    DEEP is a RNA-seq aligner for second generation sequencing (SGS) RNA-seq. The popurse of DEEP is to find out more information about transcripts from RNA-seq data. DEEP mapping fasta or fastq RNA-seq to reference with fasta formate. Reference files will contain a whole genome or more genomes, it also can be a transcipt or only a peiece of sequence with fasta formate. DEEP also accept other aligners primary align results as input for deep mode align. DEEP primary align step need a BWT index of reference, we will use BWA build BWT index. After that DEEP find overlapped max exact match (MEM) seeds of reads with BWT index, and find acceptable seed combinations, and achieve align results. DEEP output align result in SAM format, enabling interoperation with a large number of other tools that use SAM. At present, DEEP can only run on Linux operating system.

Version
============

    2.0.0

Memory usage
============
    
    DEEP can align RNA-seq which sequencing from whole human transcripts to hg19 human reference within 16GB RAM memory. If input RNA-seq data is huge, DEEP need more SATA disk to record and sort reads. 

Installation
============

    cd DEEP/src
    make

Commands and options
============

    Download BWA:
        https://sourceforge.net/projects/bio-bwa/files/
    [build BWT index]:
        ./bwa index -p <bwa_index> ref.fa
        
    [build index]: 
        ./DEEP index <input_fasta_file_route> <output_index_route>
    [align]: 
        ./DEEP complete <-p/-f> -B <bwa_index> -H <output_index_route> -1 <input_file> (-2 <paired_input_file>) -O <output_route> <other options>

        if accept primary align results from other aligners.
        ./DEEP micro -B <bwa_index> -H <output_index_route> -1 <input_file(SAM format)> -O <output_route> <other options>


    [options]:

        <basic> :
       -B <STRING>:   reference BWT index path
       -H <STRING>:   hash files path, each need file's name must be *.hash with *.ann files
       -f         :   input files with fasta format
       -q         :   input files with fastq formate
       -1 <STRING>:   input files, need *.fa,*.fq files,split with ','
       -2 <STRING>:   input pair-end files, need *.fa,*.fq files,split with ','
       -O <STRING>:   output path,need enough space

       <alignment parameter> :
       -t <INT>:      output read score threshold(90)
       -r <INT>:      search area(500000)
       -a :           print all read sites while coverage score and alignment score over set threld will be used
       -b :           print best read sites while coverage score and alignment score over set threld will be used
       -p <INT>:      thread(1)
       -m <INT>:      match score(1)
       -s <INT>:      miss score(-1)
       -g <INT>:      gap score(-3)
       -e <INT>:      min exon length(12)
       -d             deep mode(0)

Evaluation
============

    simulated data sets are available at: http://bioinf.itmat.upenn.edu/BEERS/bp1/
    evaluation programs are available at: 

Contact
============

    For advising, bug reporting and requiring help, please contact ye.wenzhu@gmail.com
