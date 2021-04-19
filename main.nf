

#!/usr/bin/env nextflow
/*
Pipeline steps:
    1. Trimming
    FastQC (post-trim) -- perform post-trim FastQC on fastq files (ensure trimming performs as expected)
    3. Mapping w/ HISAT2 -- map to genome reference file
    4. SAMtools -- convert SAM file to BAM, index BAM, flagstat BAM
    
=======
*/


genome = channel
          .fromPath(params.genome_path)

process index
{
input: 
file gm from genome

output:
path 'SM_V9_21Feb' into index_ch

script:
"""
hisat2-build SM_V9_21Feb.fa SM_V9_21Feb

"""
}



reads = channel
          .fromFilePairs(params.fasta_path)

process fastqc

{
input:
    tuple val(sampleID), path(files) from reads

script:

"""
trim_galore --fastqc --paired ${files[0]} ${files[1]}

 
"""

}

