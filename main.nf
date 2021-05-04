
/*
Pipeline steps:
    1. Trimming
    FastQC (post-trim) -- perform post-trim FastQC on fastq files (ensure trimming performs as expected)
    3. Mapping w/ HISAT2 -- map to genome reference file
    4. SAMtools -- convert SAM file to BAM, index BAM, flagstat BAM
    
=======
*/

/*
 * PREPROCESSING - Build HISAT2 index 
 */

if (params.genome_path)
{
channel.fromPath(params.genome_path)
           .ifEmpty { exit 1, "genome file not found: ${params.genome_path}" }
           .set { genome_fasta}

}
if (params.fasta_path)

{
channel.fromFilePairs(params.fasta_path)
           .ifEmpty { exit 1, "fasta file not found: ${params.fasta_path}" }
           .set {reads}

}


process index
{
tag "$gm"
input: 
file gm from genome_fasta

output:
file '*.ht2' into index_ch


script:
"""
hisat2-build ${gm} ${gm}

"""
} 

process fastqc_trim

{
input:
tuple val(sampleID), path(files) from reads

output:
tuple val(sampleID), file("*val_?.fq") into trimmed_reads_hisat2

script:

"""
trim_galore --fastqc --paired ${files[0]} ${files[1]}

"""

}



process hisat2
{

input:
val(index) from index_ch
tuple val(name), file(trimmed_reads) from trimmed_reads_hisat2

script:


"""
echo ${name}
echo ${index}
"""
}

