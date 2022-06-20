# SuffixAligner
A Long Read Aligner Based on Suffix Array Construction Algorithm for DNA Alphabets 
SuffixAligner
SuffixAligner is a python-based aligner for long noisy reads generated from third-generation sequencing machines. SuffixAligner exploits the nature of biological alphabet that has a fixed-size and a predefined lexical ordering to construct a suffix array for indexing a reference genome. FM-index is used to efficiently search the indexed reference and locate the exact matched seeds among reads and the reference. The matched seeds are arranged into windows/clusters and the ones with the maximum number of seeds are reported as candidates for mapping positions, more details about SuffixAligner can be found in: 

•	(Indexing) https://doi.org/10.1016/j.jksuci.2022.04.015

•	(Mapping) 


# System requirements

64-bit machine with python in general, numpy, and pandas libraries.

# Installation
1.	Downlowd files
2.	Run code f5

# Quick usage guide

## Indexing: 

#### Run
Find_SA_for_Genome(L=100,l=100,G_file='genome_file.fasta')    
#### Or Run
Find_SA_for_Genome()


            *[G_file] Reference Genome File     	                  [default:example1.fasta]

            *[L] Divide genome to substrings of length L     	      [default:100]

            *[l] Overlap length                              	      [default:100]


#### output
[G_file].SA.txt     the suffix array in string format 


## Mapping:

#### Run

Mapping(Type="r",G_file="example1.fasta",
            R_file="Read1.fastq",
            SA_file="example1.fasta.SA.txt",
            Start=0,End=10)     
#### Or run

Mapping(Type="s",G_file=" example1.fasta",
            SA_file="example1.fasta.SA.txt",
            Sam_file="example1.bwa.sam",
            Start=0,End=10)

#### Or Run

Mapping()


            *[Type] "r" for read, "s" for sam file             [default:r]       

            *[G_file] Reference Genome File                    [default:example1.fasta]

            *[R_file] Read File                                [default:example1.fasta]

            *[SA_file] Suffix array File                       [default:example1.fasta]

            *[Start] read number which Start from              [default:0]

            *[End] number of read which end with               [default:10] 

            *[Sam_file] Sam file which generated from other aligner  



# Notes
•	



# Read files
SuffixAligner align sequencing reads given in fastq format to reference genome given in fasta format. Also, SuffixAligner can read directly the sam files given in sam format which genertated from other aligner and search for solution for unmapped read.
# Outputs
The output of SuffixAligner is suffix array in text format from the step of indexing in the file:

            [G_file].SA.txt

The output of the step of mapping is sam file in the file:

            [R_file].[Start].[End].sam      when input is read file

#### Or

            [Sam_file].[Start].[End].sam      when input is sam file


SuffixAligner also reports the following on the screen:

•	Genome length

•	Number of read

•	Number of read in sam file

•	Number of unmapped read in sam file


