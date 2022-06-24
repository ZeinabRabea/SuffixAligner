# SuffixAligner
<p align="justify">
SuffixAligner is a python-based aligner for long noisy reads generated from third-generation sequencing machines. SuffixAligner exploits the nature of biological alphabet that has a fixed-size and a predefined lexical ordering to construct a suffix array for indexing a reference genome. FM-index is used to efficiently search the indexed reference and locate the exact matched seeds among reads and the reference. The matched seeds are arranged into windows/clusters and the ones with the maximum number of seeds are reported as candidates for mapping positions, more details about the algorithm used in the suffix array indexing stage for the SuffixAligner can be found in:          

Rabea, Z.; El-Metwally, S.;Elmougy, S. and Zakaria,M.; [A fast algorithm for constructing suffix arrays for DNA alphabets](https://www.sciencedirect.com/science/article/pii/S1319157822001434) .Journal of King Saud University - Computer and Information Sciences; 2022.  
<br>

Copyright (C) 2020-2022, and GNU GPL, by Zeinab Rabea, Sara El-Metwally,Samir Elmougy and Magdi Zakaria. </p>
           
### System Requirements
64-bit machine with python in general, numpy, and pandas libraries.

### Quick usage guide
1. Clone the [GitHub repo](https://github.com/ZeinabRabea/SuffixAligner), e.g. with `git clone https://github.com/ZeinabRabea/SuffixAligner.git`
2. For Indexing step,run `python Indexing.py -h` in the repo directory. 
3. For Mapping step,run `python Mapping.py -h` in the repo directory. 

### Genome Indexing Stage: 

1.	Run `python Indexing.py -L <Substrings_Length> -G <Genome_file> -G <Genome_file>`.
2.	Or run `python Indexing.py -G <Genome_file>`.
3.	Or run `python Indexing.py --Length <Substrings_Length> --Genome <Genome_file>`.
4.	Or run `python Indexing.py --Genome <Genome_file>`.
5.	Or run `python Indexing.py` in the repo directory.
6.	Or run `python Indexing.py -h`    


        *[-L, --Length] Divide genome to substrings of length L   [default:100]
        *[-G, --Genome] Reference Genome File                     [default:example1.fasta]
        *[-h] for help


#### The output
<Genome_file>.SA.txt         the suffix array in string format


## Mapping:

#### Run

1.	Run `python Mapping.py  -T r -G <Genome_file>  -R <Read_file>   -B <begining_read>  -E <Ending_read> -O <Output_file>`.
2.	Or run `python Mapping.py  -T r -G <Genome_file>  -R <Read_file>  -R <Read_file2>  -R <Read_file3>  -B <begining_read>  -E <Ending_read> -O <Output_file>`.
3.	Or run `python Mapping.py  -T s -G <Genome_file>  -F <sam_file>  -B <begining_read>  -E <Ending_read> -O <Output_file>`.
4.	Or run `python Mapping.py  -T s -G <Genome_file>  -F <sam_file> -F <sam_file2> -B <begining_read>   -E <Ending_read> -O <Output_file>`.
5.	Or run `python Mapping.py  -T r --Genome <Genome_file>  --Reads <Read_file>   --Begin <begining_read>  --End <Ending_read> --Output <Output_file>`.
6.	Or run `python Mapping.py  -T s --Genome <Genome_file> --Sam_File <sam_file>  --Begin <begining_read>  --End <Ending_read> --Output <Output_file>`.
7.	Or run `python Mapping.py  -T r -G <Genome_file>  -R <Read_file>`.
8.	Or run `python Mapping.py  -T s -G <Genome_file>  -F <sam_file>`.
9.	Or run `python Mapping.py`.
10.	Or run `python Mapping.py -h`.



        *[-T] take "r" align reads to genome 
        take "s" align Unmapped reads from sam file to genome     [default:r] 
        
        *[-G, --Genome] Reference Genome File                    [default:example1.fasta]
        *[-R, --Reads] Read File                                [default: Read_example1.fastq]
        *[-B, --Begin] Starting read number                     [default:0]
        *[-E, --End] Ending Read number                          [default:all number of read ] 
        *[-F, --Sam_File] Sam file which generated from other aligner  
        *[-O, --Output] Name of Sam file                         [default: <Genome_file>.sam]




# Notes

,,,,,,,,,,,,,,,,



# Read files
SuffixAligner align sequencing reads given in fastq format to reference genome given in fasta format. Also, SuffixAligner can read directly the sam files given in sam format which genertated from other aligner and search for solution for unmapped read.
# Outputs
The output of SuffixAligner is suffix array in text format from the step of indexing in the file:

            <Genome_file>.SA.txt

The output of the step of mapping is sam file in the file:

            <Genome_file>.sam     
            Or
            <Output_file>.sam 


SuffixAligner also reports the following on the screen:

•	Genome length

•	Number of read

•	Number of read in sam file

•	Number of unmapped read in sam file


