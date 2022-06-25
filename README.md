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

## Genome Indexing Stage: 

``` python Indexing.py -L [Substrings_Length] -G [Genome_file] ```
``` 
* [-L or --Length ] Divide genome to substrings of length L             [default: 100]
* [-G or --Genome ] Reference Genome File                               [default: example1.fasta]
* [-h]              Help
```

#### Notes:
- The output of Indexing stage is a suffix array file named with the ` [Genome_prefix_file_name].SA.txt `.
- If you did not specify any parameters, SuffixAligner will invoke its default parameters setting with the default reference genome file `example1.fasta` 


## Mapping:
### Case 1: Dealing with Sequencing Reads Files 

``` python Mapping.py  -T [r] -G [Genome_file]  -R [Read_File_1]  -R [Read_File_2]  -R [Read_File_3]  -B [Begining_Read_Index]  -E [Ending_Read_Index] -O [Output_File_Prefix]```

### Case 2: Resolve unmapped reads by dealing with SAM Files produced by other Aligners

``` python Mapping.py  -T [s] -G [Genome_file]  -F [SAM_File_1]  -F [SAM_File_2]  -F [SAM_File_3]  -B [Begining_Read_Index]  -E [Ending_Read_Index] -O [Output_File_Prefix]```

``` 
* [-T] Type of the input files (i.e. r for Reads, s for SAM)                  [default: r]
* [-G or --Genome ] Reference Genome File                                     [default: example1.fasta]
* [-R, --Reads] Sequencing Read files, each file is specified with -R option. [default: Read_example1.fastq]
* [-B, --Begin] Starting Read Index                                           [default:0]
* [-E, --End] Ending Read Index                                               [default:All Reads is used] 
* [-F, --Sam_File] SAM files, each file is specified with -F option. 
* [-O, --Output] Prefix name of the output SAM file produced by SuffixAligner [default: <Genome_File_Prefix>.sam]
* [-h] Help
```
#### Notes:
- Starting and Ending Reads Indicies are used when you are interested in a small subset of reads and the indicies are corresponding to the reads positions in the sequencing files.
- The default parameters setting of the SuffixAligner is considering the whole set of sequencing data. 
- The SAM option is used when you would like to improve the number of mapping reads produced by the other aligners.   
- SuffixAligner also reports the following on the screen: <br>
        * Genome length. <br>
        * Number of Sequencing Reads in both cases (i.e. reading sequencing data or SAM files). <br> 
        * Number of unmapped Reads in the SAM files produced by the other aligners.


