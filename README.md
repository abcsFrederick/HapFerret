# HapFerret

The Expectation-Maximization (EM) algorithm for haplotype inference has been superseded for accuracy by newer programs, e.g. `PHASE` and `SHAPEIT`, but remains useful for rapid analysis. 

`HapFerret` is an EM implementation that is characterized by flexibility and ease of use, notable its use of a natural format for input genotypes and output haplotypes. Genotypes may be input as lines of comma separated alleles (csv), where the alleles may be any alphanumeric; i.e. for HLA data, for the loci MogCA—DQA1—rs9273349—DQB1—rs1894407—TAP2—TAP1, a hypothetical genotype input line: 9,12 101,301 C,T 302,503 A,A 3.2,4.2 3.1,3.1 might yield an output line of haplotype calls 12-301—C—302—A—4.2-3.1; 9-101—T—503—A—3.2—3.1. Input can also be by one locus genotype per line, with the same allele format. Critically for the inference of HLA haplotypes, an arbitrary number of alleles are allowed at each locus.

A special feature is the ability to infer haplotypes for all subsequences with a specified range of lengths. This has the particular purpose of finding blocks in which haplotypes are inferred with relatively low ambiguity. The program measures the ambiguity of haplotype inference with an entropy measure, where a completely unambiguous inference has entropy 0. We propose this measure as an alternate, pragmatic criterion for haplotype blocks, that predicts the effectiveness of inferred haplotypes for association analysis. Since no coalescence calculation is used (being incompatible with the EM algorithm) this inferability is purely a function of the genotype/haplotype combinatorics, and the underlying LD. In general, inferability is a function of higher order LD; we conjecture that a necessary and sufficient condition for unambiguous haplotype inference for n loci is that D’ = 1 for LD of some order ≤ n. 

A limitation is that the EM algorithm may give an unambiguous prediction even where the inference is uncertain. To catch such errors, `HapFerret` contains a bootstrapping procedure; comparing inference between successive bootstrap replicates catches some cases of spurious precision. We show with data with known haplotypes that bootstrapped inference from HapFerret has an accuracy intermediate between standard EM and PHASE. HapFerret may be downloaded here.

### Getting Started

##### USAGE: `./hapferret`
`HapFerret` takes no command line input. It extracts the information it needs at runtime from two files with static names:
- `files_data.txt`
- `hap_search_settings.txt`

The significance and structure of these files are described below. The names of the input and output files are fixed; to run it you need to create a folder and put the input files—and a copy of ferret—there. The output files will be created there also.  So the strategy is to either create a new folder for each run—easy enough on a Mac—or to keep one folder but replace the input and move or rename the output after each run.    

#### Step 0.) Download HapFerret
Open a terminal window, and clone this repository from Github:
```bash
$ git clone https://github.com/CCBR/HapFerret.git
```
Navigate to the folder containing ferret:
```bash
$ cd HapFerret/
```
Take a peak around at the contents of the repository:
```bash
$ ls -larth
```
#### Step 1.) **Setting up HapFerret**
Now, let's create a new directory for your analysis and copy over the required files:
```bash
# Create a new directory for analysis 
$ cd .. && mkdir test_analysis
# Copy over HapFerret required files (executable and input files)
$ cp HapFerret/hapferret test_analysis/
$ touch test_analysis/files_data.txt test_analysis/hap_search_settings.txt
```
#### **Step 2.) Understanding the required input files**
The **required input files** are `files_data.txt` and `hap_search_settings.txt`. **Note**: If your genotype information file is in long-format, there is one additional required file realted to locus information. Here is more information about each of these files below.

#####	Preface
Input files and their formats are specified in a file with default (currently required) name `files_data.txt`, which must be in the directory from which HapFerret is run. An example of each data_files for each files format (wide and long) is provided below.

**File Formats:** HapFerret accepts two genotype input formats: `wide` and `long`.   

##### I.) Input data and data file: 
- `genotype data`: can be in long or wide format
- `files_data.txt`: the contents of this file change if the genotype data is in wide or long format

###### A.) **Wide format**: all genotypes for an individual on one line
`Wide` format has all genotype information for one individual on one line. Each line starts with a subject identifier, and has one genotype entry for each locus. This file must have a title line; the first entry, e.g. “PID”, identifies the column of subject IDs and is ignored.  This is followed by the names of each locus. This format is specified by fileformat “c” in the files_data.txt file.

**Example: wide-genotype input file** `hla_ceu_haps_14long.txt`
```bash
$ cat hla_ceu_haps_14long.txt
```
```
PID rs136160 rs136161 rs713753 rs4419330 rs4350853 rs136168 rs2239785
UX19193     G,C     G,G     C,C     T,C     T,T     A,A     A,G
UX18501     G,C     C,G     T,C     T,T     T,T     A,G     A,G
UX19093     G,C     C,G     C,C     T,C     T,T     A,A     A,G
UX19209     G,G     G,G     C,C     T,T     T,T     A,G     A,G
UX19144     G,C     C,G     C,C     T,T     T,T     A,G     G,G
UX19222     G,G     C,G     T,C     T,T     T,T     A,G     A,G
UX19193     G,C     G,G     C,C     T,T     T,T     A,G     A,G
UX19101     C,C     G,G     C,C     T,T     T,T     G,G     G,G
UX19101     C,C     G,G     C,C     T,C     T,T     G,A     G,G
UX19127     C,C     G,G     C,C     T,T     T,T     G,G     G,G
UX19140     C,C     G,G     C,C     T,T     T,T     G,G     G,G
UX19209     G,C     G,G     C,C     T,C     T,T     A,G     A,G
UX18522     C,G     G,G     C,C     T,C     T,G     G,A     G,A
UX19138     C,C     G,G     C,C     T,T     T,T     G,G     G,G
```
**Note:** These are standard text files, with entries divided by white space, one or more spaces or tabs.  Note that the loci are given by the standard (rs #) name, and the genotypes are `allele,allele`, with the allele given as the actual base. The locus and the allele can be identified by any alphanumeric (without spaces, and without spaces between the allele names and the comma).

**Example: wide-genotype data file** `files_data.txt`
```bash
$ cat files_data.txt
```
```
fileformat c
genotype_file hla_ceu_haps_14long.txt
```

*Description:* Each line starts with an identifier, as shown. The two lines give the file format (here “c” for comma, specifying allele, allele format for each genotype), then the genotype file name. Please see an example wide-format `files_data.txt` above.

###### B.) **Long format** (one genotype—*one locus for one individual*—on each line)
`Long` has one one genotype on each line:  PID, race, polymorphism identifier, and the two alleles seen. Entries are separated by white space; here the two alleles are separated by white space and not by a comma. This is specified by fileformat “s” in the files_data.txt file.  

**Example: long-genotype input file** `c22mge_ceu.txt`
```bash
$ cat c22mge_ceu.txt
```
```
44O00768   1   rs3752462   C   C
44O00769   1   rs3752462   C   T
44O00774   1   rs3752462   T   T
44O00779   1   rs3752462   C   C
44O00781   1   rs3752462   C   T
44O00782   1   rs3752462   C   C
44O00785   1   rs3752462   C   T
44O00787   1   rs3752462   C   C
44O00793   1   rs3752462   C   T
44O00794   1   rs3752462   C   T
44O00797   1   rs3752462   C   T
44O00803   1   rs3752462   T   T
44O00813   1   rs3752462   C   T
44O00825   1   rs3752462   C   T
44O00837   1   rs3752462   C   T
44O00840   1   rs3752462   C   T
44O00841   1   rs3752462   C   T
44O00844   1   rs3752462   C   C
44O00846   1   rs3752462   C   T
44O00849   1   rs3752462   C   T
44O00852   1   rs3752462   T   T
44O00853   1   rs3752462   C   C
44O00855   1   rs3752462   C   C
44O00859   1   rs3752462   C   C
44O00860   1   rs3752462   C   C
44O00862   1   rs3752462   C   C
44O00868   1   rs3752462   C   C
44O00869   1   rs3752462   C   C
44O00872   1   rs3752462   C   T

```
*Description*: Each line is one genotype:  PID, race, polymorphism identifier, and the two alleles seen; entries separated by white space—here the two alleles are separated by white space and not by a comma.

**Example: long-genotype data file** `files_data.txt`
```bash
$ cat files_data.txt 
```
```
fileformat s
genotype_file c22mge_ceu.txt
var_info_file c22mge_ceu.i.txt
```
*Description:* Again, each line starts with an identifier. The three lines give the file format (“s” stands for scan input) then the genotype file name, then the associated locus info file. Please see an example long-format `files_data.txt` above.

**Note**: Long format requires a locus information file (see `var_info_file` above), with each line giving a polymorphism id, followed by its genome coordinates (no chromosome info is used in the current version).  The current algorithm needs this information to know the order of the loci, but doesn’t use the actual position for the calculation. Here is an example of locus information file (see below):

**Example: locus information file** `c22mge_ceu.i.txt`

```bash
$ cat c22mge_ceu.i.txt
```
```
rs1557529 35035474
rs2157256 35037606
rs2413396 35038030
rs5750250 35038428
rs3830104 35038569
rs4820229 35038699
rs3752462 35040128
rs8141971 35041308
rs5756152 35042417
rs9610489 35043476
rs2239784 35044580
rs1005570 35045219
rs12159211 35049108
rs8136336 35052479
rs16996672 35055916
rs16996677 35057228
rs11704382 35058098
rs4820234 35059020
````
##### II.) **Settings file**
- `hap_search_settings.txt`: specify parameters for haplotype inference

This last required file provides extra parameters used for the haplotype inference. At the current moment, the filename is fixed, so the program expects a file called `hap_search_settings.txt`.

**Example: Settings file** `hap_search_settings.txt`
```bash
$ cat hap_search_settings.txt
```
```
accept_params 1
accept_params 1
mode 1
blocksequence 0
race 1
disease_data 0
target_delta 0.00003
max_iterations 1000
full_hap_call 1
subseq_hap_call 0
max_subblock = 30
calc_hap_call_entropy 1
n_bootstrap_reps 0
```
**Note**: These key, value pairs do not need to be in this order. Here is more information about each parameters:

> **accept_params** (1 or 0)  Should these parameters be used as input (1), or should the user queried to enter possible changes
> **mode** 1  (keep this setting, ignore this)
> **blocksequence** 0 (keep this setting, ignore this)
> **race** Used if genotype input file is long format, then this specifies which race—identifier in 2nd column of file—to use.  Ignored in wide format.
> **disease_data** (0 or 1)  1 if there is a disease data file; not documentation yet
> **target_delta**  The EM algorithm iterates and sums the frequency difference between iterations for each haplotype, stopping when this difference is less than this parameter
> **max_iterations** Stop the inference at this count of iterations even if target_delta hasn’t been reached
> **full_hap_call**  Infer haplotypes for the complete (ordered) set of loci in the input file
> **subseq_hap_call**  Infer haplotypes for subsequences of the input loci.  hapferret  runs along the set of loci, inferring haplotypes first for all sequences of two (contiguous) loci, then for sequences of 3, 4, etc. loci, to the limit given by max_subblock 
> **max_subblock**  Maximum length of sequences of loci to infer in subseq_hap_call.  
> **calc_hap_call_entropy** Calculated entropy—uncertainty—of the inference
> **n_bootstrap_reps**  Number of bootstrap replications, set to 0 for no bootstrapping—keep this setting for now.

### Running HapFerret
`HapFerret` takes no command line input. It extracts the information it needs at runtime from two files with static names:
- `files_data.txt`
- `hap_search_settings.txt`  

It can be run as follows:
```bash
# HapFerret Expects files_data.txt and hap_search_settings.txt 
# to be in the same working directory as the hapferret executable
$ ./hapferret
```

**Note**: If you genotype data is in long format, you will need to provide one additional file containing locus information (please see above). 


