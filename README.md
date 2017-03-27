**************************************************
    Author: Bárbara D Bitarelo

    Created: 27.03.2015

    Last modified: 27.03.2017

**************************************************

Welcome to the NCD repo! 

This repository provided scripts related to the manuscript: http://biorxiv.org/content/early/2017/03/22/119529


**What are the  NCD (Non-Central Deviation) statistics?**

![Fig1](Figures_main/Fig1_red.tiff)

NCD statistics measure the average difference between allele frequencies in a given region from a deviation point, which we call the 'target frequency (*tf*)'. So, assuming *tf* = 0.5, the more the allele frequencies are close to 0.5, the lower the NCD values. We propose two implementations of this statistic: *NCD1* (only SNPs are required and used as informative sites) and *NCD2* (SNPs and fixed differences are used as informative sites). In the manuscript, *NCD2* as used to scan the human genome.

*******************************************************



Here, we show how to:
* item  Run *NCD1* and *NCD2* (ongoing)
* item  Using examples from the manuscript (see above) which can be extended to other species

*************************************************************************


Running NCD:

this requires:
	
	* SNP input data (file in modified VCF format, see below an example)

	* Fixed differences (FD) input data (e.g. human-chimp FD bed file, as used in the manuscript)
	Note: *NCD1* only requires the first input file, whereas NCD2 requires both.

	* SGE script for parallelizing (option A)  **OR** parallelizing option without sge script.(option B)

	* Coming up: an optimized R script that does not require parallelizing jobs (this will be called option C)

Note: although NCD can be run withour parallelizing, that takes quite some time. If you have a cluster or a supercomputer, it is best to use it (option A, below). If you don't have this, you will use option B.In any case, you need first to download the input data. In the near future, we want to provide a better R function to calculate NCD that does not require parallelizing.

*************************************

Example input data (SNP data):


| CHROM | POS | ID | REF | ALT | Anc | AWS | LWK | YRI | CEU | FIN | GBR | TSI | CHB | CHS | JPT | MXL | CLM | PUR |
| ----- | --- | -- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|   1   | 834088  |  rs192235345  |  T  |  C  |  .  |  0  |  0  |  0  |  0  |  0  |  0  |  0  |  0  |  1  |  0  |  0  |  0  |  0  |

Note: a longer version of this is in example_input_files/

Example FD file

to be continued...



************************************************************************
 Getting Started: Options A and B
************************************************************************

First:

* clone this repo: go to your directory and clone:

```
git clone https://github.com/bbitarello/NCD-statistics.git
```

* go to root NCD directory

```
cd NCD-statistics/
chmod 777 * #this will avoide several issues. Make sure you are allowed to do this.
```

**Important Note**: Please note the paths for the files, logs, and tmp directories and addapt them as needed. Only 4 files need their paths to be edited:

What you should do? Replace all intances of /mnt/sequencedb/PopGen/barbara/NCD-Statistics/ by the path of the directory where you cloned this repo.

#1.Edit the paths in:
 
```
vim run_NCV_sge.sh #logs, tmpdir, input files (lines 18,21,22,23,26)
```
#2. Edit the paths in:

```
vim scripts/NCV_dir_package/scripts/run_ncv_allpops_Rscript.sge #Rscript file path

```

#3.Edit the paths in:

```
vim scripts/run_ncv_allpops_Rscript_v1.r #loading the NCV funciont in line 63
```

#4. Edit the paths in:

```
vim scripts/run_ncv_allpops_Rscript_nSGE.r #loading the NCV funciont in line 63
```

**IMPORTANT NOTE**: if you use this option, make sure you go to the file which ends with .sge in the scripts folder and edit it according to your usual SGE settings. Also, replace my email by yours.


############################################################
#
##### Option A ##### Option A ##### Option A ##### Option A
#
############################################################
#This option applies if and only if you use a SGE for job submission.
#In this case, do:

#next, cd to NCV_dir_package and run:
./run_NCV_sge.sh

#This script calls the .sge file in the /scripts folder, and that file allow submitting severaljobs each with 3Mb of sequence data to the cluster.

#These scripts are optimized to run NCV in ~900 parallel jobs at the MPI-MPG in Leipzig.
#############################################################
#
###### Option B ##### Option B ##### Option A ##### Option B
#
#############################################################
#Obtion B) No SGE, just a plain old cluster without job submission management:
#
#If you don't use SGE or anything of the sort, but are using a linux machine (of course):


#First:
#cd to 'NCV_dir_package':

#open an R session and type

install.packages('getopt')
install.packages('littler')

#close the session. Now, run:

./run_NCV.sh

#this will parallelize the jobs using 'GNU parallel'. It is not a very clever solution, but it works nevertheless.


