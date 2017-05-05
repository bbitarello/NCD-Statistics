**************************************************
    Author: Bárbara D Bitarelo

    Created: 27.03.2015

    Last modified: 05.05.2017

    Language: R

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
	
	* SNP input data as a data.table object in R

	* Fixed differences (FD) input data (e.g. human-chimp FD bed file, as used in the manuscript) as a data.table object in R.
	Note: *NCD1* only requires the first input file, whereas NCD2 requires both
        
	* An open R session



#*************************************

Example input files are provided in example_input_files/. Refer to the README.md in that directory for further explanations.


#************************************************************************

First:

* clone this repo: go to your directory and clone:

```
git clone https://github.com/bbitarello/NCD-statistics.git
```

* go to root NCD directory

```
cd NCD-statistics/
```


Second:

* open an R session and type

```
source('scripts/pramble.R') #loads several packages
source('scripts/NCD_func.R') #loads NCD functions NCD1 and NCD2
load('example_input_files/SNP_test_input_data.RData') #necessary for NCD1 and NCD2
load('example_input_files/FD_test_input_data.RData') #only necessary for NCD1
system.time(example.run<-foreach(x=1:22, .combine="rbind", .packages=c("data.table")) %dopar% NCD1(X=LWK[[x]], W=3000, S=1500)); #  
```
Note that the runtime will vary considerably depending on computational constraints. The reported value is for 11 cores in registerDoMC(11). See example_input_files/README.md 

Important: This is an example. It is a roadmap of how NCD can be used for other input data. Even though the example input data is real, it does not reproduct the findings from the paper.

