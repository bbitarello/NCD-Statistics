*************
Input Files
*************


This data is only for example purposes. The SNP data is fro Phase 3 1000 Genomes. We made this example input by taking the first 300,000 lines of each chromosome (for both SNP and FD files), so the results from this example do not reproduce those of the paper.

**SNP File**

The VCF format is pretty much standard these days. The (NCD) scripts are suitable, howver, for an R object in data.table format, i.e, it assumes you have already handled reading in vcf into R and pre-processing it to meet the requirements of the NCD function.
Note: I do intend to provide scripts for going from vcf to NCD input format in the near future, but for now, the use will have to do it.

In an R session:

```
load(data.table)
load('SNP_test_input.RData')
LWK[[1]]
```

The above command should show you the first and last five lines of the SNP input file for the LWK population. Analogously, LWK[[22]] contains the data for chr 22, etc. I.e, the object 'LWK' is a list of 22 elements, where each element is a chromosome. 
The output of the above command should look like this:



| CHR | POS | ID  | REF | ALT | AF1 | AF2 | AF3 | MAF |
| --- | --- | --  | --- | --- | --- | --- | --- | --- |
|1 | 10505  | 1\|10505 | A | T | 1.000000 | 0.000000 | NA | 0.000000 |
|1 | 10506  | 1\|10506 | C | G | 1.000000 | 0.000000 | NA | 0.000000 |
|1 | 10511  | 1\|10511 | G | A | 1.000000 | 0.000000 | NA | 0.000000 |
|1 | 10539  | 1\|10539 | C | A | 1.000000 | 0.000000 | NA | 0.000000 |
|1 | 10542  | 1\|10542 | C | T | 1.000000 | 0.000000 | NA | 0.000000 |
---
|1 | 9682490 | 1\|9682490 | A | T | 1.000000 | 0.0000000 | NA | 0.0000000 |
|1 | 9682512 | 1\|9682512 | C | T | 1.000000 | 0.0000000 | NA | 0.0000000 |
|1 | 9682534 | 1\|9682534 | C | T | 0.934343 | 0.0656566 | NA | 0.0656566 |
|1 | 9682555 | 1\|9682555 | G | A | 0.989899 | 0.0101010 | NA | 0.0101010 |
|1 | 9682631 | 1\|9682631 | G | A | 0.525253 | 0.4747470 | NA | 0.4747470 |

Where CHR and POS define the chromosomal position of the SNP, ID is an ID for the position, defined based on the first two columns, REF is the human reference allele (hg19), ALT is the alterante allele, AF1:AF3 are the allele frequencies of the alleles, MAF is the minor allele frequency.

Note: there are actually a few additional columns in this example input file, but they can be excluded or ignored.

**FD Input Data**

Remember: if you only want *NCD1* you shall not need this, ever. However, unless there is really no closely-related outgroup data available, we strongly recommend *NCD2* instead (see paper for more on that).

In an R session:


```
load(data.table)
load('FD_test_input.RData')
FD_list[[1]]
```

The above command should show you the first and last five lines of the FD input file. Analogously, FD_list[[22]] contains the data for chr 22, etc. I.e, the object 'FD_list' is a list of 22 elements, where each element is a chromosome.
The output of the above command should look like this:


| CHR |POS | REF | Chimp_REF | ID |
| --- | -- | --- | --------- | -- |
|1 | 834174 | T  | C | 1\|834174 |
|1 | 834187 | C  | T | 1\|834187 |
|1 | 834214 | T  | C | 1\|834214 |
|1 | 834215 | G  | A | 1\|834215 |
|1 | 834240 | G  | A | 1\|834240 |
---
| 1 | 31645394 | T | C 1|31645394 |
| 1 | 31645397 | G | A 1|31645397 |
| 1 | 31645398 | A | C 1|31645398 |
| 1 | 31645413 | C | T 1|31645413 |
| 1 | 31645454 | A | G 1|31645454 |

Where CHR and POS define the chromosomal position, REF is the human reference allele (hg19), Chimp_REF is the chimp reference allel (pantro2), and ID is an ID for the position, defined based on the first two columns.

