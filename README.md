# DNAm-qc-scores

Here we present three new scores for evaluating the quality of DNA methylation raw data from Illumina-based arrays. 


## DB-Score
In order to identify samples with a doubtful quality we mathematically summarized the overall DNA methylation distribution in the so called Distribution-Score. Therefore, the number of CpGs with beta values between 0.3 and 0.7 is divided by the number of CpGs with beta values smaller or equal 0.3 and higher or equal 0.7.

>db.score(input)
>- input: Name of the table with the TargetID and beta values (samples are listed per column)
>
### Example
db.score(input = data)



## CNV-Score
In order to differentiate between those cases which technical failed from those with possible meaningful biological effects we implemented another quality score named Copy Number Variation Score (CNV-Score). As basis for the calculation of the CNV-Score the R-package Conumee is used, which enables copy number variation calling based on DNA methylation data. Thereby, CpGs within a predefined range are summarized in bins (represented as points in the plot) and used to visualize gains and losses throughout the whole genome. Within these plots the distribution of the bins over the y-scale can be used as an indicator for a technical failure and therefore differentiate between those cases with a possible meaningful biological background

>cnv.score(input = Mset, array_type = "450k")
>-input: A Mset containing methylated and unmethylated signals (preferably generated with the minfi >package)
>-array_type: Choose "EPIC" or "450k"

### Example
cnv.score(input = Mset, array_type = "450k")
