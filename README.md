# DNAm-qc-scores

Here we present three new quality control scores (QC-Scores) for evaluating the quality of DNA methylation raw data from Illumina-based arrays.  We suggest following cut-offs for the QC-Scores, which are based on experience.  


| QC-Score |     Cut-off          | Interpretation |
|:---------	|:------------------	|:---------------------|
|DB-Score 	| <1 		| A DB-Score above 1 indicates a deviation of the bimodal distribution of beta values and thus a technical problem or a possible meaningful biological background. |
| CNV-Score	| <0.25 		| A CNV-Score above 0.25 is an indication of a technical error or a problem with the DNA itself. |
| HL-Score high and low | <20% | Describe the percentage of deviant CpGs within the stable methylated loci and indicates a technical error or a problem with the DNA itself. |
| HL-Score difference | <20% | Describe the differences between the HL-Scores high and low. A high difference is a measurement of a possible meaningful background. |


## DB-Score
In order to identify samples with a doubtful quality we mathematically summarized the overall DNA methylation distribution in the so called Distribution-Score. Therefore, the number of CpGs with beta values between 0.3 and 0.7 is divided by the number of CpGs with beta values smaller or equal 0.3 and higher or equal 0.7. The DB-Score can be calculated using the db.score() function.

![Examples Histograms and DB-Score](DB.Score/Figure_DB-Score_Histogram.png)
**Figure Histograms:** Outstanding examples of samples with good (left), doubtful (middle) and bad (right) quality based on their DB-Scores.


>db.score(input)
>- input: Name of the table including one column with the TargetID followed by the beta values (samples are listed per column)
>
### Example
db.score(input = data)



## CNV-Score
To differentiate between those cases which technical failed from those with possible meaningful biological effects we implemented another quality score named Copy Number Variation Score (CNV-Score). As basis for the calculation of the CNV-Score the R-package conumee is used, which enables copy number variation calling based on DNA methylation data. Thereby, CpGs within a predefined range are summarized in bins (represented as points in the plot) and used to visualize gains and losses throughout the whole genome. Within these plots the distribution of the bins over the y-scale can be used as an indicator for a technical failure and therefore differentiate between those cases with a possible meaningful biological background. The CNV-Score is calculated using the cnv.score() function. Tables including the ranges of the bins for EPIC and 450k are needed for calculation and are provided in the corresponding folder.


![Outstanding examples of samples with good (A) and bad (B) CNV-Scores](CNV.Score/Figure_CNV-Plots.png)
**Figure CNV-Plots:** Outstanding examples of samples with good (A) and bad (B) CNV-Scores.

>cnv.score(input, array_type)
>
>- input: A Mset containing methylated and unmethylated signals (preferably generated with the minfi package)
>- array_type: Choose "EPIC" or "450k"

### Example
cnv.score(input = Mset, array_type = "450k")


## HL-Scores
As a further indication for a possible meaningful biological background the Heatmap-Lane-Scores (HL-Scores) were established. Therefore, CpG loci with a constant DNA methylation pattern over various tissues, cancers and preparations methods for samples were identified. This CpG loci in turn were further differentiated in stable hypermethylated and hypomethylated loci (450K: 279 hyper- and 313 hypomethylated, EPIC: 249 hyper- and 299 hypomethylated). Based on this stable hyper- or hypomethylated loci three HL-Scores were implemented: HL-Score high, HL-Score low and HL-Score difference (absolute difference between HL-Score high and low). All three scores are calculated using the hl.score() function. A list including the stable loci is needed for the calculation and provided in the corresponding folder. 

> hl.score(input)
> 
>- input: Name of the table including one column with the TargetID followed by the beta values (samples are listed per column)
>
### Example
hl.score(input = data)
