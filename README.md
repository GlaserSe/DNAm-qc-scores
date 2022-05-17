# DNAm-qc-scores

Here we present three new scores for evaluating the quality of DNA methylation raw data from Illumina-based arrays. 


## DB-Score
In order to identify samples with a doubtful quality we mathematically summarized the overall DNA methylation distribution in the so called Distribution-Score. Therefore, the number of CpGs with beta values between 0.3 and 0.7 is divided by the number of CpGs with beta values smaller or equal 0.3 and higher or equal 0.7.

>db.score(input)
>- input: Name of the table with the TargetID and beta values (samples are listed per column)
>
# Example
db.score(input = data)
