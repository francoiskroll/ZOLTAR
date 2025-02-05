# ZOLTAR

https://francoiskroll.shinyapps.io/zoltar/

Shiny app for predictive pharmacology

Francois Kroll, 2022 @Rihel lab, UCL.

Please cite us if you use ZOLTAR in your work.

François Kroll, Joshua Donnelly, Güliz Gürel Özcan, Eirinn Mackay, Jason Rihel  
_eLife_  
**Behavioural pharmacology predicts disrupted signalling pathways and candidate therapeutics from zebrafish mutants of Alzheimer’s disease risk genes**  
https://doi.org/10.7554/eLife.96839.2

:email: francois@kroll.be  
[@francoiskroll.bsky.social](https://bsky.app/profile/francoiskroll.bsky.social)

___

## Tutorial

> Just want to try the app? You can jump to step 2• by downloading the sample data available in the [in the app](https://francoiskroll.shinyapps.io/zoltar/).

### 0• Do your experiment

Currently, ZOLTAR expects at least:
* day0 (incomplete)
* night0
* day1
* night1
* day2
* night2  
* day3 (incomplete)

Each day being 14 hr and each night being 10 hr. ZOLTAR will analyse day1, day2, and night1, night2. Night0 is excluded as habituation period. The experiment can be longer but additional days/nights are not analysed.

Rihel et al., 2010 started treatment with small molecules in the morning of 4 dpf, then behaviour tracking at 11 PM evening of 4 dpf. We prefer to start the experiment at 5 dpf as it gives a bit more time to any 'late' larvae to inflate its swim bladder. In our case, the experiment thus stops in the morning of 8 dpf. No feeding is required if not extending the experiment past that point.  

### 1• Prepare your data in 'middur' format

> Do you already have a DATA.txt file created with a MATLAB script called perl_batch_192.m? ZOLTAR will also accept this file. If you can, I think preparing the middur file from the RAWs.csv as described below is slightly better. If you choose to use DATA.txt files, assume I mean DATA.txt when I mention middur.csv below.

Preparing your data in the correct format requires the [FramebyFrame R package](https://github.com/francoiskroll/FramebyFrame).  

You will need the RAWs.csv file prepared using `vpSorter(...)`. Please follow the instructions there.  

This RAWs.csv file stores all the frame-by-frame Δ pixel data of your experiment. The Rihel et al., 2010 database of behavioural fingerprints uses minute-by-minute 'middur' data, so we need to match this format. Middur is a parameter defined by Viewpoint. It is, for each larva and each time bin (here, 1 min), the time spent above _Freeze_ threshold (typically set at 3 Δ px) and below _Burst_ threshold (typically set at 200 Δ px). If you already have the minute-by-minute middur data from your ZebraBox—sorry!—it is easier to have a common pipeline for everyone so you will still have to follow this.  

Once your YYMMDD_BX_RAWs.csv file is ready, convert it to middur data using:  

```
rawToMiddur(ffpath='full/path/to/RAWs.csv',
            freezing=3,
            burst=200,
            exportOrNo=TRUE)
```

Adjust `freezing` and `burst` to the values you usually use in ZebraLab.  

This will generate YYMMDD_BX_middur.csv in the ffpath folder.  

You will also need a genotype.txt file to label each well with a condition. The `genotypeGenerator(...)` from the FramebyFrame package may be helpful. 

Your genotype file should be called YYMMDD_BXgenotype.txt. It should match the corresponding middur.csv file. For example, the genotype file for experiment 230214_04_middur.csv must be called 230214_04genotype.txt.  

### 2• Launch the ZOLTAR app

https://francoiskroll.shinyapps.io/zoltar/

### 3• Give it your behavioural data

Drag-and-drop or browse to input your middur.csv file(s) and genotype file(s).  

If you have replicate experiments, you can input multiple middur.csv files. ZOLTAR will match each middur.csv file with its genotype file based on the YYMMDD_BX, so make sure they match (see above). ZOLTAR will then use the average fingerprint of replicate experiments for predictions.  

ZOLTAR will read the group names from your genotype files. Tell it which group is the treatment group and which is the control group.

### 4• Press Go!

You should see your results gradually appear in a total of ~ 7 min.

## Version history

### v1
First online version.

### v2
Two changes to the analysis:
* When calculating enrichments, ZOLTAR uses a version of the database (drugDbSUM.csv) where each drug is present as only one fingerprint. If a drug had multiple replicate fingerprint, they were averaged. This is to avoid having significant enrichments that are simply driven by multiple fingerprints of the same drug. Additional support for an hypothesis should come from more unique drugs with the same annotation having high or low cosines, not from a single drug with a particularly high or low cosine being in the database many times. Tab "Drug fingerprints ranked" still lists the complete database, where a drug can be present as multiple fingerprints.
* ZOLTAR now directly sums (absolute) cosines, not ranks, to measure enrichments. Previously, in a situation where there were no high cosines for a given query fingerprint (e.g. maximum cos = 0.3) and a situation where there were high cosines (e.g. maximum cos = 0.8), the drug with the maximum cosine (0.3 or 0.8) would both get the same rank, so would be worth the same in the enrichment analysis. I now think it makes more sense to use directly the cosine, as the drug with cosine 0.8 should provide more support for the hypothesis (enrichment of its annotations) than the drug with cosine 0.3.
  
There are other aesthetics/formatting changes. For example, the windows that appear when clicking on a row in the tables is bigger. There are also more detailed descriptions of each tab.