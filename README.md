# ZOLTAR

https://francoiskroll.shinyapps.io/zoltar/

Shiny app for predictive pharmacology

Francois Kroll, 2022 @Rihel lab, UCL.

Please cite us if you use ZOLTAR in your work.

François Kroll, Joshua Donnelly, Güliz Gürel Özcan, Eirinn Mackay, Jason Rihel  
**Behavioural pharmacology predicts disrupted signalling pathways and candidate therapeutics from zebrafish mutants of Alzheimer’s disease risk genes**  
_eLife_, 2024  
https://doi.org/10.7554/eLife.96839.3

If you use a dataset included here, please also cite:

Jason Rihel, David A Prober, Anthony Arvanites, Kelvin Lam, Steven Zimmerman, Sumin Jang, Stephen J Haggarty, David Kokel, Lee L Rubin, Randall T Peterson, Alexander F Schier  
**Zebrafish behavioral profiling links drugs to biological targets and rest/wake regulation**  
_Science_, 2010  
https://10.1126/science.1183090

Thank you!

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

## About the dataset

Feel free to make use of the ZOLTAR dataset in your own work! Here are some pointers.

`drugDb.csv` is the main dataset with the behavioural data. Columns are:
* `name` original name of the compound from Rihel et al., 2010
* `cleanm` simplified name where most of the stereochemistry and salt information was removed. This allowed to assign the same PubChem CID (see below) to compounds that are either the same but recorded differently (e.g. ‘Dopamine hydrochloride’ and ‘Dopamine HCl) and to compounds that varied by their stereochemistry or salts (e.g. atropine sulfate vs. atropine methyl nitrate, `cleanm` is "atropine" for both)
* `cid` PubChem CID
* `tid` Therapeutic Target Database ID
* `library` source drug library, column from Rihel et al. 2010
* `concentration` concentration used for the experiment. We relied on the source library to assign the concentration.
* `MolecularWeight` in g/mol
* `Complexity` unsure what this column represents, some version of a complexity measure, perhaps based on Tanimoto score. One should re-calculate it if it is useful. It is not used in ZOLTAR.
* `inSmall` whether this compound was labelled as "shortlisted" by Rihel et al., 2010. A compound was shortlisted if it affected at least one behavioural parameter with a large effect size and/or affected the same parameter in the same direction across the two days/nights.
* `structCluster` as "Complexity", unsure. I think it groups compounds with similar structures. Will need to search Rihel et al., 2010 in greater details or re-calculate. It is not used in ZOLTAR.
* `pharma` pharmaceutical activity/drug targets labelled by Rihel et al., 2010. Those annotations are not used by ZOLTAR. Instead, we extracted drug targets from Therapeutic Target Database (TTD) and STITCH, we would recommend using these newer/more detailed annotations.

Following columns are the datapoints of each behavioural fingerprint. Units are in z-scores compared to controls.

Each behavioural parameter has four datapoints: night1, day1, night2, day2.

The behavioural parameters are:
* `sleep` total sleep (originally in hours)
* `sleepBout` number of sleep bouts
* `sleepLength` average sleep bout length (originally in minutes)
* `sleepLatency` minutes until first sleep bout
* `averageActivity` seconds active by minute
* `averageWaking` seconds active by minute, excluding sleep (inactive minutes)

The last two columns `nightmean_averageWaking` and `daymean_averageWaking` are means of day1 & day2 or night1 & night2 averageWaking. They are not used by ZOLTAR as they are a bit redundant.

`drugDbSUM.csv` is a simplified version of this database where each compound (unique PubChem CID, so ignoring different stereochemistry or salt) is present as a single average fingerprint. ZOLTAR uses this dataset when calculating enrichments so that a single compound could not drive a significant enrichment simply because it was present as many fingerprints. This dataset could also be used to get a more representative unique fingerprint for each compound. Column `nfingerprints` records how many fingerprints (rows of `drugDb.csv`) were averaged. For columns `library`, `concentration`, `MolecularWeight`, `Complexity`, `structCluster`, `pharma`, the different rows are concatenated with '/'.

`compounds.csv` has the same columns `name`, `cleanm`, `cid`, `tid`, `MolecularWeight`, `Complexity` as `drugDb.csv` but simply lists only once each unique compound (unique PubChem CID).

`TTDindications.csv` information from TTD about the indications for each compound (i.e. what the compound is used for). A compound can be present multiple times if it has multiple indications. Only compounds with a TTD ID are listed so the total count is much lower than the number of unique PubChem CID in `compounds.csv`.

`TTDtargets.csv` information from TTD about the drug targets for each compound (i.e. the protein with which the compound interacts with). A compound can be present multiple times if it has multiple targets. `TARGETID` is an ID for each target protein from TTD.

`TTDkegg.csv` information from TTD about the KEGG pathways each compound interacts with through each of its target. Each row is a unique TTD ID (compound) – target – KEGG pathway link.

In the TTD datasets, only compounds with a TTD ID are listed so the total count of compounds is much lower than the number of unique PubChem CID in `compounds.csv`.

## Version history

### v1
First online version.

### v2
Two changes to the analysis:
* When calculating enrichments, ZOLTAR uses a version of the database (drugDbSUM.csv) where each drug is present as only one fingerprint. If a drug had multiple replicate fingerprint, they were averaged. This is to avoid having significant enrichments that are simply driven by multiple fingerprints of the same drug. Additional support for an hypothesis should come from more unique drugs with the same annotation having high or low cosines, not from a single drug with a particularly high or low cosine being in the database many times. Tab "Drug fingerprints ranked" still lists the complete database, where a drug can be present as multiple fingerprints.
* ZOLTAR now directly sums (absolute) cosines, not ranks, to measure enrichments. Previously, in a situation where there were no high cosines for a given query fingerprint (e.g. maximum cos = 0.3) and a situation where there were high cosines (e.g. maximum cos = 0.8), the drug with the maximum cosine (0.3 or 0.8) would both get the same rank, so would be worth the same in the enrichment analysis. I now think it makes more sense to use directly the cosine, as the drug with cosine 0.8 should provide more support for the hypothesis (enrichment of its annotations) than the drug with cosine 0.3.
  
There are other aesthetics/formatting changes. For example, the windows that appear when clicking on a row in the tables is bigger. There are also more detailed descriptions of each tab.