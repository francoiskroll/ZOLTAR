# ZOLTAR

https://francoiskroll.shinyapps.io/zoltar/

Shiny app for predictive pharmacology

Francois Kroll, 2022 @Rihel lab, UCL.

[![alt text][1.2]][1] [@francois_kroll](https://twitter.com/francois_kroll)

:email: francois@kroll.be

<!-- icons with padding -->
[1.1]: http://i.imgur.com/tXSoThF.png (twitter icon with padding)

<!-- icons without padding -->
[1.2]: http://i.imgur.com/wWzX9uB.png (twitter icon without padding)

<!-- links to your social media accounts -->
[1]: https://twitter.com/francois_kroll

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

### 2• Launch the ZOLTAR app

https://francoiskroll.shinyapps.io/zoltar/

### 3• Give it your behavioural data

Drag-and-drop or browse to input your middur.csv file(s) and genotype file(s).  

If you have replicate experiments, you can input multiple middur.csv files. ZOLTAR will match each middur.csv file with its genotype file based on the YYMMDD_BX, so make sure they match. ZOLTAR will then use the average fingerprint of replicate experiments for predictions.  

ZOLTAR will read the group names from your genotype files. Tell it which group is the treatment group and which is the control group.

### 4• Press Go!

You should see your results gradually appear in a total of ~ 7 min.