![software logo](codeFolder/www/phylomad.temp.png)


## Basic features

- Assessment of phylogenetic substitution model adequacy. For assessment, PhyloMAd uses maximum likelihood phylogenetic inference (IQtree), parametric bootstrap simulations, and several common test statistics. Diagnostics of risk of biased inferences are provided based on previous extensive simulation studies.

- Assessment of substitution saturation using measures of entropy. Diagnostics of risk due to excessive saturation are provided based on an extensive simulations study.

- Assessment of the type of the strength and number of phylogenetic tree topology signals in a locus or a set of gene trees. Signal is defined as the support that appears in the data for a given species or gene tree, while the null is the signal of trees that can emerge from the model used for inference.

- Assessment of phylogenetic clock model adequacy. This assess the expectation that the number of substitutions along branches as expressed by a clock model (times and molecular rates) is consistent with the substitutions observed when using clock-free estimation.

## Introduction

This repository contains PhyloMAd, a software for assessing several components of phylogenetic analysis. Briefly, these methods of assessment provide an independent account of the power of a phylogenetic model and sequence data for describing the evolutionary process. These methods of assessment are different from model selection, where a set of candidate models are compared to each other under the assumption that at least the best-fitting model is an adequate description of the data.

PhyloMAd supports several of the most common nucleotide and amino-acid substitution models (e.g., JC, HKY, and GTR), including gamma-distributed rates across sites with 4 discrete categories. Clock models supported include those implemented in BEAST 2.

## Download and package installation

PhyloMAd requires that the R statistical computing language is installed. R is freely available from the R project website.

https://www.r-project.org/

Download and unzip the PhyloMAd repository by clicking the following link:

https://github.com/duchene/phylomad/archive/master.zip

This can also be done by pressing the *Clone or download* button in github or by typing the following in a bash shell (if available and the machine has git installed).

```coffee
git clone https://github.com/duchene/phylomad.git
```

The rest of this tutorial focuses on the Graphical User Interface. Refer to the manual for information about running PhyloMAd through the command-line.

Double-click the PhyloMAd file runMac.command. Opening the package for the first time might take several minutes and will require internet connection to install R dependencies.

If you have difficulty opening the program, you might want to try opening it directly through R. By typing the following in the R command line:

```coffee
setwd("pathToPhyloMAd")
source("phylomad.Rscript")
```

Similarly, if you have access to a bash shell, you can set your directory to the PhyloMAd folder and execute the R script by hand.

```coffee
cd pathToPhyloMAd
Rscript phylomad.Rscript
```

PhyloMAd will open a terminal window when opened and log the progress, which is useful for debugging or monitoring progress.

## Brief tutorial

This is a brief tutotrial for assessing substitution model adequacy in a single locus. Refer to the manual for more detailed settings.

### User interface

Once you have opened PhyloMAd, select the model you wish to assess from the buttons at the left of the window.

After selecting a model from the box to the left, you will see a corresponding set of tabs in the main screen ranging from *Data* to *Other options and START*.

You can go through the tabs before modifying the settings to make sure you understand them. Further information about each of the options can be found in the manual.

1. For a quick trial of the software, press *Browse* in the Data tab under the *Select nucleotide or amino acid alignment* header.

Browse through your files into the *phylomad/codeFolder/exampleData/* folder, and select the data set called *example.covarion.longtips.nex*. If you have selected the file successfully, the PhyloMAd window should say *Upload complete*.

![browseButton](codeFolder/www/browseButton.png)

2. You can now read through the other options without changing anything. Browse through each of the tabs entitled *Model*, *Test statistics*, and *Output*, making sure you understand what the settings mean, and referring to the manual for any further details.

![tabs](codeFolder/www/tabs.png)

3. Once you get to the *Other options and START* tab, you can reduce the number of simulations to 40, since this is only a test run. Check whether your machine has multiple cores and increase the number in the corresponding section to 4 cores.

Now press the *START ASSESSMENT* button. 

![startAssessment](codeFolder/www/startAssessment.png)

You can monitor the assessment in the shell (mac) or by checking the log file that is created in the main PhyloMAd folder (windows). Also check the output folder and make sure that data are being stored and removed as successive simulations are analyzed.

### Interpreting output

In this example, the output folder was kept as the default. This means that the output will have been saved in the *outputFolder*, which is in the main PhyloMAd folder.

Inside the output folder, there will be a separate folder containing the output for each locus for which the model was assessed. In the case of this example, there will be a single folder containing three files. 

![outputFiles](codeFolder/www/outputFiles.png)

One is a PDF with histograms of the values of the test statistics calculated from simulated data, with the value calculated for the empirical data set shown in a red line. These graphics are useful for qualitative interpretation of the results. This qualitative interpretation is necessary for many of the test statistics until they are better understood.

![exampleHist](codeFolder/www/exampleHist.png)

Next is a file with an interpretation of the test for five test statistics that are known to be sensitive to multiple sources of wrong inferences of tree topology and branch lengths. These recommendations come from a detailed previous simulations study. The test allows for interpretation of the results, according to the thresholds investigated in that study. Note that while the model might be adequate regarding each of these test statistics, it might still suffer an unknown type of bias. This means that the output should be placed in context of what each of these statistics describe. See manual for more details.

![chisqResults](codeFolder/www/chisqResult.png)

In this example, the chi-squared statistic test tells us that the amount of compositional heterogeneity is unlikely to be of concern.

The last file is a table that shows three values (rows) for each test statistic selected (columns). These three values include:

- The tail area probability, which is the proportion of values calculated for simulations that are lower than the value calculated for the empirical data set; 

- The statistic calculated for the empirical data set;

- The distace between the empirical data set and the mean of the simulated distribution in terms of the number of standard deviations of the distribution.

![resultsTable](codeFolder/www/resultsTable.png)

This example shows that the multinomial statistic is not available for analysis. This is because this data set comes from a scenario with extremely long terminal branches. The high divergence among taxa has led to every site having a different pattern in every data set, so every simulated data set has the same multinomial likelihood as the empirical data set. This test statistic is therefore not useful in this scenario.

## Support and bug reports

For support contact David A. Duchêne:
david.duchene[at]sydney.edu.au

## Licence

Copyright 2017 by the PhyloMAd authors. The software PhyloMAd is distributed without warranty of any kind or support for previous versions. The authors will not be responsible for any damage resulting from the use of this software. The source and documentation are distributed under the GNU General Public License except where stated otherwise. See the copy of the GNU GPL license found in the *codeFolder* for further details.
