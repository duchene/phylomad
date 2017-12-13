![software logo](codeFolder/www/phylomad.temp.png)


## Introduction

This repository contains PhyloMAd, a software for easily-accessible assessment of phylogenetic model adequacy. Briefly, an assessment of model adequacy provides an independent account of the power of a model for describing the empirical data. This is different to model selection, where a set of candidate models are compared to each other, under the assumption that at least the best-fitting model is an adequate description of the data.

Substitution models supported by PhyloMAd inlcude the JC, HKY, and GTR models, including gamma-distributed rates across sites with 4 discrete categories. Clock models supported include those implemented in BEAST 2.

## Download and package installation

PhyloMAd requires that the R statistical computing language is installed. R is freely available from the R project website.

https://www.r-project.org/

Download and unzip the PhyloMAd repository by clicking the following link:

https://github.com/duchene/phylomad/archive/master.zip

This can also be done by pressing the *Clone or download* button in github or by typing the following in a bash shell. The later option assumes that the machine has git installed.

```coffee
git clone https://github.com/duchene/phylomad.git
```

Double-click the PhyloMAd file according to the platform (runMac.command or runWin.vbs). Opening the package for the first time might take several minutes and will require internet connection if any R dependencies are not yet installed.

If you have difficulty opening the program, you might want to try opening a bash shell, setting your directory to the PhyloMAd folder and executing the R script by hand.

```coffee
cd pathToPhyloMAd
Rscript phylomad.Rscript
```

In mac machines, PhyloMAd will open a terminal window when opened and log the progress. In windows, a log file with the progress will be saved in the main PhyloMAd folder. This log file is mainly useful for debugging or monitoring progress, so it is safe to delete it when the program is closed.

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

Next is a file with an interpretation of the test for five test statistics that are known to be sensitive to multiple sources of wrong inferences of tree topology and branch lengths. These recommendations come from a detailed previous simulations study. The test allows for interpretation of the results, according to the thresholds investigated in that study. 

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
