![software logo](codeFolder/www/phylomad.temp.png)

For support contact David A. Duchêne
david.duchene[at]sydney.edu.au

23 May 2017

Introduction
------------

This repository contains PhyloMAd, a software for easily-accessible assessment of phylogenetic model adequacy.

Copyright 2017 by the PhyloMAd authors. The software PhyloMAd is distributed without warranty of any kind or support for previous versions. The authors will not be responsible for any damage resulting from the use of this software. The source and documentation are distributed under the GNU General Public Licence except where stated otherwise. See http://www.opensource.org/licenses for details.

Substitution models supported inlcude the JC, HKY, and GTR models, including gamma-distributed rates across sites with 4 discrete categories. Clock models supported include those implemented in BEAST 2.

Download and package installation
---------------------------------

PhyloMAd requires that the R statistical computing language is installed. R is freely available from the R project website.

https://www.r-project.org/

Download and unzip the PhyloMAd repository from GitHub. This can be done by pressing the *Clone or download* button above or by typing the following in a bash shell. The later option assumes that the machine has git installed.

```coffee
git clone https://github.com/duchene/modadclocks.git
```

Double-click the PhyloMAd file according to the platform (runMac.command or runWin.vbs). If needed, this will install all the required R packages. This might take several minutes.

If you have difficulty opening the program, you might want to try opening a bash shell, setting your directory to the PhyloMAd folder and executing the R script by hand.

```coffee
cd pathToPyloMAd
Rscript phylomad.Rscript
```

Once you have opened PhyloMAd, select the model you wish to assess from the buttons at the left of the screen, and press the *Install required packages* button. This will install any other packages that might not yet be installed.

Basic assessment of substitution model adequacy
-----------------------------------------------



Advice on model adequacy based on simulations
---------------------------------------------




Basic assessment of clock model adequacy
----------------------------------------



Graphical output of results
---------------------------
It is strongly recommended to use qualitative checks of models using graphical analyses. 

