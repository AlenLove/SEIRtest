# MVP [![](https://img.shields.io/badge/Issues-%2B-brightgreen.svg)](https://github.com/XiaoleiLiuBio/rMVP/issues/new) [![](http://www.r-pkg.org/badges/version/rMVP?color=red)](https://cran.r-project.org/package=rMVP) [![](https://img.shields.io/badge/GitHub-1.0.6-blueviolet.svg)]() ![](http://cranlogs.r-pkg.org/badges/grand-total/rMVP?color=green) [![](https://cranlogs.r-pkg.org/badges/rMVP)](https://cran.r-project.org/package=rMVP) <a href="https://hits.seeyoufarm.com"/><img src="https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2Fxiaolei-lab%2FrMVP"/></a>

## A [M](https://)emory-efficient, [V](https://)isualization-enhanced, and [P](https://)arallel-accelerated Tool for Genome-Wide Association Study

## SEIR: Selector-Embedded Iterative Regression
   A fast and extremely powerful R package for GWAS
  
   Written by Mengjin Zhu & GuangLiang Zhou
  
   Last update: Sep 27, 2021
## Installation
   install.packages("devtools")

   devtools::install_github("AlenLove/SEIRtest")
## Demo Data
   SEIR R version only support numeric data type
## Usage
   library(SEIR)

#### # genotype information data
    data(map)
#### # genotype data
    data(geno)
#### # phenotype data
    data(phe)
#### # association analysis
SEIR(Y=phe,X=geno,GM=map,CV=NULL,maxStep=10,
    selector="stepwise",fun.selection="fastLm",
    extraction="min",X.check="FALSE",chunk.num = NULL,
    file.output=FALSE,plot.style="SEIR",cutOff=0.05)
#### # more features   
More parameters explained [here]
### Issues
Bugs could be reported [here]
## The license notice
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details
## Author
Mengjin zhu (glzhou@webmail.hzau.edu.cn)
