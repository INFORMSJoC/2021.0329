[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# 2021-0329

This archive is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

The scripts in this repository are a snapshot of the scripts
that were used in the research reported on in the paper 
Convergence Analysis of Stochastic Kriging-Assisted Simulation with Random Covariates by C. Li, S. Gao and J. Du. 

## Cite

To cite this software, please cite the [paper](https://doi.org/) using its DOI and the software itself using the following DOI.

Below is the BibTex for citing this version of the code.

```
@article{li2022,
  author =        {C. Li, S. Gao and J. Du},
  publisher =     {INFORMS Journal on Computing},
  title =         {Convergence Analysis of Stochastic Kriging-Assisted Simulation with Random Covariates},
  year =          {2022},
  doi =           {unknown currently},
  url =           {https://github.com/INFORMSJoC/2021-0329},
}  
```

## Description

The goal of this repository is to illustrate the the convergence behaviors of stochastic kriging-assisted simulation with random covariates.

## Requirements

For this project, we use Matlab. The user should respect the license of the used sotfware.

## Repository Structure

### Scripts

*   **"scripts/d1"** folder:  testing the rates of maximal IMSE and IPFS under the random sampling and comparing the random design with the Adaptive MSE procedure in the **one**\-dimensional benchmark function examples,
*   **"scripts/d2"** folder:  testing the rates of maximal IMSE and IPFS under the random sampling and comparing the random design with the Adaptive MSE procedure in the **two**\-dimensional benchmark function examples,
*   **"scripts/d3"** folder:  testing the rates of maximal IMSE and IPFS under the random sampling and comparing the random design with the Adaptive MSE procedure in the **three**\-dimensional benchmark function examples,
*   **"scripts/d10"** folder: testing the rates of maximal IMSE and IPFS under the random sampling and comparing the random design with the Adaptive MSE procedure in the **ten**\-dimensional benchmark function examples.
*   **"scripts/queue"** folder: testing the rates of maximal IMSE and IPFS under the random sampling and showing how the theory can help decide the sample size to achieve a target precision in the **M/M/1 queue** example.

### Results

The **online supplement** of this paper is also in this subdirectory. The following tells the corresponding script for each result of the numerical experiment.

#### Figures in manuscript:

*   “**Figure 1**” of the manuscript (i.e., results/Figure1(1).png and results/Figure1(2).png of this repository): obtained by running "**scripts/d1/AdapMSE.m**" (main part: **Lines 23-35 and 163-264**) in the "scripts/d1" folder

*   “**Figure 2**” of the manuscript (i.e., results/Figure2(1).png and results/Figure2(2).png of this repository): obtained by running "**scripts/d1/AdapMSE2.m**" (main part: **Lines 23-35 and 164-268**) in the "scripts/d1" folder

*   “**Figure 3**” of the manuscript (i.e., results/Figure3(1).png and results/Figure3(2).png of this repository): obtained by running "**scripts/d3/Adap13D.m**" (main part: **Lines 23-35 and 73-174**) in the "scripts/d3" folder

*   “**Figure 4**” of the manuscript (i.e., results/Figure4(1).png and results/Figure4(2).png of this repository): obtained by running "**scripts/d10/Adap10D.m**" (main part: **Lines 22-35 and 164-257**) in the "scripts/d10" folder

*   “**Figures 5 and 6**” of the manuscript (i.e., results/Figure5.png and results/Figure6.png of this repository): obtained by running "**scripts/queue/manuscript/AdapQ.m**" in the "scripts/queue/manuscript" folder

#### Figures in online supplement:

*   “**Figure 3**” of the online supplement (i.e., results/OS52Case1(1).png and results/OS52Case1(2).png of this repository): obtained by running "**scripts/d1/AdapMSE.m**" (main part: **Lines 66-95 and 265-328**) in the "scripts/d1" folder

*   “**Figure 4**” of the online supplement (i.e., results/OS52Case2(1).png and results/OS52Case2(2).png of this repository): obtained by running "**scripts/d1/AdapMSE.m**" (main part: **Lines 36-65 and 98-162**) in the "scripts/d1" folder

*   “**Figure 5**” of the online supplement (i.e., results/OS52Case3(1).png and results/OS52Case3(2).png of this repository): obtained by running "**scripts/d2/Adap12D.m**" (main part: **Lines 38-67 and 73-140**) in the "scripts/d2" folder

*   “**Figure 6**” of the online supplement (i.e., results/OS52Case4(1).png and results/OS52Case4(2).png of this repository): obtained by running "**scripts/d3/Adap13D.m**" (main part: **Lines 36-66 and 175-241**) in the "scripts/d3" folder

*   “**Figure 7**” of the online supplement (i.e., results/OS52Case5(1).png and results/OS52Case5(2).png of this repository): obtained by running "**scripts/d1/AdapMSE2.m**" (main part: **Lines 66-89 and 269-338**) in the "scripts/d1" folder

*   “**Figure 8**” of the online supplement (i.e., results/OS52Case6(1).png and results/OS52Case6(2).png of this repository): obtained by running "**scripts/d1/AdapMSE2.m**" (main part: **Lines 36-65 and 93-160**) in the "scripts/d1" folder

*   “**Figure 9**” of the online supplement (i.e., results/OS52Case7(1).png and results/OS52Case7(2).png of this repository): obtained by running "**scripts/d10/Adap10D.m**" (main part: **Lines 66-90 and 261-329**) in the "scripts/d10" folder

*   “**Figure 10**” of the online supplement (i.e., results/OS52Case8(1).png and results/OS52Case8(2).png of this repository): obtained by running "**scripts/d10/Adap10D.m**" (main part: **Lines 36-65 and 93-161**) in the "scripts/d10" folder

#### Table in online supplement:

*   Table of **Online Supplement's Section 5.3**: obtained by running "**scripts/queue/online supplement/Queue_pred.m**" in the "scripts/queue/online supplement" folder



