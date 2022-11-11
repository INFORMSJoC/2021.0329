[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# 2021-0329

This archive is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

The scripts in this repository are a snapshot of the scripts
that were used in the research reported on in the paper 
Convergence Analysis of Stochastic Kriging-Assisted Simulation with Random Covariates by C. Li, S. Gao and J. Du. 

## Cite

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

## Scripts

*   **"scripts/d1"** folder:  testing the rates of maximal IMSE and IPFS under the random sampling and comparing the random design with the Adaptive MSE procedure in the **one**\-dimensional benchmark function examples,
*   **"scripts/d2"** folder:  testing the rates of maximal IMSE and IPFS under the random sampling and comparing the random design with the Adaptive MSE procedure in the **two**\-dimensional benchmark function examples,
*   **"scripts/d3"** folder:  testing the rates of maximal IMSE and IPFS under the random sampling and comparing the random design with the Adaptive MSE procedure in the **three**\-dimensional benchmark function examples,
*   **"scripts/d10"** folder: testing the rates of maximal IMSE and IPFS under the random sampling and comparing the random design with the Adaptive MSE procedure in the **ten**\-dimensional benchmark function examples.
*   **"scripts/queue"** folder: testing the rates of maximal IMSE and IPFS under the random sampling and showing how the theory can help decide the sample size to achieve a target precision in the **M/M/1 queue** example.

## Results

The **README.md** in the "results" subdirectory tells the detailed script to obtain each result of the numerical experiment. The **online supplement** of this paper is also in this subdirectory.

