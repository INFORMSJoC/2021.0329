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
@article{CacheTest,
  author =        {C. Li, S. Gao and J. Du},
  publisher =     {INFORMS Journal on Computing},
  title =         {Convergence Analysis of Stochastic Kriging-Assisted Simulation with Random Covariates},
  year =          {2022},
  doi =           {test},
  url =           {https://github.com/INFORMSJoC/JoCTemplate},
}  
```

## Description

The goal of this software is to illustrate the the convergence behaviors of stochastic kriging-assisted simulation with random covariates.

## Replicating

*   “scripts/d1”:  testing the rates of maximal IMSE and IPFS and comparing the random design with the **Adaptive MSE procedure** in the **one**\-dimensional benchmark function examples,
*   “scripts/d2”:  testing the rates of maximal IMSE and IPFS and comparing the random design with the **Adaptive MSE procedure** in the **two**\-dimensional benchmark function examples,
*   “scripts/d3”:  testing the rates of maximal IMSE and IPFS and comparing the random design with the **Adaptive MSE procedure** in the **three**\-dimensional benchmark function examples,
*   “scripts/d10”: testing the rates of maximal IMSE and IPFS and comparing the random design with the **Adaptive MSE procedure** in the **ten**\-dimensional benchmark function examples.
*   “scripts/queue”: testing the rates of maximal IMSE and IPFS and showing how the theory can help decide the sample size to achieve a target precision in the **MM1 queue**\-dimensional benchmark function examples.

Figure 1 in the paper shows the results of the multiplication test with different
values of K using `gcc` 7.5 on an Ubuntu Linux box.

![Figure 1](results/mult-test.png)

Figure 2 in the paper shows the results of the sum test with different
values of K using `gcc` 7.5 on an Ubuntu Linux box.

![Figure 1](results/sum-test.png)


## Ongoing Development

This code is being developed on an on-going basis at the author's
[Github site](https://github.com/tkralphs/JoCTemplate).

## Support

For support in using this software, submit an
[issue](https://github.com/tkralphs/JoCTemplate/issues/new).
