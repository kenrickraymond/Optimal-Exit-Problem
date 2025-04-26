# A Numerical Scheme for the Optimal Liquidation Problem Under Jump Diffusion Dynamics on High-Frequency Data

This repository contains resources relevant to the project **"A Numerical Scheme for the Optimal Liquidation Problem Under Jump Diffusion Dynamics on High-Frequency Data."**

## Project Overview
The project focuses on solving the optimal stopping time problem for an investor who aims to liquidate an asset based on a stop-loss level. The value function for this problem is modeled as a free-boundary problem, solved in a manner similar to pricing an American option. Given the dynamics of the stock price, which includes jumps, no closed-form solutions existed at the time. Therefore, I used finite differences and quadrature methods to obtain an estimate.

## Resources

- [Research PDF](https://github.com/kenrickraymond/Optimal-Exit-Problem/blob/main/A_Numerical_Scheme_for_the_Optimal_Liquidation_Problem_Under_Jump_Diffusion_Dynamics_on_High_Frequency_Data.pdf): Includes a literature review, methodology, and a numerical example.

- [Python Code to Obtain Cointegrated Pair](https://github.com/kenrickraymond/Optimal-Exit-Problem/blob/main/dataRetrieverLatest.ipynb): Python code used for retrieving and analyzing cointegrated pairs.

- [R Code (Jupyter Notebook) for Methodology Implementation](https://github.com/kenrickraymond/Optimal-Exit-Problem/blob/main/Optimal%20Exit%20Problem.ipynb): R code implemented in a Jupyter notebook to execute the methodology for solving the optimal stopping problem.
