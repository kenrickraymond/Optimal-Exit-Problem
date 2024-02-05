This repository contains resources relevant to the project "A Numerical Scheme for the Optimal Liquidation Problem Under Jump Diffusion Dynamics on High-Frequency Data." Essentially, an investor is interested in the optimal stopping time based on a stop-loss level. The value function can be represented as a free-boundary problem, which was solved analogously to how an American option is priced. Due to the dynamics of the stock price, which includes a jump, no closed-form solutions existed at the time hence I resorted to finite differences and quadrature to obtain an estimate.

* [A research PDF](https://github.com/kenrickraymond/Optimal-Exit-Problem/blob/main/A_Numerical_Scheme_for_the_Optimal_Liquidation_Problem_Under_Jump_Diffusion_Dynamics_on_High_Frequency_Data.pdf) including a literature review, methodology, and numerical example.

* [Code written in Python to obtain cointegrated pair](https://github.com/kenrickraymond/Optimal-Exit-Problem/blob/main/dataRetrieverLatest.ipynb)

* [Code written in R (Jupyter Notebook) to implement methodology](https://github.com/kenrickraymond/Optimal-Exit-Problem/blob/main/Optimal%20Exit%20Problem.ipynb)
