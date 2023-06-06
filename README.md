# Nonparametric Derivative Estimation with Weighted Difference Quotients

This repository contains code for my preliminary exam, which implements the derivative estimation method proposed by the following paper.
- Paper reference: **Yu Liu and Kris De Brabanter**. [Smoothed Nonparametric Derivative Estimation using
Weighted Difference Quotients](https://www.jmlr.org/papers/volume21/19-246/19-246.pdf). _Journal of Machine Learning Research_, 21(1): 2438-2482, 2020.
- We provide both Python3 and R implementation of the proposed derivative estimation method. As for reproducing the simulation studies and their associated figures in the paper, we use our R implementation.

### Abstract for My Prelim Presentation

Precise derivative estimation is of great significance in various applications, ranging from exploring curve structures and implementing biased-corrected inference procedures in Statistics to analyzing human movements kinematically and conducting studies beyond the field of Statistics. In this presentation, we will discuss a nonparametric derivative estimation method for random design proposed in the paper by Liu and De Brabanter (2020). Specifically, we will examine their data-driven derivative estimation framework, which combines weighted difference quotients with local polynomial regression. In addition to presenting the asymptotic properties of the proposed derivative estimators in the paper, we will also fill the theoretical gaps by establishing the consistency results of the final proposed derivative estimators. Finally, we extend their simulation studies by comparing the proposed derivative estimators with other classical derivative estimation methods in terms of both estimation accuracy and computational efficiency.
