# WeightedPLS
This is a modification of the Rotman Baycrest PLS software to allow for the use of sampling weights (e.g. non-participation weights or post-stratification weights) with the PLS algorithm. Original software here https://www.rotman-baycrest.on.ca/index.php?section=84

The modifications are based on weighted versions of the mean, variance, and covariance outlined in this paper https://www.sciencedirect.com/science/article/abs/pii/S0263237316300688

The primary code changes are in rri_xcor.m
