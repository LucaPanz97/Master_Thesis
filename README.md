# Master_Thesis

# Goal

Robust clustering for high-dimensional data represents a significant challenge, as existing robust clustering methods suffer from the curse of dimensionality when the number of variables is large, while existing approaches for high-dimensional data are, in general, not robust. We propose a solution to this challenge, by incorporating high-dimensional covariance matrix estimators into a fast and efficient algorithm for robust constrained clustering, the TCLUST methodology, which has been extensively shown to perform well on contaminated low-dimensional data. We employ three different covariance matrix estimators to enhance TCLUST: the Minimum Regularized Covariance Determinant estimator, the linear shrinkage estimator of Ledoit-Wolf (which is a special case of the previous one) and the sparse CovGlasso estimator. This enhancement enables the algorithm to effectively handle the complexity and variability of high-dimensional data while protecting against the harmful effects of outliers.

# Real-World Challenge: Handwritten Digits Recognition  

After testing and comparing our methodologies on high-dimensional Gaussian simulated data, we address a real-world challenge: handwritten digits recognition. We focus on two distinct problems within handwritten digits recognition: clustering digits 0, 1 and 4 with contamination, where we tackle the challenges of data contamination and high dimensionality. Second, we deal with clustering digits 3, 5 and 8 with contamination, which represents the most complex scenario in handwritten digits recognition due the high degree of similarity between these digits. In addition to the previously mentioned challenges, we also address the issue of limited class separation.

# Conclusions

CovGlasso in TCLUST is proved to be our most robust and effective solution for addressing handwritten digits recognition problems.

# Future Steps



# Authors

Luca Panzeri    *Politecnico di Milano - MSc in Statistical Learning, Mathematical Engineering*

Davide Zaltieri    *Politecnico di Milano - MSc in Statistical Learning, Mathematical Engineering*

# Acknowledgments

Special thanks to Prof. Cappozzo for his guidance and support during this work.
