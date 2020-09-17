


# DistancePreservingMatrixSketch

**Content**
This is the code repository for the research publication "A Distance-preserving Matrix Sketch" by Hengrui Luo, Giovanni Nattino and [Matthew T. Pratola](http://www.matthewpratola.com/). 
The manuscript of this paper can be accessed at https://arxiv.org/abs/1908.08864. 

 - In [Python folder](https://github.com/hrluo/SparseAdditiveGaussianProcess/tree/master/Python), we provided a set of illustrative code that serves as a proof of concept, and also a set of robust code that can be executed for large datasets.
 - In [R folder](https://github.com/hrluo/SparseAdditiveGaussianProcess/tree/master/R), we provided user-friendly R code with more features like adaptive Metropolis-Hasting step-width. However, we point out that its computational efficiency is not comparable to Python version.

**Abstract**
In this paper we introduce a novel model for Gaussian process (GP) regression in the fully Bayesian setting. Motivated by the ideas of sparsification, localization and Bayesian additive modeling, our model is built around a recursive partitioning (RP) scheme. Within each RP partition, a sparse GP (SGP) regression model is fitted. A Bayesian additive framework then combines multiple layers of partitioned SGPs, capturing  both global trends and local refinements with efficient computations. The model addresses both the problem of efficiency in fitting a full Gaussian process regression model and the problem of prediction performance associated with a single SGP. Our approach mitigates the issue of pseudo-input selection and avoids the need for complex inter-block correlations in existing methods.  The crucial trade-off becomes choosing between many simpler local model
components or fewer complex global model components, which the practitioner can sensibly tune. Implementation is via a Metropolis-Hasting Markov chain Monte-Carlo algorithm with Bayesian back-fitting.
We compare our model against popular alternatives on simulated and real datasets, and find the performance is competitive, while the fully Bayesian procedure enables the quantification of model uncertainties. 

**Citation**
We provided both iPynb illustrative code, Python/R production code for reproducible and experimental purposes under [LICENSE](https://github.com/hrluo/SparseAdditiveGaussianProcess/blob/master/LICENSE).
Please cite our paper using following BibTeX item:

    @article{luo2019sparse,
        title={Sparse Additive Gaussian Process Regression},
        author={Hengrui Luo and Giovanni Nattino and Matthew T. Pratola},
        year={2019},
        eprint={1908.08864},
        archivePrefix={arXiv},
        primaryClass={math.ST}
    }

Thank you again for the interest and please reach out if you have further questions.
