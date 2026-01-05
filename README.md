[![arXiv][arxiv-shield]][arxiv-url]
[![DOI][doi-shield]][doi-url]

# [SubConD-Small-EVs-SVs-App][arxiv-url]

We consider the approximation of the smallest eigenvalue of a large parameter-dependent Hermitian matrix, with parameters lying in a continuum compact domain. The approach is based on approximating the smallest
eigenvalue by the one obtained by projecting the large matrix onto a suitable (small) subspace. The singular value case, which is of interest when the matrix is non-Hermitian, is also considered with an original approach. This folder contains the Matlab codes to reproduce the numerical experiments shown in the paper [Uniform Approximation of Eigenproblems of a Large-Scale Parameter-Dependent Hermitian Matrix][arxiv-url].

##Code info:

* Each script, following their denomination, reproduces the corresponding Figure or Table data in [preprint][arxiv-url]; 

Algorithm 2 Func.:

* **approx\_smallesteig\_all**: corresponds to Algorithm 2 in [preprint][arxiv-url];
* **approx\_smallesteig\_all\_hybrid**: corresponds to first using Algorithm 2 and then the [**SSCM**][Ref1], see [Section 6.2.1][arxiv-url];
* **lamin\_error\_all**: given a parameter $\mu$ the function evaluates $H^{(j)}(\mu)$ (or $H_{r}^{(j)}(\mu)$), see  [preprint][arxiv-url]; 
* **subspace\_SCMM**: script to run the [**SSCM**][Ref1] method;


Algorithm 3 Func.:

* **SSSV\_EIGOPT**: stands for subspace smallest singular value (SSSV) with [EIGOPT][Ref2] and corresponds to Algorithm 3 in [preprint][arxiv-url];
* **SSSV**: corresponds to use Algorithm 3 with a discrete grid;
* **sigma\_error\_all**: given a parameter $\mu$ the function evaluates $S^{(j)}(\mu)$, see [preprint][arxiv-url]; 

eigopt: the folder containing the functions that perform the optimization over the continuum domain $\mathcal{D}$, the software is developed by the authors of [EigOpt][Ref2];

Plot\_Functions: a folder that contains functions to visualise the numerical results. In particular: **func_plots** and **LB_ULB\_UB** are used to produce Figure 2 of [preprint][arxiv-url];


## Citing
If you use this project for academic work, please consider citing our
[publication][arxiv-url]:

    N. Guglielmi, M. Manucci and E. Mengi
    Uniform Approximation of Eigenproblems of a Large-Scale Parameter-Dependent Hermitian Matrix
    ArXiv e-print 2409.05791, 2024.
    
## License
Distributed under the MIT License. See `LICENSE` for more information.

## Other references

* [Subspace Acceleration for Large-Scale Parameter-Dependent Hermitian Eigenproblems][Ref1]
* [Numerical Optimization of Eigenvalues of Hermitian Matrix Functions][Ref2]

## Contacts

* Mattia Manucci - [mattia.manucci@simtech.uni-stuttgart.de](mattia.manucci@simtech.uni-stuttgart.de)
* Emre Mengi - [emengi@ku.edu.tr](emengi@ku.edu.tr)
* Nicola Guglielmi - [nicola.guglielmi@gssi.it](nicola.guglielmi@gssi.it)



[doi-shield]: https://img.shields.io/badge/DOI-10.5281%20%2F%20zenodo.13254480-blue.svg?style=for-the-badge
[doi-url]: https://doi.org/10.5281/zenodo.13254480
[link-drive]: https://drive.google.com/file/d/1Y-RDkTQOvLaeccjgYzoq_qsD4tQLqN5u/view?usp=sharing
[arxiv-shield]: https://img.shields.io/badge/arXiv-2204.13474-b31b1b.svg?style=for-the-badge
[arxiv-url]: https://arxiv.org/abs/2409.05791

[Ref1]: https://doi.org/10.1137/15M1017181
[Ref2]: https://doi.org/10.1137/130933472
