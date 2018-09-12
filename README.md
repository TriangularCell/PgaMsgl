# PgaMsgl

PgaMsgl is an *R* package on Proximal Gradient Algorithm for Multi-variate Sparse Group Lasso. It is mainly implemented with *Rcpp* and *RcppEigen*.

### Installing

Download the whole package and **Build** within R.
Later on I will release it on CRAN.

### Running the tests

The package includes two demo datasets for testing usage. One is with relative low dimension, and another is with relative high dimension.

For example, to run the low dimension dataset:

```
data(lowD)
system.time(lowD_result <- PgaMsgl(lowD$X, lowD$Y, lowD$B0, model="L121", lowD$Gm, lowD$mi, lowD$mg, lowD$mc))
```

### Authors

**Yiming Qin** - *Initial work* - [TriangularCell](https://github.com/TriangularCell)

### License

This project is licensed under the GPL-3.0 License - see the [LICENSE.md](LICENSE.md) file for details

### Acknowledgments

* Yaohua Hu
* Jing Qin
* Xinlin Hu
* Sifan Liu

### References

Hu, Yaohua, et al. "[Group sparse optimization via lp, q regularization.](http://www.jmlr.org/papers/volume18/15-651/15-651.pdf)" J Mach Learn Res 18 (2017): 1-52.
