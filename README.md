[![DOI](https://zenodo.org/badge/681533852.svg)](https://zenodo.org/badge/latestdoi/681533852)
# Fortran Fast math
A collection of functions for fast number crunching using Fortran.

In order to get the maximum performance of this library, compile with "-O3 -march=native".

# Available functions

* fast (and precise) sum for 1D arrays - possibility of including a mask.
    fsum: fastest method and at worst, 1 order of magnitud more precise than the intrinsic sum. It groups chunks of values in a temporal working batch which is summed up once at the end.
    fsum_pair: Highest precision. It has a precission equivalent to a quadratic sum (for real32 summing with real64, and fo real64 summing with real128). runtime can be slightly slower or just as fast as the intrinsic sum.

* fast (and precise) dot product for 1D arrays - possibility of including a mask.
    fprod: fastest method and at worst, 1 order of magnitud more precise than the intrinsic dot_product. It groups chunks of products in a temporal working batch which is summed up once at the end.

## Trigonometric
* fast cosinus
* fast sinus
* fast tangent
* fast acosinus
* fast tanh

## Other
* fast erf

# API documentation

To generate the API documentation for `fast_math` using
[ford](https://github.com/Fortran-FOSS-Programmers/ford) run the following
command:

```shell
ford ford.yml
```

# TODO
* Contribution guidelines
* Finalize naming convention
* Integrate github CICD
* Polish autodoc

# Acknowledgement

* Compilation of this library was possible thanks to [Transvalor S.A.](https://www.transvalor.com/en/homepage) research activities. 
* Part of this library is based on the work of [Perini and Reitz](https://doi.org/10.1016/j.combustflame.2018.04.013), that was funded through the Sandia National Laboratories by the U.S. Department of Energy, Office of Vehicle Technologies, program managers Leo Breton, Gupreet Singh.

Contribution of open-source developers:

[jalvesz](https://github.com/jalvesz)

[perazz](https://github.com/perazz)