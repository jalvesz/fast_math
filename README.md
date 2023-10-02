[![DOI](https://zenodo.org/badge/681533852.svg)](https://zenodo.org/badge/latestdoi/681533852)
# Fortran Fast math
A collection of functions for fast number crunching using Fortran.

In order to get the maximum performance of this library, compile with "-O3 -march=native" (or equivalent).

# Available functions

| function | name(s)               | shapes     | types            | 
|----------|-----------------------|------------|------------------|
| sum      | `fsum` `fsum_pair`(1) |        `1d`|`real32` `real64` |
| dot      | `fprod`(2)            |        `1d`|`real32` `real64` |
| cos      | `fcos`                | `elemental`|`real32` `real64` |
| sin      | `fsin`                | `elemental`|`real32` `real64` |
| tan      | `ftan`                | `elemental`|`real32` `real64` |
| tanh     | `ftanh`               | `elemental`|`real32` `real64` |
| acos     | `facos`               | `elemental`|`real32` `real64` |
| erf      | `ferf`                | `elemental`|`real32` `real64` |
| log      | `flog_p3` `flog_p5`   | `elemental`|         `real64` |
| rsqrt(3) | `frsqrt`              | `elemental`|`real32` `real64` |

* (1) fast (and precise) sum for 1D arrays - possibility of including a mask.
    `fsum`: fastest method and at worst, same or 1 order of magnitud more precise than the intrinsic sum. runtime can vary between 3X and 9X the intrinsic. It groups chunks of values in a temporal working batch which is summed up once at the end.
    `fsum_pair`: Highest precision. It has a precission equivalent to a quadratic sum (for real32 summing with real64, and fo real64 summing with real128). runtime can vary between 0.9X to 1.6X the intrinsic sum.

* (2) fast (and precise) dot product for 1D arrays - possibility of including a 3rd weighting array.
    `fprod`: fastest method and at worst, 1 order of magnitud more precise than the intrinsic dot_product. runtime can vary between 3X and 8X the intrinsic. It groups chunks of products in a temporal working batch which is summed up once at the end (based on `fsum`).
* (3) rsqrt: reciprocal square root $f(x)=1/sqrt(x)$
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
* The [fortran lang community](https://fortran-lang.discourse.group/) discussions such as [Some Intrinsic SUMS](https://fortran-lang.discourse.group/t/some-intrinsic-sums/5760) and [fastGPT](https://fortran-lang.discourse.group/t/fastgpt-faster-than-pytorch-in-300-lines-of-fortran/5385)

Contribution of open-source developers:

[jalvesz](https://github.com/jalvesz)

[perazz](https://github.com/perazz)