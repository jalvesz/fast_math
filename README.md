[![DOI](https://zenodo.org/badge/681533852.svg)](https://zenodo.org/badge/latestdoi/681533852)
# Fortran Fast math
A collection of functions for fast number crunching using Fortran.

In order to get the maximum performance of this library, compile with "-O3 -march=native" (or equivalent).

# Available functions

| function | name(s)               | shapes     | types            | 
|----------|-----------------------|------------|------------------|
| sum      | `fsum` `fsum_kahan`(1) |        `1d`|`real32` `real64` |
| dot      | `fprod` `fprod_kahan`(2)|        `1d`|`real32` `real64` |
| cos      | `fcos`                | `elemental`|`real32` `real64` |
| sin      | `fsin`                | `elemental`|`real32` `real64` |
| tan      | `ftan`                | `elemental`|`real32` `real64` |
| tanh     | `ftanh`               | `elemental`|`real32` `real64` |
| acos     | `facos`               | `elemental`|`real32` `real64` |
| atan     | `fatan`               | `elemental`|`real32` `real64` |
| erf      | `ferf`                | `elemental`|`real32` `real64` |
| log      | `flog_p3` `flog_p5`   | `elemental`|         `real64` |
| rsqrt(3) | `frsqrt`              | `elemental`|`real32` `real64` |

* (1) fast (and precise) sum for 1D arrays - possibility of including a mask.
    `fsum`: fastest method and at worst, same or 1 order of magnitud more precise than the intrinsic sum. It groups chunks of values in a temporal working batch which is summed up once at the end.
    `fsum_kahan`: Highest precision. It has a precission close to a quadratic sum (for real32 summing with real64, and fo real64 summing with real128). It also uses the chunks principle with an elemental kahan operator applied on top.

* (2) fast (and precise) dot product for 1D arrays - possibility of including a 3rd weighting array.
    `fprod`: fastest method and at worst, 1 order of magnitud more precise than the intrinsic dot_product. runtime can vary between 3X and 8X the intrinsic. It groups chunks of products in a temporal working batch which is summed up once at the end (based on `fsum`).
    `fprod_kahan`: Same idea as `fsum_kahan` but on top of chunked products.
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

# Elapsed time examples and precision
Warning: The following values are just references as to see how different can they be between different compilers. Actual speed-ups(downs) should be measured under the true use conditions to account for (lack-off) inlinement, etc etc. Results obtained using a Intel(R) Xeon(R) Gold 6226R CPU @ 2.90GHz   2.89 GHz.
<details>
<summary>(Click to unfold) WSL2 gfortran 13.2 > fpm test --flag "-cpp -O3 -march=native -flto"</summary>

|      sum r32 | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    intrinsic |           1.0300 |     1.00 |      3.1511E-06 |
|        kahan |           0.1200 |     8.58 |      9.5367E-08 |
|        chunk |           0.0900 |    11.44 |      1.0824E-07 |
 
|      sum r64 | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    intrinsic |           1.2100 |     1.00 |      5.6974E-15 |
|        kahan |           0.4300 |     2.81 |      1.3278E-16 |
|        chunk |           0.1100 |    11.00 |      2.3359E-16 |
 
| sum r32 mask | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    intrinsic |           1.2300 |     1.00 |      1.6180E-06 |
|        kahan |           4.3400 |     0.28 |      8.3327E-08 |
|        chunk |           0.3800 |     3.24 |      8.8394E-08 |
 
| sum r64 mask | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    intrinsic |           3.8600 |     1.00 |      2.9463E-15 |
|        kahan |           4.1950 |     0.92 |      6.8723E-17 |
|        chunk |           0.4200 |     9.19 |      1.1879E-16 |
 
|      dot r32 | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    intrinsic |           1.2100 |     1.00 |      3.2994E-06 |
|        kahan |           0.2300 |     5.26 |      9.8348E-08 |
|        chunk |           0.1200 |    10.08 |      1.1307E-07 |
 
|      dot r64 | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    intrinsic |           1.1900 |     1.00 |      5.9648E-15 |
|        kahan |           0.4400 |     2.70 |      1.2812E-16 |
|        chunk |           0.0900 |    13.22 |      2.2760E-16 |

|        trigo | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|     fsin r32 |           0.5280 |     7.54 |      3.0972E-07 | 
|     fsin r64 |           0.9320 |     8.76 |      3.9779E-16 | 
|    facos r32 |           0.3080 |    20.87 |      2.9135E-05 | 
|    facos r64 |           0.5960 |    15.90 |      2.1557E-14 | 

|       hyperb | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    ftanh r32 |           1.7280 |    10.07 |      7.4200E-08 | 
|    ftanh r64 |           1.9360 |     9.32 |      1.3282E-09 | 
|     ferf r32 |           0.4760 |    31.71 |      9.6432E-08 | 
|     ferf r64 |           0.7640 |    18.42 |      9.6298E-08 | 

|        rsqrt | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|   frsqrt r32 |           0.3480 |     1.06 |      9.4399E-04 | 
|   frsqrt r64 |           0.6320 |     2.23 |      8.6268E-04 | 
</details>

<details>
<summary>(Click to unfold) WSL2 nvfortran 23.9 > fpm test --flag "-Mpreprocess -fast -Minline"</summary>

|      sum r32 | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    intrinsic |           0.1000 |     1.00 |      1.1295E-07 |
|        kahan |           1.2500 |     0.08 |      9.8169E-08 |
|        chunk |           0.0700 |     1.43 |      7.0930E-08 |
 
|      sum r64 | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    intrinsic |           0.1400 |     1.00 |      3.8969E-16 |
|        kahan |           1.6300 |     0.09 |      1.2623E-16 |
|        chunk |           0.2500 |     0.56 |      1.8996E-16 |
 
| sum r32 mask | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    intrinsic |           0.1700 |     1.00 |      2.0742E-07 |
|        kahan |           5.5650 |     0.03 |      8.1956E-08 |
|        chunk |           0.2550 |     0.67 |      5.8889E-08 |
 
| sum r64 mask | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    intrinsic |           0.3600 |     1.00 |      3.8136E-16 |
|        kahan |           5.7750 |     0.06 |      6.2839E-17 |
|        chunk |           0.4400 |     0.82 |      8.5598E-17 |
 
|      dot r32 | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    intrinsic |           0.1400 |     1.00 |      1.1426E-07 |
|        kahan |           1.9700 |     0.07 |      9.7811E-08 |
|        chunk |           0.1700 |     0.82 |      7.1764E-08 |
 
|      dot r64 | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    intrinsic |           0.2400 |     1.00 |      3.9246E-16 |
|        kahan |           1.8700 |     0.13 |      1.3178E-16 |
|        chunk |           0.4100 |     0.59 |      1.9129E-16 |

|        trigo | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|     fsin r32 |           0.0160 |   726.00 |      1.0325E-07 | 
|     fsin r64 |           0.0280 |   388.86 |      5.0118E-17 | 
|    facos r32 |           0.0120 |   466.67 |      1.0563E-06 | 
|    facos r64 |           0.0200 |   390.60 |      3.7996E-15 | 

|       hyperb | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    ftanh r32 |           0.0240 |   676.67 |      5.3264E-08 | 
|    ftanh r64 |           0.0080 |  1595.00 |      1.3282E-09 | 
|     ferf r32 |           0.0040 |  4851.00 |      9.1205E-08 | 
|     ferf r64 |           0.0320 |   549.62 |      9.6298E-08 | 

|        rsqrt | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|   frsqrt r32 |          16.5480 |     0.02 |      9.4387E-04 | 
|   frsqrt r64 |          15.8280 |     0.09 |      8.6745E-04 | 

</details>

<details>
<summary>(Click to unfold) WSL2 ifort 2021.10.0 > fpm test --flag "-fpp -O3 -xHost -ipo"</summary>

|      sum r32 | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    intrinsic |           0.0700 |     1.00 |      6.2262E-08 |
|        kahan |           0.2400 |     0.29 |      9.4564E-08 |
|        chunk |           0.1000 |     0.70 |      7.0930E-08 |
 
|      sum r64 | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    intrinsic |           0.0800 |     1.00 |      1.9862E-16 |
|        kahan |           0.5200 |     0.15 |      1.2867E-16 |
|        chunk |           0.1400 |     0.57 |      2.0384E-16 |
 
| sum r32 mask | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    intrinsic |           0.2000 |     1.00 |      2.0568E-07 |
|        kahan |           0.2150 |     0.93 |      7.7122E-08 |
|        chunk |           0.1450 |     1.38 |      6.7770E-08 |
 
| sum r64 mask | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    intrinsic |           0.2150 |     1.00 |      1.9040E-16 |
|        kahan |           0.4400 |     0.49 |      7.0610E-17 |
|        chunk |           0.3700 |     0.58 |      8.5154E-17 |
 
|      dot r32 | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    intrinsic |           0.0700 |     1.00 |      6.2031E-08 |
|        kahan |           0.2100 |     0.33 |      1.0544E-07 |
|        chunk |           0.0500 |     1.40 |      7.1526E-08 |
 
|      dot r64 | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    intrinsic |           0.2200 |     1.00 |      6.3782E-16 |
|        kahan |           0.4600 |     0.48 |      2.4047E-16 |
|        chunk |           0.1200 |     1.83 |      1.8829E-16 |

|        trigo | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|     fsin r32 |           0.3560 |     1.26 |      1.9746E-07 | 
|     fsin r64 |           0.9280 |     1.38 |      7.5661E-17 | 
|    facos r32 |           0.3200 |     2.01 |      3.0743E-06 | 
|    facos r64 |           0.6520 |     3.36 |      6.3642E-15 | 

|       hyperb | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    ftanh r32 |           0.3960 |     3.70 |      1.1537E-08 | 
|    ftanh r64 |           0.6760 |     5.17 |      1.3282E-09 | 
|     ferf r32 |           0.3360 |     2.50 |      1.0924E-07 | 
|     ferf r64 |           0.8440 |     2.18 |      9.6298E-08 | 

|        rsqrt | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|   frsqrt r32 |           0.2600 |     1.31 |      9.4032E-04 | 
|   frsqrt r64 |           0.6360 |     2.27 |      8.7360E-04 | 
</details>

<details>
<summary>(Click to unfold) Windows ifx 2023.2.0 > fpm test --flag "-fpp -O3 -xHost -ipo"</summary>

|      sum r32 | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    intrinsic |           0.3200 |     1.00 |      8.4376E-07 |
|        kahan |           1.0300 |     0.31 |      8.7321E-08 |
|        chunk |           0.4800 |     0.67 |      8.7082E-08 |
 
|      sum r64 | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    intrinsic |           1.1200 |     1.00 |      5.7371E-15 |
|        kahan |           0.9400 |     1.19 |      1.9507E-16 |
|        chunk |           0.5600 |     2.00 |      1.9418E-16 |
 
| sum r32 mask | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    intrinsic |           2.2700 |     1.00 |      1.5584E-06 |
|        kahan |           4.4750 |     0.51 |      9.1434E-08 |
|        chunk |           4.7200 |     0.48 |      8.7559E-08 |
 
| sum r64 mask | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    intrinsic |           2.1750 |     1.00 |      2.9075E-15 |
|        kahan |           4.7550 |     0.46 |      1.0636E-16 |
|        chunk |           4.0250 |     0.54 |      1.0525E-16 |
 
|      dot r32 | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    intrinsic |           0.2600 |     1.00 |      7.9530E-07 |
|        kahan |           1.3800 |     0.19 |      6.8307E-08 |
|        chunk |           0.4900 |     0.53 |      6.9737E-08 |
 
|      dot r64 | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    intrinsic |           0.6200 |     1.00 |      2.9848E-15 |
|        kahan |           1.4200 |     0.44 |      1.8197E-16 |
|        chunk |           0.5800 |     1.07 |      1.8330E-16 |

|        trigo | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|     fsin r32 |           3.4640 |     0.47 |      1.3924E-07 | 
|     fsin r64 |           3.2320 |     1.31 |      1.0296E-15 | 
|    facos r32 |           1.3960 |     5.22 |      3.1710E-05 | 
|    facos r64 |           1.4080 |     6.28 |      5.2928E-13 | 

|       hyperb | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    ftanh r32 |           2.8280 |     1.22 |      2.3012E-08 | 
|    ftanh r64 |           2.6280 |     2.97 |      1.3282E-09 | 
|     ferf r32 |           3.8600 |     1.57 |      3.0995E-07 | 
|     ferf r64 |           3.9600 |     5.67 |      9.6298E-08 | 

|        rsqrt | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|   frsqrt r32 |           1.6640 |     0.19 |      9.4038E-04 | 
|   frsqrt r64 |           1.4320 |     0.96 |      8.7360E-04 |
</details>

# Acknowledgement

* Compilation of this library was possible thanks to [Transvalor S.A.](https://www.transvalor.com/en/homepage) research activities. 
* Part of this library is based on the work of [Perini and Reitz](https://doi.org/10.1016/j.combustflame.2018.04.013), that was funded through the Sandia National Laboratories by the U.S. Department of Energy, Office of Vehicle Technologies, program managers Leo Breton, Gupreet Singh.
* The [fortran lang community](https://fortran-lang.discourse.group/) discussions such as [Some Intrinsic SUMS](https://fortran-lang.discourse.group/t/some-intrinsic-sums/5760) and [fastGPT](https://fortran-lang.discourse.group/t/fastgpt-faster-than-pytorch-in-300-lines-of-fortran/5385)

Contribution of open-source developers:

[jalvesz](https://github.com/jalvesz)

[perazz](https://github.com/perazz)
