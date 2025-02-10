[![DOI](https://zenodo.org/badge/681533852.svg)](https://zenodo.org/badge/latestdoi/681533852)
# Fortran Fast math
A collection of functions for fast number crunching using Fortran.

In order to get the maximum performance of this library, compile with "-O3 -march=native -flto" (or equivalent). Note: For the elemental functions, inlinement is key to extract maximum performance. It can be achieved either by use of the `-flto`(gcc)/`-ipo`(intel) flag or using the `include` mechanism.

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
* Polish autodoc

# Elapsed time examples and precision
Warning: The following values are just references as to see how different can they be between different compilers. Actual speed-ups(downs) should be measured under the true use conditions to account for (lack-off) inlinement, etc etc.
<details>
<summary>(Click to unfold) Windows gfortran 14.1 > fpm test --flag "-O3 -march=native -mtune=native"</summary>
CPU: Intel(R) Core(TM) i7-8565U CPU @ 1.80GHz   1.99 GHz

|      sum r32 | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    intrinsic |           1.2100 |     1.00 |      3.3794E-06 |
|        kahan |           0.1800 |     6.72 |      1.0425E-07 |
|        chunk |           0.1100 |    11.00 |      1.1265E-07 |
 
|      sum r64 | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    intrinsic |           1.3000 |     1.00 |      5.9269E-15 |
|        kahan |           0.3100 |     4.19 |      1.7286E-16 |
|        chunk |           0.1500 |     8.67 |      2.1416E-16 |
 
| sum r32 mask | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    intrinsic |           4.1250 |     1.00 |      1.5687E-06 |
|        kahan |           0.1600 |    25.78 |      9.1493E-08 |
|        chunk |           0.1600 |    25.78 |      8.8453E-08 |
 
| sum r64 mask | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    intrinsic |           4.0350 |     1.00 |      2.9428E-15 |
|        kahan |           0.3750 |    10.76 |      1.2179E-16 |
|        chunk |           0.2450 |    16.47 |      1.2768E-16 |
 
|      dot r32 | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    intrinsic |           1.0600 |     1.00 |      3.2735E-06 |
|        kahan |           0.1500 |     7.07 |      9.8348E-08 |
|        chunk |           0.1000 |    10.60 |      1.1587E-07 |

|      dot r64 | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    intrinsic |           1.2100 |     1.00 |      5.8091E-15 |
|        kahan |           0.3300 |     3.67 |      1.8407E-16 |
|        chunk |           0.2000 |     6.05 |      2.0528E-16 |

|        trigo | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|     fsin r32 |           2.8840 |    13.82 |      3.4749E-07 |
|     fsin r64 |           3.1040 |    12.17 |      4.0784E-16 |
|    facos r32 |           1.6600 |    28.64 |      2.9135E-05 | 
|    facos r64 |           1.6800 |     6.89 |      2.9274E-14 | 
|    fatan r32 |           1.6720 |    23.36 |      1.7730E-06 | 
|    fatan r64 |           2.5120 |     3.94 |      6.6869E-06 | 

|       hyperb | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    ftanh r32 |           2.1640 |     8.61 |      5.9480E-08 | 
|    ftanh r64 |           2.3480 |     7.16 |      1.3282E-09 | 
|     ferf r32 |           2.3600 |    27.21 |      7.9573E-08 | 
|     ferf r64 |           4.1200 |    15.60 |      9.6298E-08 | 

|        rsqrt | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|   frsqrt r32 |           1.7720 |     0.26 |      9.4039E-04 | 
|   frsqrt r64 |           2.2280 |     0.64 |      8.9297E-04 | 
</details>

<details>
<summary>(Click to unfold) Windows ifx 2025.0.4 > fpm test --flag "/O3 /Qxhost"</summary>
CPU: Intel(R) Core(TM) i7-8565U CPU @ 1.80GHz   1.99 GHz

|      sum r32 | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    intrinsic |           0.4300 |     1.00 |      3.8308E-07 |
|        kahan |           0.1700 |     2.53 |      6.0938E-08 |
|        chunk |           0.0100 |    43.00 |      6.0938E-08 |
 
|      sum r64 | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    intrinsic |           0.3500 |     1.00 |      1.5061E-15 |
|        kahan |           0.1800 |     1.94 |      1.3033E-16 |
|        chunk |           0.0200 |    17.50 |      1.3886E-16 |
 
| sum r32 mask | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    intrinsic |           0.3000 |     1.00 |      2.0369E-07 |
|        kahan |           0.2200 |     1.36 |      5.2360E-08 |
|        chunk |           0.1750 |     1.71 |      5.2515E-08 |
 
| sum r64 mask | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    intrinsic |           0.3500 |     1.00 |      3.7423E-16 |
|        kahan |           0.2900 |     1.21 |      8.3862E-17 |
|        chunk |           0.2800 |     1.25 |      9.4422E-17 |
 
|      dot r32 | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    intrinsic |           0.3400 |     1.00 |      3.9539E-07 |
|        kahan |           0.1600 |     2.12 |      6.7639E-08 |
|        chunk |           0.1600 |     2.12 |      6.6906E-08 |
 
|      dot r64 | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    intrinsic |           0.7100 |     1.00 |      1.4730E-15 |
|        kahan |           0.1500 |     4.73 |      1.2270E-16 |
|        chunk |           0.1700 |     4.18 |      1.2459E-16 |

|        trigo | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|     fsin r32 |           3.0960 |     0.26 |      2.0412E-08 | 
|     fsin r64 |           2.7080 |     1.01 |      3.5190E-17 | 
|    facos r32 |           1.6440 |     0.46 |      1.3946E-05 | 
|    facos r64 |           1.7560 |     1.51 |      2.0708E-11 | 
|    fatan r32 |           2.6880 |     0.28 |      4.4950E-06 | 
|    fatan r64 |           1.9000 |     1.73 |      6.6869E-06 | 

|       hyperb | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    ftanh r32 |           2.3200 |     0.48 |      1.0284E-08 | 
|    ftanh r64 |           2.3080 |     2.19 |      1.3282E-09 | 
|     ferf r32 |           3.3160 |     0.23 |      7.5974E-07 | 
|     ferf r64 |           2.9760 |     0.89 |      9.6298E-08 | 

|        rsqrt | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|   frsqrt r32 |           1.7280 |     0.21 |      9.4033E-04 | 
|   frsqrt r64 |           1.6520 |     0.90 |      8.7360E-04 |
</details>

<details>
<summary>(Click to unfold) WSL2 nvfortran 24.3 > fpm test --flag "-Mpreprocess -fast"</summary>
CPU: Intel(R) Core(TM) i7-8565U CPU @ 1.80GHz   1.99 GHz

|      sum r32 | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    intrinsic |           0.2100 |     1.00 |      1.1295E-07 |
|        kahan |           0.3200 |     0.66 |      9.8169E-08 |
|        chunk |           0.1400 |     1.50 |      7.1764E-08 |
 
|      sum r64 | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    intrinsic |           0.3300 |     1.00 |      3.8969E-16 |
|        kahan |           0.3200 |     1.03 |      1.8086E-16 |
|        chunk |           0.2200 |     1.50 |      9.0372E-17 |
 
| sum r32 mask | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    intrinsic |           0.2400 |     1.00 |      2.0742E-07 |
|        kahan |           0.3050 |     0.79 |      8.9645E-08 |
|        chunk |           0.1550 |     1.55 |      5.8651E-08 |
 
| sum r64 mask | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    intrinsic |           0.4150 |     1.00 |      3.8136E-16 |
|        kahan |           0.5000 |     0.83 |      1.2734E-16 |
|        chunk |           0.2850 |     1.46 |      2.4869E-17 |
 
|      dot r32 | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    intrinsic |           0.2500 |     1.00 |      1.1426E-07 |
|        kahan |           0.2600 |     0.96 |      9.7811E-08 |
|        chunk |           0.1400 |     1.79 |      7.2122E-08 |
 
|      dot r64 | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    intrinsic |           0.2600 |     1.00 |      3.9246E-16 |
|        kahan |           0.3800 |     0.68 |      1.9229E-16 |
|        chunk |           0.1900 |     1.37 |      9.0927E-17 |

|        trigo | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|     fsin r32 |           0.0600 |   190.80 |      1.0325E-07 | 
|     fsin r64 |           0.0320 |   357.25 |      5.0118E-17 | 
|    facos r32 |           0.0280 |   221.43 |      1.0563E-06 | 
|    facos r64 |           0.0160 |   546.75 |      3.7996E-15 | 
|    fatan r32 |           0.0240 |   300.50 |      5.4993E-06 | 
|    fatan r64 |           0.0400 |   244.40 |      6.6869E-06 | 

|       hyperb | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|    ftanh r32 |           0.0280 |   510.71 |      5.5308E-08 | 
|    ftanh r64 |           0.0360 |   348.56 |      1.3282E-09 | 
|     ferf r32 |           0.0400 |   496.90 |      9.1205E-08 | 
|     ferf r64 |           0.0360 |   532.44 |      9.6298E-08 | 

|        rsqrt | <time> [ns/eval] | Speed-Up | relative error  |
|--------------|------------------|----------|-----------------|
|   frsqrt r32 |          16.3120 |     0.03 |      9.4387E-04 | 
|   frsqrt r64 |          16.7680 |     0.11 |      8.6745E-04 |
</details>

# Acknowledgement

* Compilation of this library was possible thanks to [Transvalor S.A.](https://www.transvalor.com/en/homepage) research activities. 
* Part of this library is based on the work of [Perini and Reitz](https://doi.org/10.1016/j.combustflame.2018.04.013), that was funded through the Sandia National Laboratories by the U.S. Department of Energy, Office of Vehicle Technologies, program managers Leo Breton, Gupreet Singh.
* The [fortran lang community](https://fortran-lang.discourse.group/) discussions such as [Some Intrinsic SUMS](https://fortran-lang.discourse.group/t/some-intrinsic-sums/5760) and [fastGPT](https://fortran-lang.discourse.group/t/fastgpt-faster-than-pytorch-in-300-lines-of-fortran/5385)

Contribution of open-source developers:

[jalvesz](https://github.com/jalvesz)

[perazz](https://github.com/perazz)
