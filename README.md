# Fortran Fast math
A collection of functions for fast number crunshing using Fortran.

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

# TODO
* Contribution guidelines
* Finalize naming convention
* Integrate github CICD
* Polish autodoc

# Acknowledgement

Compilation of this library was possible thanks to [Transvalor S.A.](https://www.transvalor.com/en/homepage) research activities and, the contribution of open-source developers:

[jalvesz](https://github.com/jalvesz)