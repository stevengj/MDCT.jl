# MDCT module for Julia

[![Build Status](https://travis-ci.org/stevengj/MDCT.jl.svg?branch=master)](https://travis-ci.org/stevengj/MDCT.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/pl074ibwbl445tal?svg=true)](https://ci.appveyor.com/project/StevenGJohnson/mdct-jl)

[![MDCT](http://pkg.julialang.org/badges/MDCT_0.3.svg)](http://pkg.julialang.org/?pkg=MDCT&ver=0.3)
[![MDCT](http://pkg.julialang.org/badges/MDCT_0.4.svg)](http://pkg.julialang.org/?pkg=MDCT&ver=0.4)
[![MDCT](http://pkg.julialang.org/badges/MDCT_0.5.svg)](http://pkg.julialang.org/?pkg=MDCT&ver=0.5)
[![MDCT](http://pkg.julialang.org/badges/MDCT_0.6.svg)](http://pkg.julialang.org/?pkg=MDCT&ver=0.6)

This module computes the modified discrete cosine transform (MDCT) in
the Julia language and the inverse transform (IMDCT), using the fast
type-IV discrete cosine tranform (DCT-IV) functions in the
[FFTW.jl package](https://github.com/JuliaMath/FFTW.jl).

Definitions of the MDCT and IMDCT can be found, for example in the
[Wikipedia MDCT
article](http://en.wikipedia.org/wiki/Modified_discrete_cosine_transform).
The MDCT is a linear transformation that takes 2N inputs and produces
N outputs, which is designed to be applied to a sequence of
50%-overlapping blocks of a longer sequence (e.g. audio samples).
Because this is non-square (fewer outputs than inputs), the IMDCT is
not an "inverse" transformation in the usual sense; it only recovers
the original data when IMDCTs of overlapping blocks are added (by
"time-domain aliasing cancellation").

## Installation

Within Julia, just use the package manager to run
`Pkg.add("MDCT")` to install the files.

## Usage

To use the MDCT functions, simply do

    using MDCT
    Y = mdct(X)
    Z = imdct(Y)

where `X` is any numeric `AbstractVector` (1d array).  Currently, the
length of `X` must be a multiple of 4.

For example, suppose we make a random vector `X` of length 1000 and
consider 50%-overlapping blocks of length 100 (`X[1:100]`,
`X[51:150]`, `X[101:200]`, and so on).  If we perform the MDCT of two
such blocks, then the IMDCT, and then add the overlapping halves of
the IMDCT outputs, we recover that portion of the original data:

    X = rand(1000)
    Y1 = mdct(X[1:100])
    Y2 = mdct(X[51:150])
    Z1 = imdct(Y1)
    Z2 = imdct(Y2)
    norm(Z1[51:100] + Z2[1:50] - X[51:100])

where the last line computes the difference between the overlapped
IMDCT sum and the original data, and should be around
10<sup>&minus;15</sup> (floating-point roundoff error).

## Author

This module was written by [Steven G. Johnson](http://math.mit.edu/~stevenj/).
