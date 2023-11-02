# EZmock

![GitHub](https://img.shields.io/github/license/cheng-zhao/EZmock.svg)
![Codacy grade](https://img.shields.io/codacy/grade/2952c648b1934e49bd1cb8acfb4fc1fa.svg)

## Table of Contents

-   [Introduction](#introduction)
-   [Compilation](#compilation)
-   [References](#references)

## Introduction

**E**ffective **Z**el'dovich approximation **mock** generator (EZmock<sup>[\[1\]](#ref1)</sup>) is a technique for constructing approximate mock catalogues for large-scale structures of the Universe. It relies on the Zel'dovich approximation to generate the dark matter density field at a desired redshift, and then populates tracers (e.g., galaxies and quasars) based on an effective bias description. The product is also dubbed &ldquo;EZmock&rdquo; (plural: EZmocks).

EZmock is among the fastest methods that are able to reproduce the clustering of dark matter haloes from an *N*-body simulation with &#8818;&thinsp;5% precision down to the scale of *k* &sim; 0.4&thinsp;*h*&thinsp;Mpc<sup>&minus;1</sup> (see [\[2\]](#ref2) for a comparison of mock generation methods). It was used to produce the 1000 multi-tracer mock catalogues for the final extended Baryon Oscillation Spectroscopic Survey (eBOSS) data release<sup>[\[3\]](#ref3)</sup>.

This repository consists of 3 interfaces of the EZmock algorithm:
-   A C library with APIs for the construction of EZmocks;
-   A python package that wraps the C library functions;
-   A C program that generates EZmocks using the core library functions.

All the interfaces are parallelised with the shared-memory Open Multi-Processing ([OpenMP](https://www.openmp.org)), and rely on the [FFTW](http://www.fftw.org) library for Fast Fourier Transforms (FFT). The C codes are compliant with the ISO C99 and IEEE POSIX.1-2008 standards. The software packages are written by Cheng Zhao (&#36213;&#25104;) based on the initial [fortran version](https://github.com/chia-hsun-chuang/ezmock) developed by [Dr. Chia-Hsun Chuang](https://github.com/chia-hsun-chuang).

As a whole, the software packages in this repository are released under the [GPLv3 license](LICENSE.txt), due to its reliance on the FFTW library, though many of the source files are distributed under the [MIT license](LICENSE_MIT.txt), as indicated in their header lines. If you use this tool in research work that results in publications, please cite the papers [\[1\]](#ref1) and [\[3\]](#ref2).

## Compilation

The building of the C interfaces of EZmock is based on the make utility. Customisable compilation options can be set in the file [`options.mk`](options.mk).

Once the setting is done, the following command should compile the C program and the C library:

```bash
make
```

Alternatively, different interfaces can be compiled via:
```bash
make EZmock        # compile the C program
make libEZmock.so  # compile the dynamic C library
make libEZmock.a   # compile the static C library
```

Once the C library is compiled, it can be installed at the specific location with
```bash
make install
```

The following commands clear the compiled files:
```bash
make clean     # clear the C program
make cleanlib  # clear the C libraries
make cleanall  # clean everything
```

If the C libraries are installed in the root directory of this repository, the python package can be built in the [`python`](python/) folder via
```bash
python setup.py build_ext -i
```

Do not forget to add the path with `libEZmock.so` to your environment variable `LD_LIBRARY_PATH` before importing the python package.


## References

<span id="ref1">\[1\]</span> Chuang C.-H., Kitaura F.-S., Prada F., Zhao C., Yepes G., 2015, [EZmocks: extending the Zel'dovich approximation to generate mock galaxy catalogues with accurate clustering statistics](https://doi.org/10.1093/mnras/stu2301), *MNRAS*, 446(3), 2621&ndash;2628. \[[arXiv:1409.1124](https://arxiv.org/abs/1409.1124)\] \[[ADS Abstract](https://ui.adsabs.harvard.edu/abs/2015MNRAS.446.2621C)\]

<span id="ref2">\[2\]</span> Chuang C.-H., et al., 2015, [nIFTy cosmology: Galaxy/halo mock catalogue comparison project on clustering statistics](https://doi.org/10.1093/mnras/stv1289), *MNRAS*, 452(1), 686&ndash;700. \[[arXiv:1412.7729](https://arxiv.org/abs/1412.7729)\] \[[ADS Abstract](https://ui.adsabs.harvard.edu/abs/2015MNRAS.452..686C)\]

<span id="ref3">\[3\]</span> Zhao C., et al., 2021, [The completed SDSS-IV extended Baryon Oscillation Spectroscopic Survey: 1000 multi-tracer mock catalogues with redshift evolution and systematics for galaxies and quasars of the final data release](https://doi.org/10.1093/mnras/stab510), *MNRAS*, 503(1), 1149&ndash;1173. \[[arXiv:2007.08997](https://arxiv.org/abs/2007.08997)\] \[[ADS Abstract](https://ui.adsabs.harvard.edu/abs/2021MNRAS.503.1149Z)\]

