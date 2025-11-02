# SPME — Smooth Particle Mesh Ewald in Fortran

A compact Fortran 2008 implementation of **Smooth Particle Mesh Ewald (SPME)** and **direct-space Ewald** with **OpenMP** and **FFTW3**. Suitable for research, benchmarking, and MSINDO coupling.

---

## Installation

**Requirements**
- Fortran compiler with OpenMP (e.g., `gfortran`, `ifort`, `nvfortran`)
- [FFTW3](http://www.fftw.org) in double precision

**Build**
```bash
make
