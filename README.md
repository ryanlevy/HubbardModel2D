# HubbardModel2D
A C++ version of the 2D Hubbard Model with arbitrary lattice. The Hubbard Model is defined by

<img src="images/hubbardmodel.png" width="60%">

## Requirements
* C++11
* Eigen3
* For diagonalization 
  * By default the code uses Lanczos using [ietl](https://github.com/garrison/ietl)
    * ietl requires Boost and LAPACK
  * For more exact eigenvalue/vector calculations there is [Spectra](https://github.com/yixuan/spectra). See the note below.


## Notes
Defaults to 4x4 lattice with 2 ↑ and 2 ↓ electrons at U/t=4. The boundary conditions are set by the bond file implicitly. 
Unfortunately the code doesn't take into account symmetries; therefore it doesn't scale well.

To confirm data see [E. Dagotto, A. Moreo, F. Ortolani, D. Poilblanc, J. Riera, Phys. Rev. B 45 (1992) 10741–10760.](https://link.aps.org/doi/10.1103/PhysRevB.45.10741)
The default values, `Ne=4` and `U/t=4` should have a ground state of `-11.53029`

### Input
There is not currently a default lattice, so a list of bonds is required. Place `bond.dat` in the same folder as the executable when running. Two examples of 4x1 (how to have the code calculate 1D) and 4x4 are included. The number of sites*2 needs to be set as the BASISSIZE in `Hubbard2D.hpp`

One can programatically change the number of electrons and the `U/t` value using the following command line input:
```
$ ./HubbardModel2D -ne [#num electrons] -U [# for U/t]
```
Thefore the default numbers are `./HubbardModel2D -ne 2 -U 4`

### Spectra
---
If you use spectra, this code requires a patch of spectra, (for non lanczos)
```
cp eigenref.patch /path/to/spectra/
git apply eigenref.patch
```
