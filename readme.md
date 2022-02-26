![image](https://user-images.githubusercontent.com/73105740/142720205-ddb5a6ad-4d6c-4a1a-8171-9ea1af5eece4.png)
## Density Function Theory Program - kspy
  In the 1970's when I was a post-doctoral student DFT wasn't talked about as a quantum chemistry method, so it's only now with a renewed interest in the topic that it's come to my attention. I wanted to write a DFT program in the hope that along the way I'd learn something about it. The problem with writing a DFT program is the same as with writing an HF program - evaluating the integrals. With HF you will probably end up using either McMurchie-Davidson or Obara-Saika recurrance schemes (See Josh Going's blog article [A (Hopefully) gentle guide to computer implementation of molecular integrals](https://joshuagoings.com/assets/integrals.pdf)), for DFT we're going to have to use numerical methods. This means integrating on grid of points...

### The Grid
  I could have used a cube as a space grid and then taken Riemann sums to evaluate (there's a great YouTube series by James Johns where he develops a matlab HF program and shows how to convert it to DFT. In this he uses Riemann sums to evaluate the integrals in the DFT program.) However, I decided to try for a 'proper' atom centered spherical grid approach. A useful paper was [PMW Gill, BG Johnson and JA Poples 'A standard grid for density functional theory'](https://rsc.anu.edu.au/~pgill/papers/038SG1.pdf), although I didn't use this SG-1 grid the paper helped understand the techniques involved. The grid I settled on was a coarse grid of (10,15) radial points for period 1 and period 2 elements respectively. The radial grid is a Mura-Knowles radial grid [ME Mura and PJ Knowles 'Improved radial grids for quadrature in density-functional calculations' JCP 104, 9848 (1996); DOI:10.1063/1.471749](https://aip.scitation.org/doi/10.1063/1.471749). The 'coarse' angular grid is of Lebedev orders (11, 15) for period 1 and period 2 respectively. This translates into 50 and 86 points respectively arranged on a spherical shell (VI Lebedev, and DN Laikov, Doklady Mathematics, 'A Quadrature formula for the sphere of the 131st algebraic order of accuracy' Vol. 59, No. 3, (1999)). There are various sources for this data given in the external links of the wikipedia article on Lebedev integration.  A pruning scheme is employed to systematically reduce the number of angular points in regions where dense angular quadrature is not necessary, such as near the nuclei where the charge density is approximately spherically symmetric and at long distance from the nucleus. The pruning scheme I employed was the Treutler-Aldrich scheme [O Treutler and R Ahlrich, 'Efficient molecular numerical integration schemes' JCP 102, 346 (1995); DOI:10.1063/1.469408](https://aip.scitation.org/doi/pdf/10.1063/1.469408). The partitioning of the atomic centered grids to a molecular grid follows a Becke scheme after Stratmann [RE Stratmann, GE Scuseria and MJ Frisch, 'Achieving Linear scaling in exchange-correlation density functional quadratures' CPL 257, 3-4 (1996); DOI:10.1016/009-2614(96)00600-8](https://www.sciencedirect.com/science/article/abs/pii/0009261496006008?via%3Dihub). Finally I have implemented a final radius adjustment during the partition (Becke suggests doing this) using the Bragg radius. A second 'close' grid is also included which is a (50, 75) radial and (29, 29) angular, the latter representing 302 points on each shell. The grid routines are in ks_grid.py.

### The HF Integrals
  To get the DFT SCF started we need an initial density. To do this I use a HF overlap matrix S, and an initial Fock matrix composed of the sum of the 1-electron kinetic and coulomb integrals (core Hamiltonian - T+V). This Fock is then orthogonalised (F<sup>'</sup>) as (S<sup>-0.5</sup>)<sup>T</sup>FS<sup>-0.5</sup>, eigensolve the resulting orthogonal Fock for orbital coefficients C orthogonal, transform back to atomic basis as S<sup>-0.5</sup>C<sup>'</sup>, use resulting ao coefficients to compute a density matrix D<sub>&mu;&nu;</sub> = c<sub>&mu;i</sub>c<sub>i&nu;</sub> where i is over occupied orbitals. This initial density can be used with initial Fock and 2-electron repulsion integrals to form the coulomb integral J (we don't want the HF exchange integral K for DFT). To get these integrals I've used a modified version of Harpy's Cython integral package *aello*. I've removed the angular and electric field integrals and additionally the 2-electron repulsions are returned as a tensor rather than a linear array. These are in ks_aello.pyx.

### Molecule and Basis Sets
  The molecule definition is contained in a *mol* object which is itself comprised of objects from an atom class. Each instance of the atom class contains the atom symbol, atomic number and the coordinates of the atom center (array[3]). The molecule is hard coded as H<sub>2</sub>O. The basis is contained in an *orb* object which is itself comprised of objects from a gaussian class. Each instance of the gaussian class contains the atom the Gaussian is centered on, the momentum(array[3]), the exponents (array[primatives]), the coefficients (array[primatives]), the normalisation (array[primatives]) and a copy of the atom center coordinates (array[3]). The momenta are given as s [0,0,0] px [1,0,0] py [0,1,0] and pz [0,0,1]. The basis used is a simple STO-3G so we only require s and p orbitals. The primatives exponent and coefficient values are hard-coded in the __main__ section. (I use the psi4 format of the basis sets from BSE which have some (small) differences from the nwchem format versions as used by eg pyscf. This might lead to numerical differences in values when using high precision).

### The Functional
  We are using the LDA (Local Density Approximation) with the VWN III - defined [here](http://user.it.uu.se/~anakr367/files/courses/elec_struc/report.pdf)

### The Density
  The density is determined by first evaluating each Gaussian orbital over the grid (&rho;<sub>g</sub>), then assembling those into a complete basis over the grid (&rho;<sub>b</sub> = &Sigma;&rho;<sub>g</sub>) and finally computing &rho; as the product  &rho;<sub>b</sub>*&rho;<sub>b</sub>D. 

### Convergence
  The first runs of the program without diis resulted in an oscillating scf which failed to converge. Implementing a diis scheme resulted in good convergence. The diis scheme is taken from harpy (diis.py) with some minor changes. This is implemented as a class in ks_util.

### Output
        ks output
    molecule is             water
    geometry is             OH bonds - 0.96      HOH angle - 104.42
    basis is                STO-3G {psi4 format}
    analytic integration    aello cython - McMurchie-Davidson scheme
    numerical integration   (Mura-Knowles, Lebedev)
                            radial prune is Aldrich-Treutler
                            Becke Partition scheme is Stratmann 
                            Radial adjust is Treutler
                            order: period 1 (10,11) and period 2 (15,15)
    mesh                    close
    functional              VWN3
    diis                    True   buffer size is  6

    scf control             maximum cycles are  50         convergence tolerance  1e-06

     cycle     1 electron         coulomb         exchange          electrons
                                                                                       ΔE         diis norm
    -------------------------------------------------------------------------------------------------------
        0     -127.35805722     54.98251023     -9.97719421         10.0000 

        1     -117.49240988     42.82393549     -8.69633734         10.0000 
                                                                                     8.110350      0.928093 
        2     -125.72307277     51.61343056     -9.53620544         10.0000 
                                                                                     5.270233      0.764866 
        3     -122.66235858     47.64875083     -9.10024936         10.0000 
                                                                                     1.324036      0.060737 
        4     -122.43605098     47.39593167     -9.07645589         10.0000 
                                                                                     0.089054      0.007794 
        5     -122.40701219     47.36370308     -9.07330448         10.0000 
                                                                                     0.011201      0.001079 
        6     -122.40308403     47.35935264     -9.07288295         10.0000 
                                                                                     0.001425      0.000168 
        7     -122.40234724     47.35853669     -9.07280381         10.0000 
                                                                                     0.000335      0.000003 
        8     -122.40236065     47.35855154     -9.07280525         10.0000 
                                                                                     0.000006      0.000000 
        9     -122.40236065     47.35855154     -9.07280525         10.0000 

    final energies (Hartree)
    ------------------------
    one electron        -122.4023606488
    coulomb               47.3585515443  
    exchange              -9.0728052453  
    nuclear repulsion      9.1882584177   
    total electronic     -84.1166143498 

    final total energy   -74.9283559321 

### Test
  The python code for the functional used here was copied into a pyscf eval_xc subroutine (user defined functional) and run with the same grid. The pyscf value is -74.92835585398149 Hartree <sup>*</sup>, this agrees better than the convergence criteria of 1e-6. The value using libxc from pyscf is -74.92835585398149 Hartree.
  
  **<sup>*</sup>** the basis definition with pyscf (nwchem) differs slightly from the basis here (psi4).

### Installation
  Copy all files to a directory, run

  	python3 setup.py build_ext --inplace install --user

then 

    python3 ks_main.py

    
### Additions 
12 December 2021 - added some simple post SCF properties.

    molecular orbitals
    ------------------
     0  -18.29159 occupied  
     1   -0.85022 occupied  
     2   -0.40199 occupied  
     3   -0.16895 occupied  
     4   -0.07683 homo      
     5    0.29772 lumo      
     6    0.40637 virtual   

    mulliken populations
    --------------------
     0    1.99686   1s  
     1    1.84808   2s  
     2    2.00000   2px 
     3    1.08764   2py 
     4    1.45282   2pz 
     5    0.80729   1s  
     6    0.80729   1s  

    atomic charge
    ---------------
     0   -0.38542 O   
     1    0.19271 H   
     2    0.19271 H   

    dipole momemts (Debye)
    ----------------------
     x= 0.00000  y= 0.00000  z= 1.73423 

