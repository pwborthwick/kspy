from __future__ import division
from ks_grid import GRID
from ks_util import dii_subspace, out
import numpy as np

units = {'bohr->angstrom' : 0.52917721092}
elements = {'H':1, 'O':8}

class atom(object):
    #molecule class

    def __init__(self, symbol, center):

        self.symbol = symbol
        self.number = elements[symbol]
        self.center = np.array(center)

class orbital(object):
    #basis class

    def __init__(self, atom, momentum, exponents, coefficients, normalization, atoms):

        self.atom     = atom
        self.momentum = np.array(momentum)
        self.ex       = np.array(exponents)
        self.co       = np.array(coefficients)
        self.normal   = np.array(normalization)
        self.center   = atoms[atom].center

def nuclear_repulsion(mol):
    #nuclear repulsion

    eNuc = 0.0

    atoms = len(mol)
    for i in range(atoms):
        for j in range(i+1, atoms):
            r = np.linalg.norm(mol[i].center-mol[j].center)
            eNuc += mol[i].number*mol[j].number / r

    return eNuc

def lda_functional(name, rho):
    #exchange-correlation functional VWN3

    if name == 'VWN3':

        #slater exchange
        alpha = 2/3
        rho[rho == 0] = 1e-200
        Cx = -(9/8)*alpha*pow(3 /np.pi,1/3)
        ex = Cx*pow(rho,1/3)

        #VWN III parameterisation
        A= 0.0310907 ; x0=-0.409286 ; b= 13.0720 ; c = 42.7198
        rs = pow( 3.0 / ( 4.0 * np.pi * rho ) , 1.0/3.0 )
        zeta = 0.0
        x = np.sqrt(rs)
        X = lambda y: y*y + b*y + c
        Q = np.sqrt(4*c - b*b)

        ec = A * ( np.log( pow(x,2.0) / X(x) ) + 2.0 * b * np.arctan( Q/(2.0*x + b) ) * pow(Q,-1.0) -\
                  b * x0 * ( np.log( pow(x-x0,2.0) / X(x) ) + 2.0 * (b + 2.0 * x0) * np.arctan( Q/(2.0*x + b) ) * \
                  pow(Q,-1.0) ) * pow(X(x0),-1.0) ) 
        
        vc = ec - (A/3) * (c * (x-x0) - b*x*x0)/((x-x0)*(x*x + b*x + c))
        vx = 4*ex/3

        return ex + ec, vx + vc

def evaluate_gto(gto, p):
    #compute the value of gaussian density at (x,y,z)

    A = (p - gto.center) ; L = np.prod(A**gto.momentum, axis=1).reshape(-1,1)

    phi = np.sum(L*gto.normal*gto.co*np.exp(-gto.ex*np.sum(A*A, axis=1).reshape(-1,1)), axis=1)

    return phi.reshape(-1,1)

def evaluate_atomic_orbital(basis, p):
    #evaluate the GTO of the atomic orbitals of the molecule

    ao = []
    for i in basis:
        ao.append(evaluate_gto(i, p))

    return np.hstack(ao)

def evaluate_rho_lda(d, ao, weights):
    #evaluate the density over grid shells

    d = d + d.T
    ao_density = np.einsum('pr,rq->pq', ao, d, optimize=True)
    ao_density = np.einsum('pi,pi->p', ao, ao_density, optimize=True)

    #set small values to zeros
    ao_density[abs(ao_density) < 1.0e-15] = 0

    return ao_density

def evaluate_vxc(vxc, ao, weights):
    #construct exchange-correlation matrix

    weighted_ao = np.einsum('pi,p->pi', ao, 0.5*weights*vxc, optimize=True)
    xc = np.einsum('rp,rq->pq', ao, weighted_ao, optimize=True)

    return xc + xc.T

def evaluate_exc(exc, rho, weights):
    #evaluate exchange-correlation energy

    return np.einsum('p,p->', rho*weights, exc, optimize=True)

if __name__ == '__main__':

    mesh = 'close' ; functional = 'VWN3'; DIIS = True ; DIIS_SIZE = 6

    #define the molecule atoms first then basis (sto-3g)
    mol = []
    mol.append(atom('O', [0.0,0.0,0.0])) ; mol.append(atom('H', [0,-0.757 ,0.587])) ; mol.append(atom('H', [0,0.757,0.587]))
    for m in mol:
        m.center /= units['bohr->angstrom']

    orb = []
    orb.append(orbital(0, [0,0,0], [130.7093214, 23.80886605, 6.443608313], [0.1543289672962566, 0.5353281422870151, 0.44463454218921483],   \
                                   [27.551167822078394, 7.681819989204459, 2.882417873168662], mol))
    orb.append(orbital(0, [0,0,0], [5.033151319, 1.169596125, 0.38038896],  [-0.09996722918837482, 0.399512826093505, 0.7001154688886181],   \
                                   [2.394914882501622, 0.8015618386293724, 0.34520813393821864], mol))
    orb.append(orbital(0, [1,0,0], [5.033151319, 1.169596125, 0.38038896],  [0.15591627500155536, 0.6076837186060621, 0.39195739310391],     \
                                   [10.745832634231427, 1.7337440707285054, 0.4258189334467701], mol))
    orb.append(orbital(0, [0,1,0], [5.033151319, 1.169596125, 0.38038896],  [0.15591627500155536, 0.6076837186060621, 0.39195739310391],     \
                                   [10.745832634231427, 1.7337440707285054, 0.4258189334467701], mol))
    orb.append(orbital(0, [0,0,1], [5.033151319, 1.169596125, 0.38038896],  [0.15591627500155536, 0.6076837186060621, 0.39195739310391],     \
                                   [10.745832634231427, 1.7337440707285054, 0.4258189334467701], mol))
    orb.append(orbital(1, [0,0,0], [3.425250914, 0.6239137298, 0.168855404], [0.15432896729459913, 0.5353281422812658, 0.44463454218443965], \
                                   [1.7944418337900938, 0.5003264922111158, 0.1877354618463613], mol))
    orb.append(orbital(2, [0,0,0], [3.425250914, 0.6239137298, 0.168855404], [0.15432896729459913, 0.5353281422812658, 0.44463454218443965], \
                                   [1.7944418337900938, 0.5003264922111158, 0.1877354618463613], mol))

    #output details of molecule
    out([mol, DIIS, DIIS_SIZE, functional, mesh], 'initial')
    #use a reduced version of Harpy's cython integrals
    from ks_aello import aello
    s, t, v, eri = aello(mol, orb)

    #orthogonal transformation matrix
    from scipy.linalg import fractional_matrix_power as fractPow
    x = fractPow(s, -0.5)

    #inital fock is core hamiltonian
    h_core = t + v

    #orthogonal Fock
    fo = np.einsum('rp,rs,sq->pq', x, h_core, x, optimize=True )

    #eigensolve and transform back to ao basis
    eo , co = np.linalg.eigh(fo)
    c = np.einsum('pr,rq->pq', x, co, optimize=True)

    #build our initial density
    nocc = np.sum([a.number for a in mol])//2

    d = np.einsum('pi,qi->pq', c[:, :nocc], c[:, :nocc], optimize=True)

    #SCF conditions
    cycles = 50 ; tolerance = 1e-6
    out([cycles, tolerance], 'cycle')

    #get grid
    grid, weights = GRID(mol, mesh)

    #evaluate basis over grid
    ao = evaluate_atomic_orbital(orb, grid)

    last_cycle_energy = 0.0

    #diis initialisation
    if DIIS: diis = dii_subspace(DIIS_SIZE)

    #SCF loop
    for cycle in range(cycles):

        #build the coulomb integral
        j = 2.0 * np.einsum('rs,pqrs->pq', d, eri, optimize=True)

        #evalute density over mesh
        rho = evaluate_rho_lda(d, ao, weights)

        #evaluate functional over mesh
        exc, vxc = lda_functional(functional, rho)

        out([cycle, np.einsum('pq,pq->', d, (2.0*h_core), optimize=True), \
                    np.einsum('pq,pq->', d, ( j), optimize=True), \
                    evaluate_exc(exc, rho, weights),np.sum(rho*weights) ],'scf')

        #evaluate potential
        vxc = evaluate_vxc(vxc, ao, weights)

        f = h_core + j + vxc

        if (cycle != 0) and DIIS:
            f = diis.build(f, d, s, x)

        #orthogonal Fock and eigen solution
        fo = np.einsum('rp,rs,sq->pq', x, f, x, optimize=True )

        #eigensolve
        eo , co = np.linalg.eigh(fo)
        c = np.einsum('pr,rq->pq', x, co, optimize=True)

        #construct new density
        d = np.einsum('pi,qi->pq', c[:, :nocc], c[:, :nocc], optimize=True)

        #electronic energy
        eSCF = np.einsum('pq,pq->', d, (2.0*h_core + j), optimize=True) + evaluate_exc(exc, rho, weights)

        if abs(eSCF - last_cycle_energy) < tolerance: break
        if DIIS: vector_norm = diis.norm
        else:    vector_norm = ''
        out([cycle, abs(eSCF - last_cycle_energy), vector_norm],'convergence')

        last_cycle_energy = eSCF


    out([eSCF, np.einsum('pq,pq->', d, (2.0*h_core), optimize=True), \
                np.einsum('pq,pq->', d, ( j), optimize=True), \
                evaluate_exc(exc, rho, weights), nuclear_repulsion(mol) ], 'final')

    out([eo, c, np.sum(rho*weights), d, s, mol, orb], 'post')

