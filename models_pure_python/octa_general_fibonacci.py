# octahedron model
# Note: model title and parameter table are inserted automatically
r"""

General octahedron shape with truncature.

Pure python model with orientation average computed using Lebedev quadrature on the unit sphere.
Need installation of Lebedev library: 
https://pypi.org/project/pylebedev/ 
using: pip install pylebedev

This model contains the 3D amplitude of the Form Factor amp3D(qx,qy,qz).
And the 3D form factor intensity int3D(qx,qy,qz). Scaled to 1 at q close to zero.
The average form factor Iq(q) is computed from amp3D(qx,qy,qz) using Lebedev quadrature.
(qx,qy,qz) are the 3D coordinates of the scattering vectr and q its modulus.

Singularities of the 3D amplitude form factor are encoutered when:
Qx=QY, QY=QZ or QZ=QX
They are not numerically corrected here yet !!!
Division by zero is avoided simply by making sure that the quadrature points do not match some singularity direction.
Because Lebedev quadrature is based on octagonal symmetry, it contains points along the symmetry axis of the ocatahedron shape.
In practice, this is done by applying two small consecutive rotations (around z-axis and x-axis) of all the quadrature points coordinates.
Works well for a regular octahedron. 
To be tested when b2a_ratio and c2a_ratio are different from one. 
----------

The general octahedron is defined by three distances (length_a, length_b and length_c) from the center of the octahedron to its six
vertices. 
Coordinates of the six vertices are : 
(length_a,0,0) (-length_a,0,0) 
(0,length_b,0) (0,-length_b,0) 
(0,0,length_c) (0,0,-length_c)

In addition, truncation at all vertices by a square facet is introduced, with a truncature parameter t.

.. math::

The model is using 4 parameters for the general octahedron shape.
length_a, the two ratios b2a_ratio and c2a_ratio and t:
 b2a_ratio = length_b/length_a
 c2a_ratio = length_c/length_a
 1/2 < t < 1

Volume of the truncated general octahedron is:
V = (4/3)*length_a^3*b2a_ratio)*c2a_ratio*(1-3*(1-t)^3)

Lengths of edges for the general octahedron are :
    A_edge^2 = length_a^2+length_b^2
    B_edge^2 = length_a^2+length_c^2
    C_edge^2 = length_b^2+length_c^2

A regular octahedron corresponds to :
b2a_ratio = c2a_ratio = 1

Truncature parameter is between 1/2 and 1.
For a regular (no truncaton) octahedron :
    t=1
    
For a cuboctahedron :
    t=1/2


For a regular truncated octahedron :
    b2a_ratio = c2a_ratio = 1   
    A_edge = B_edge = C_edge = length_a*sqrt(2)
    length_a = length_b = length_c = A_edge/sqrt(2)
    V = (4/3) * length_a^3 = (sqrt(2)/3) * A_edge^3  * (1-3(1-t)^3)

Amplitude of the 3D form factor AP is calculated using a dimensionless scattering vector (Qx,Qy,Qz):
    qx = Qx * length_a
    qy = Qy * length_b
    qz = Qz * length_c

normalisation to 1. of AP at the limit Qx=Qy=Qz=0. Division by a Factor 4/3.
                
    AP = 6./(1.-3(1.-t)^3)*(AA+BB+CC);

    AA, BB and CC derive from each other by a circular permutation of QX,QY,QZ:

    AA = 1./(2*(Qy*Qy-Qz*Qz)*(Qy*Qy-Qx*Qx))*((Qy-Qx)*sin(Qy*(1.-t)-Qx*t)+(Qy+Qx)*sin(Qy*(1.-t)+Qx*t))+
        1./(2*(Qz*Qz-Qx*Qx)*(Qz*Qz-Qy*Qy))*((Qz-Qx)*sin(Qz*(1.-t)-Qx*t)+(Qz+Qx)*sin(Qz*(1.-t)+Qx*t));

    BB = 1./(2*(Qz*Qz-Qx*Qx)*(Qz*Qz-Qy*Qy))*((Qz-Qy)*sin(Qz*(1.-t)-Qy*t)+(Qz+Qy)*sin(Qz*(1.-t)+Qy*t))+
        1./(2*(Qx*Qx-Qy*Qy)*(Qx*Qx-Qz*Qz))*((Qx-Qy)*sin(Qx*(1.-t)-Qy*t)+(Qx+Qy)*sin(Qx*(1.-t)+Qy*t));

    CC = 1./(2*(Qx*Qx-Qy*Qy)*(Qx*Qx-Qz*Qz))*((Qx-Qz)*sin(Qx*(1.-t)-Qz*t)+(Qx+Qz)*sin(Qx*(1.-t)+Qz*t))+
        1./(2*(Qy*Qy-Qz*Qz)*(Qy*Qy-Qx*Qx))*((Qy-Qz)*sin(Qy*(1.-t)-Qz*t)+(Qy+Qz)*sin(Qy*(1.-t)+Qz*t));

       

$\rho_\text{p}$ is the scattering length inside the volume, $\rho_\text{solvent}$
is the scattering length of the solvent, and (if the data are in absolute
units) *scale* represents the volume fraction (which is unitless).


Validation
----------

Validation of the code is made by comparing with the results of octahedron_truncated model available here:
https://marketplace.sasview.org/models/152/octahedron.
This model performs the orientation average with ancillary c code using Gauss-Legendre quadrature.

References
----------
1. A. L. Patterson. “The Diffraction of X-Rays by Small Crystalline Particles”. In: Phys. Rev. 56 (10 Nov.
1939), pp. 972–977. doi: 10.1103/PhysRev.56.972.

2. J. Schelten R. W. Hendricks and W. Schmatz. “Studies of voids in neutron-irradiated aluminium single
crystals: II. Small-angle neutron scattering”. In: The Philosophical Magazine: A Journal of Theoretical
Experimental and Applied Physics 30.4 (1974), pp. 819–837. doi: 10.1080/14786437408207237.

3. Wei-Ren Chen et al. “Scattering functions of Platonic solids”. 
In:Journal of Applied Crystallography - J APPL CRYST44 (June 2011).
DOI:10.1107/S0021889811011691


Authorship and Verification
----------------------------

* **Authors: Marianne Imperor-Clerc (marianne.imperor@cnrs.fr)
* **Last Modified by MI: 2d September 2025**
* **Last Reviewed by:**
"""

import numpy as np
import time
from numpy import inf
import matplotlib.pyplot as plt



name = "octa_general_fibonacci"
title = "General octahedron pure python model using Fibonacci quadrature"
description = """
            AP = 6./(1.-3(1.-t)**3)*(AA+BB+CC);
            AA = 1./(2*(Qy*Qy-Qz*Qz)*(Qy*Qy-Qx*Qx))*((Qy-Qx)*sin(Qy*(1.-t)-Qx*t)+(Qy+Qx)*sin(Qy*(1.-t)+Qx*t))+
                                1./(2*(Qz*Qz-Qx*Qx)*(Qz*Qz-Qy*Qy))*((Qz-Qx)*sin(Qz*(1.-t)-Qx*t)+(Qz+Qx)*sin(Qz*(1.-t)+Qx*t));

            BB = 1./(2*(Qz*Qz-Qx*Qx)*(Qz*Qz-Qy*Qy))*((Qz-Qy)*sin(Qz*(1.-t)-Qy*t)+(qz+qy)*sin(qz*(1.-t)+Qqy*t))+
                                1./(2*(Qx*Qx-Qy*Qy)*(Qx*Qx-Qz*Qz))*((Qx-Qy)*sin(Qx*(1.-t)-Qy*t)+(Qx+Qy)*sin(Qx*(1.-t)+Qy*t));

            CC = 1./(2*(Qx*Qx-Qy*Qy)*(Qx*Qx-Qz*Qz))*((Qx-Qz)*sin(Qx*(1.-t)-Qz*t)+(Qx+Qz)*sin(Qx*(1.-t)+Qqz*t))+
                                1./(2*(Qy*Qy-Qz*Qz)*(Qy*Qy-Qx*Qx))*((Qy-Qz)*sin(Qy*(1.-t)-Qz*t)+(Qy+Qz)*sin(Qy*(1.-t)+Qz*t));

normalisation to 1. of AP at q = 0. Division by a Factor 4/3.

    	    Qx = qx * length_a;
    	    Qy = qy * length_b;
    	    Qz = qz * length_c;
    	    0 < t < 1
            
"""
category = "shape:parallelepiped"

#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["sld", "1e-6/Ang^2", 130, [-inf, inf], "sld",
               "Octahedron scattering length density"],
              ["sld_solvent", "1e-6/Ang^2", 9.4, [-inf, inf], "sld",
               "Solvent scattering length density"],
              ["length_a", "Ang", 500, [0, inf], "volume",
               "half height along a axis"],
              ["b2a_ratio", "", 1, [0, inf], "volume",
               "Ratio b/a"],
              ["c2a_ratio", "", 1, [0, inf], "volume",
               "Ratio c/a"],
              ["t", "", 0.99, [0.0, 1.0], "volume",
               "truncation along axis"],
              ["npoints_fibonacci",     "",           300, [1, 1e6],   "", 
                       "Number of points on the sphere for the Fibonacci integration"]
              ]


def volume(length_a, b2a_ratio, c2a_ratio, t):
    return (4/3)*length_a**3*b2a_ratio*c2a_ratio*(1.-3*(1.-t)**3)

# def A(qx,qy,qz,length_a, b2a_ratio, c2a_ratio, t): # Ancillary function. Term AA in the amplitude.
#     length_b=length_a*b2a_ratio
#     length_c=length_a*c2a_ratio
#     qnx=qx*length_a # conversion to dimensionless coordinate
#     qny=qy*length_b # conversion to dimensionless coordinate
#     qnz=qz*length_c # conversion to dimensionless coordinate
    
#     AA = 1./((qny*qny-qnz*qnz)*(qny*qny-qnx*qnx))*((qny-qnx)*np.sin(qny*(1.-t)-qnx*t)+(qny+qnx)*np.sin(qny*(1.-t)+qnx*t))+1./((qnz*qnz-qnx*qnx)*(qnz*qnz-qny*qny))*((qnz-qnx)*np.sin(qnz*(1.-t)-qnx*t)+(qnz+qnx)*np.sin(qnz*(1.-t)+qnx*t))
#     return AA

def _rotate_unit_vectors(u, tetaz=0.01076, tetax=0.0203):
    """Small global rotation to avoid exact symmetry directions (safe default)."""
    cz, sz = np.cos(tetaz), np.sin(tetaz)
    cx, sx = np.cos(tetax), np.sin(tetax)
    Rz = np.array([[cz, -sz, 0.0],
                   [sz,  cz, 0.0],
                   [0.0, 0.0, 1.0]])
    Rx = np.array([[1.0, 0.0, 0.0],
                   [0.0, cx, -sx],
                   [0.0, sx,  cx]])
    R = Rx.dot(Rz)
    return u.dot(R.T)

def A(qx,qy,qz,length_a, b2a_ratio, c2a_ratio, t): # Ancillary function. Term AA in the amplitude.
    length_b=length_a*b2a_ratio
    length_c=length_a*c2a_ratio
    qnx=qx*length_a # conversion to dimensionless coordinate
    qny=qy*length_b # conversion to dimensionless coordinate
    qnz=qz*length_c # conversion to dimensionless coordinate

    # protect denominators against exact zeros by adding tiny epsilon
    eps_den = 1e-18
    denom1 = (qny*qny - qnz*qnz) * (qny*qny - qnx*qnx)
    denom2 = (qnz*qnz - qnx*qnx) * (qnz*qnz - qny*qny)
    denom1 = denom1 + eps_den * (np.abs(denom1) < eps_den)  # add epsilon where denom is zero
    denom2 = denom2 + eps_den * (np.abs(denom2) < eps_den)  # add epsilon where denom is zero

    term1 = ((qny - qnx) * np.sin(qny*(1.-t) - qnx*t) + (qny + qnx) * np.sin(qny*(1.-t) + qnx*t)) / denom1
    term2 = ((qnz - qnx) * np.sin(qnz*(1.-t) - qnx*t) + (qnz + qnx) * np.sin(qnz*(1.-t) + qnx*t)) / denom2

    AA = term1 + term2
    return AA

def amp3D(qx,qy,qz,length_a, b2a_ratio, c2a_ratio, t): 
    # code taking into account the circular permutation between AA, BB and CC terms in AP
    # amp3D is scaled to 1 near the origin for a regular shape
    AA = A(qx,qy,qz,length_a, b2a_ratio, c2a_ratio, t)

    BB = A(qy,qz,qx,length_a, b2a_ratio, c2a_ratio, t)

    CC = A(qz,qx,qy,length_a, b2a_ratio, c2a_ratio, t)

	# normalisation to 1. of AP at q = 0. Division by a Factor 4/3.
    # and global 1/2 coefficient
    AP = 3./(1.-3*(1.-t)**3)*(AA+BB+CC)

    return AP 

def int3D(qx,qy,qz,length_a, b2a_ratio, c2a_ratio, t): # code taking into account the circular permutation between AA, BB and CC terms in AP
    # int3D is scaled to 1 near the origin for a regular shape
    amp = amp3D(qx,qy,qz,length_a, b2a_ratio, c2a_ratio, t)
    intensity = amp**2

    return intensity




######################################## Fibonacci integration ################################################

def fibonacci_sphere(npoints_fibonacci: int):
    """
    Generates npoints quasi-uniformly distributed on the unit sphere
    in Cartesian coordinates (x,y,z) and their associated weights.
    Parameters
    ----------
    npoints : int
        Number of points to generate.
    Returns
    -------
    points : ndarray, shape (npoints, 3)
        Cartesian coordinates of the points on the unit sphere.
    weights : ndarray, shape (npoints,)
        Weights associated with each point for integration on the sphere.
    """

    indices = np.arange(0, npoints_fibonacci, dtype=float) + 0.5
    phi = np.arccos(1 - 2*indices/npoints_fibonacci)        
    theta = 2 * np.pi * indices / ((1 + 5**0.5) / 2)  
    
    x = np.cos(theta) * np.sin(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(phi)

    weight = np.full(len(x),1/npoints_fibonacci)
    
    return np.column_stack((x, y, z)), weight


def plot_fibonacci_sphere(npoints_fibonacci=500, figsize=(7,7)):
    """
    3D representation of Fibonacci points on the unit sphere.
    Parameters
    ----------
    npoints : int
        Number of points to generate and display.
    figsize : tuple
        Size of the figure.

    """
    pts,w = fibonacci_sphere(npoints_fibonacci)

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(pts[:,0], pts[:,1], pts[:,2], s=10, alpha=0.6)

    # Axes et esthétique
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_title(f"Fibonacci points on the unit sphere, ({npoints_fibonacci} total points)")
    ax.set_box_aspect([1,1,1])  # sphère non déformée
    plt.show()


# def Iq(q, sld, sld_solvent, nsides:int, Rave, L, npoints_fibonacci:int=1000):
#     """
#     Integration over the unit sphere with Fibonacci points.
#     Each point has an equal weight = 1/npoints.

#     Parameters
#     ----------
#     q : float ou array
#         Norm of the scattering vector
#     sld, sld_solvent : 
#         Contrast of scattering length density
#     nsides, Rave, L : 
#         Geometrical parameters of the prism
#     npoints : int
#         Number of Fibonacci points on the sphere
#     Returns
#     -------
#     Iq : ndarray
#         Scattering intensity averaged over all orientations
#     time_fib : float
#         Execution time of the function in seconds
#     n_points_total : int
#         Total number of points used in the quadrature
#     """

#     time_start = time.time()
#     n_points_total = npoints_fibonacci
#     nsides = int(nsides)
#     q = np.atleast_1d(q)  # vecteur q

    
#     q_unit,w = fibonacci_sphere(npoints_fibonacci)   # shape (npoints, 3)

#     # Projections
#     qa = q[:, np.newaxis] * q_unit[:, 0][np.newaxis, :]
#     qb = q[:, np.newaxis] * q_unit[:, 1][np.newaxis, :]
#     qc = q[:, np.newaxis] * q_unit[:, 2][np.newaxis, :]

#     # Compute intensity
#     intensity = Iqabc(qa, qb, qc, nsides, Rave, L)  # shape (nq, npoints)

#     # Uniform average over the sphere
#     integral = np.sum(w[np.newaxis, :] * intensity, axis=1)
#     integral = np.mean(intensity, axis=1)
#     time_end = time.time()
#     time_fib = time_end - time_start
   
#     print(f'Execution time Fibonacci with {npoints_fibonacci} points: {time_fib:.4f} seconds')

#     return integral * (sld - sld_solvent)**2


# Iq.vectorized = True



def Iq(q, sld, sld_solvent,length_a=500, b2a_ratio=1, c2a_ratio=1, t=0.99,  npoints_fibonacci:int=400, rotate:bool=True):
    """
    Parameters:
        q:      input scattering vectors, units 1/Ang
    Returns:
        I(q):   1D scattering intensity at q, units 1/cm
    """


    time_start = time.time()
    q = np.atleast_1d(q)  # vecteur q
    q_unit,w = fibonacci_sphere(npoints_fibonacci)   # shape (npoints, 3)
    if rotate:
        q_unit = _rotate_unit_vectors(q_unit)

    # # w are the weights and q_unit are the quadrature points on the unit sphere
    # # multiplication by q, the modulus of the scattering vector 
    # qquad_x=q*q_unit[:,0] # all x coordinates of the quadrature points multiplied by q
    # qquad_y=q*q_unit[:,1] # all y coordinates of the quadrature points multiplied by q
    # qquad_z=q*q_unit[:,2] # all z coordinates of the quadrature points multiplied by q

    # tetaz=0.01076 # add a small rotation around z-axis to avoid points having qx=qy, qy=qz or qx=qz
    # qrotz_x=np.cos(tetaz)*qquad_x-np.sin(tetaz)*qquad_y
    # qrotz_y=np.sin(tetaz)*qquad_x+np.cos(tetaz)*qquad_y
    # qrotz_z=qquad_z

    # tetax=0.0203 # add a second small rotation around x-axis to avoid points having qx=qy, qy=qz or qx=qz
    # qx=qrotz_x
    # qy=np.cos(tetax)*qrotz_y-np.sin(tetax)*qrotz_z
    # qz=np.sin(tetax)*qrotz_y+np.cos(tetax)*qrotz_z

    # build qx,qy,qz arrays with correct broadcasting -> shape (nq, npoints)
    qx = q[:, None] * q_unit[None, :, 0]
    qy = q[:, None] * q_unit[None, :, 1]
    qz = q[:, None] * q_unit[None, :, 2]

    # compute intensity grid using existing amp3D (vectorized)
    amp_grid = amp3D(qx, qy, qz, length_a, b2a_ratio, c2a_ratio, t)  # shape (nq, npoints)
    intensity_grid = np.abs(amp_grid)**2

    integral = np.sum(intensity_grid * w[None, :], axis=1) #summation over all points

    time_end = time.time()
    time_fib = time_end - time_start
    print(f'Execution time Fibonacci with {npoints_fibonacci} points: {time_fib:.4f} seconds')


    return integral*0.0001*(sld-sld_solvent)**2*volume(length_a, b2a_ratio, c2a_ratio, t)


Iq.vectorized = True # Here Iq works only for a single float value of q parameter