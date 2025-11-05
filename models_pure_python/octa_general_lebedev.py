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
from numpy import inf
from pylebedev import PyLebedev # class in the pyLebedev library

import matplotlib.pyplot as plt
import time

# List of the different quadrature orders in Lebedev
orderlist=[3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 35, 41, 47, 53, 59, 65, 71, 77, 83, 89, 95, 101, 107, 113, 119, 125, 131]

leblib = PyLebedev()

# list of the different quadrature orders in Lebedev
# orderlist goes from 0 to 31
orderlist=[3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 35, 41, 47, 53, 59, 65, 71, 77, 83, 89, 95, 101, 107, 113, 119, 125, 131]

name = "octa_general_lebedev"
title = "General octahedron pure python model using Lebedev quadrature"
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
              ["norder_index",     "",        20, [0, 31],   "",       
                  "Order index for the Lebedev quadrature, see orderlist"]
              ]


def volume(length_a, b2a_ratio, c2a_ratio, t):
    return (4/3)*length_a**3*b2a_ratio*c2a_ratio*(1.-3*(1.-t)**3)

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

    intensity=amp3D(qx,qy,qz,length_a, b2a_ratio, c2a_ratio, t)**2

    return intensity 




######################################## Lebedev integration ##################################################

def lebedev_sphere(norder_index: int):
    """
    Generates quasi-uniformly distributed points on the unit sphere
    in Cartesian coordinates (x,y,z) and their associated weights.
    Parameters
    ----------
    norder : int
        index of the order of the Lebedev quadrature (see orderlist).
    Returns
    -------
    points : ndarray, shape (npoints, 3)
        Cartesian coordinates of the points on the unit sphere.
    weights : ndarray, shape (npoints,)
        Weights associated with each point for integration on the sphere.

    """

    order = orderlist[int(norder_index)]
    q_unit, w = leblib.get_points_and_weights(order)  # shape (npoints,3)

    return q_unit, w


def plot_lebedev_sphere(norder_index, figsize=(7,7)):
    """
    3D representation of Lebedev points on the unit sphere.
    Parameters
    ----------
    norder : int
        Index in orderlist for the Lebedev order.
    figsize : tuple
        Size of the figure.
    """
    pts, w = lebedev_sphere(norder_index)  # shape (npoints,3)
    order = orderlist[int(norder_index)]
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(pts[:,0], pts[:,1], pts[:,2], s=10, alpha=0.6)

    # Axes et esthétique
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_title(f"Lebedev points on the unit sphere  (order of the polynom: {order}, {len(pts)} total points)")
    ax.set_box_aspect([1,1,1])  # sphère non déformée
    plt.show()




def Iq(q, sld, sld_solvent,length_a=500, b2a_ratio=1, c2a_ratio=1, t=0.99,  norder_index:int=30, rotate:bool=True):
    """
    Parameters:
        q:      input scattering vectors, units 1/Ang
    Returns:
        I(q):   1D scattering intensity at q, units 1/cm
    """


    time_start = time.time()
    q = np.atleast_1d(q)  # vecteur q
    q_unit, w = lebedev_sphere(norder_index)  # shape (npoints,3)

    if rotate:
        q_unit = _rotate_unit_vectors(q_unit)


    # build qx,qy,qz arrays with correct broadcasting -> shape (nq, npoints)
    qx = q[:, None] * q_unit[None, :, 0]
    qy = q[:, None] * q_unit[None, :, 1]
    qz = q[:, None] * q_unit[None, :, 2]

    # compute intensity grid using existing amp3D (vectorized)
    amp_grid = amp3D(qx, qy, qz, length_a, b2a_ratio, c2a_ratio, t)  # shape (nq, npoints)
    intensity_grid = np.abs(amp_grid)**2

    integral = np.sum(intensity_grid * w[None, :], axis=1) #summation over all points

    time_end = time.time()
    time_lebedev= time_end - time_start
    print(f'Execution time Lebedev with {len(q_unit)} points: {time_lebedev:.4f} seconds')


    return integral*0.0001*(sld-sld_solvent)**2*volume(length_a, b2a_ratio, c2a_ratio, t)


Iq.vectorized = True # Here Iq works only for a single float value of q parameter