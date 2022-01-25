# octahedron model
# Note: model title and parameter table are inserted automatically
r"""

This model provides the form factor, $P(q)$, for an octahedron.
This model is constructed in a similar way as the rectangular prism model.

Definition
----------

The octahedron is defined by three dimensions along the two-fold axis which contain the 6 vertices.
length_a, length_b and length_c are the distances from the center of the octahedron to its vertices.
Coordinates of the six vertices are :
(length_a,0,0)
(-length_a,0,0)
(0,length_b,0)
(0,-length_b,0)
(0,0,length_c)
(0,0,-length_c)

.. math::

The model is using length_a and the two ratios b2a_ratio and c2a_ratio :
 b2a_ratio = length_b/length_a
 c2a_ratio = length_c/length_a

Volume of the octahedron is:
V = (4/3) * length_a * (length_a*b2a_ratio) * (length_a*c2a_ratio)

Lengths of edges are equal to :
    A_edge^2 = length_a^2+length_b^2
    B_edge^2 = length_a^2+length_c^2
    C_edge^2 = length_b^2+length_c^2

For a regular octahedron : 
    b2a_ratio = c2a_ratio = 1   
    A_edge = B_edge = C_edge = length_a*sqrt(2)
    length_a = length_b = length_c = A_edge/sqrt(2)
    V = (4/3) * length_a^3 = (sqrt(2)/3) * A_edge^3  
    
Amplitude of the form factor AP is calculated with a scaled scattering vector (qx,qy,qz) :
    AP = (3./4.)*(A+B)/(qx*qx-qy*qy)
    with :
        A = 8.*(qy*sin(qy)-qz*sin(qz))/(qy*qy-qz*qz)
        B = 8.*(qz*sin(qz)-qx*sin(qx))/(qx*qx-qz*qz)
    and : 
            Qx = q * sin_theta * cos_phi;
    	    Qy = q * sin_theta * sin_phi;
    	    Qz = q * cos_theta;
    	    qx = Qx * length_a;
    	    qy = Qy * length_b;
    	    qz = Qz * length_c; 
       
$\theta$ is the angle between the $z$ axis and the
c axis of the octahedron ($length_c$), and $\phi$ is the angle between the
scattering vector (lying in the $xy$ plane) and the $y$ axis.

The normalized form factor in 1D is obtained averaging over all possible
orientations. Same code as for a rectangular prism.

.. math::
    P(q) =  \frac{2}{\pi} \int_0^{\frac{\pi}{2}} \,
    \int_0^{\frac{\pi}{2}} A_P^2(q) \, \sin\theta \, d\theta \, d\phi

And the 1D scattering intensity is calculated as

.. math::
    I(q) = \text{scale} \times V \times (\rho_\text{p} -
    \rho_\text{solvent})^2 \times P(q)

where $V$ is the volume of the octahedron, $\rho_\text{p}$
is the scattering length inside the volume, $\rho_\text{solvent}$
is the scattering length of the solvent, and (if the data are in absolute
units) *scale* represents the volume fraction (which is unitless).

For 2d data the orientation of the particle is required, described using
angles $\theta$, $\phi$ and $\Psi$ as in the diagrams below, for further details
of the calculation and angular dispersions see :ref:`orientation` .
The angle $\Psi$ is the rotational angle around the long *C* axis. For example,
$\Psi = 0$ when the *B* axis is parallel to the *x*-axis of the detector.

For 2d, constraints must be applied during fitting to ensure that the inequality
$A < B < C$ is not violated, and hence the correct definition of angles is preserved. The calculation will not report an error,
but the results may be not correct.

.. figure:: img/parallelepiped_angle_definition.png

    Definition of the angles for oriented core-shell parallelepipeds.
    Note that rotation $\theta$, initially in the $xz$ plane, is carried out first, then
    rotation $\phi$ about the $z$ axis, finally rotation $\Psi$ is now around the axis of the cylinder.
    The neutron or X-ray beam is along the $z$ axis.

.. figure:: img/parallelepiped_angle_projection.png

    Examples of the angles for oriented rectangular prisms against the
    detector plane.



Validation
----------

Validation of the code is made using numerical checks.


References
----------
1. Wei-Ren Chen et al. “Scattering functions of Platonic solids”. 
In:Journal of AppliedCrystallography - J APPL CRYST44 (June 2011).
DOI:10.1107/S0021889811011691
2. Croset, Bernard, "Form factor of any polyhedron: a general compact
formula and its singularities" In: J. Appl. Cryst. (2017). 50, 1245–1255 
https://doi.org/10.1107/S1600576717010147


Authorship and Verification
----------------------------

* **Authors: Marianne Imperor-Clerc (marianne.imperor@universite-paris-saclay.fr)
             Alexandra Beikert (abeikert@gmx.de)**
* **Last Modified by MI: 9 October 2020**
* **Last Reviewed by:**
"""

import numpy as np
from numpy import inf

name = "octahedron_tz"
title = "Truncated Octahedron."
description = """
    I(q)= 64(abc)^2/(Qx^2-Qy^2)(1/(1-(1-t)^3)^2)(AA+BB)^2
    P.S. :  The expression of AA and BB are detailed in c code!!
"""
category = "shape:parallelepiped"

#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["sld", "1e-6/Ang^2", 128, [-inf, inf], "sld",
               "Octahedron scattering length density"],
              ["sld_solvent", "1e-6/Ang^2", 9.4, [-inf, inf], "sld",
               "Solvent scattering length density"],
              ["length_a", "Ang", 35, [0, inf], "volume",
               "half height along a axis"],
              ["b2a_ratio", "", 1, [0, inf], "volume",
               "Ratio b/a"],
              ["c2a_ratio", "", 1, [0, inf], "volume",
               "Ratio c/a"],
              ["t", "", 0.99, [0.0, 1.0], "volume",
               "truncation along z axis"],
              ["theta", "degrees", 0, [-360, 360], "orientation",
               "c axis to beam angle"],
              ["phi", "degrees", 0, [-360, 360], "orientation",
               "rotation about beam"],
              ["psi", "degrees", 0, [-360, 360], "orientation",
               "rotation about c axis"],
             ]

source = ["lib/gauss76.c", "octahedron_tz.c"]
have_Fq = True

