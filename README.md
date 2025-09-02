#### SASView Form Factor models for octahedron shape including truncation
Models are computed from analytical expression of the Form factor.
* Octa_general_Lebedev.ipynb contains usage exemple of the model 
* octa-general.py is a pure python model where orientation average is computed using Lebedev quadrature
* Lebedev quadrature is implemented using the open source library:
  https://pypi.org/project/pylebedev/
* folder ./models_with_c_code/ contains other SASView models with ancillary c code where orientation average is computed using Gauss-Legendre quadrature. 
