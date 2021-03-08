# ASyT
ASyT is a symplectic tracking code for accelerators that use CSABAm, TEAPOTm and SIMPOLEm integrators.

## For compilation use:
gfortran -O3 symplcetic_tracking.f08 -o abc
 
## In lxplus use:
lxplus7
## or:
scl enable devtoolset-7 bash

## References:
K. Skoufaris, Ch. Skokos, Y. Papaphilippou and J. Laskar, "Accelerator tracking via high order symplectic integration method with only foreword integration steps".
J. Laskar and P. Robutel, "High order symplectic integrators for perturbed hamiltonian systems".
H. Burkhardt, R. De Maria, M. Giovannozzi, and T. Risselada, "Improved Teapot method and tracking with thick quadrupoles for the LHC and its upgrade".
