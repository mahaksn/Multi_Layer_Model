1. Generate the mesh  
   Run `meshRectangle.m` to create the nodes and edges for the computational domain.  
   - You can adjust the *domain size* and *mesh size* within this script.

2. Run the main simulation  
   Execute `main_2L_y2D.m` to perform the 2D simulations with transport coupling.  
   - This script calls the following main functions:
     - `AssembleGlobalMatrices_corrected2.m` — assembles the FEM matrices using quadrature points.  
     - `ReactKineInt.m` — computes the reaction kinetics.  
     - `NLBoundFluxIn.m` — applies nonlinear boundary flux terms.
   - These functions use the following helper files for basis functions and quadrature points:
     - `basis_linear_1D.m`
     - `basis_linear_2D.m`
     - `RefEdgeQuad.m`
     - `RefElemQuad.m`
