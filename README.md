1. Generate the mesh  
   Run `meshRectangle.m` to create the nodes and edges for the computational domain.  
   - This script generates the following files:  
     `pS_2D.m`, `pC_2D.m`, `tS_2D.m`, `tC_2D.m`,  
     `bottom_pC_1D.m`, `bottom_pS_1D.m`, `bottom_tC_1D.m`, `bottom_tS_1D.m`,  
     `top_pS_1D.m`, `top_pC_1D.m`, `top_tC_1D.m`, `top_tS_1D.m`.  
   - These files contain the edges and nodes for the entire domain, as well as the top and bottom boundary edges of the rectangular domain.  
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
