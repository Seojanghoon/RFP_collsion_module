# RFP Collision Module

## Overview

The RFP (Rosenbluth Fokker-Planck) collision operator module is designed to model collision processes in a gyrokinetic context. In this module, the particle distribution function ( f ) is represented using a Discontinuous Galerkin (DG) basis, while the Rosenbluth potentials are expressed using a Finite Element Method (FEM) basis. This module has been developed as a supplementary component for a gyrokinetic code. As such, it is not standalone and cannot be compiled or executed independently.

The primary goal of this repository is to provide the core algorithm and flow of logic underpinning the RFP collision operator. Although this version is not directly compilable, it includes detailed information necessary to implement a fully functional collision operator. For access to the complete gyrokinetic code, please contact the primary author at [seojh@kfe.re.kr](mailto:seojh@kfe.re.kr).

Should the requesting organization wish to gain access to the full code, a collaboration agreement with the Korea Institute of Fusion Energy (KFE) is required. This agreement will outline the scope of collaboration and security provisions related to sensitive information contained in the full code.

## Overall Workflow

The following steps outline the key processes involved in the execution of the collision module:

1. **Initialize Gyrokinetic Class in `main.cpp`:**
   Instantiate the `Gyrokinetic` class in `main.cpp`.

2. **Call `gyrokinetic.procs_mult`:**
   In `main.cpp`, invoke the `gyrokinetic.procs_mult` method.

3. **In `gyrokinetic.procs_mult` (located in `gyrokinetic.cpp`):**

   * `collision_dg->h0g0fM_update(sp_dg_coeff, step, 0);`
     Performs a periodic update for the lowest-order part of the distribution function ( f ).

   * `collision_dg->RK_implicit_col(sp_dg_coeff, step);`
     Executes the main collision operation.

4. **In `RK_implicit_col` (located in `collision_dg.cpp`):**

   * Validates the collision operator parameters ( n ), ( U ), and ( T ).

   * `RK_implicit_col_single`:
     The core collision algorithm is performed here. All operations necessary for this step are described in detail within this module.

For further clarification or questions about the implementation, please contact the main author at [seojh@kfe.re.kr](mailto:seojh@kfe.re.kr).
