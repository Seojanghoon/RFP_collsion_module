# RFP Collision Module

## Overview

The RFP (Rosenbluth-Fokker-Planck) collision operator module is designed to model collision processes within a gyrokinetic framework. In this implementation, the particle distribution function *f* is represented using a Discontinuous Galerkin (DG) basis, while the Rosenbluth potentials are expressed using a Finite Element Method (FEM) basis.

This module has been developed as a supplementary component for a gyrokinetic code and is not intended to be standalone. Consequently, it cannot be compiled or executed independently of the parent codebase.

The primary objective of this repository is to provide the core algorithms and logic flow underlying the RFP collision operator. While this version is not directly compilable, it contains comprehensive documentation necessary for implementing a fully functional collision operator. 

For access to the complete gyrokinetic code, please contact the principal investigator at [seojh@kfe.re.kr](mailto:seojh@kfe.re.kr). Organizations seeking access to the full codebase must establish a collaboration agreement with the Korea Institute of Fusion Energy (KFE). This agreement will define the scope of collaboration and establish security provisions for handling sensitive information contained within the complete code.

## Repository Structure

### 1. **include/**
Contains C++ header files corresponding to the source files in the `source/` directory.

### 2. **source/**
Contains the core C++ implementation files:

- **main.cpp**  
  Main entry point and top-level class instantiation.

- **gyrokinetic.cpp**  
  Setup procedures and main simulation loop.

- **collision_dg.cpp**  
  Core collision operator implementation.

- **integration.cpp, integration_dg.cpp**  
  Phase space integration routines for collision operations.

- **basis.cpp, basis_dg.cpp**  
  Discontinuous Galerkin basis function definitions for the particle distribution function *f*.

- **flux.cpp**  
  Initialization routines for the Maxwellian equilibrium distribution function.

## Computational Workflow

The execution sequence of the collision module follows these steps:

### 1. Initialization
The `Gyrokinetic` class is instantiated in `main.cpp`.

### 2. Multi-Process Execution
The `gyrokinetic.procs_mult` method is invoked from `main.cpp`.

### 3. Collision Update Sequence
Within `gyrokinetic.procs_mult` (defined in `gyrokinetic.cpp`), the following operations are performed:

- **`collision_dg->h0g0fM_update(sp_dg_coeff, step, 0);`**  
  Executes periodic updates for the lowest-order component of the distribution function *f*.

- **`collision_dg->RK_implicit_col(sp_dg_coeff, step);`**  
  Performs the primary collision operation.

### 4. Core Collision Algorithm
Within `RK_implicit_col` (defined in `collision_dg.cpp`):

- Parameter validation is performed for density (*n*), flow velocity (*U*), and temperature (*T*).
- **`RK_implicit_col_single`** is executed, containing the complete collision algorithm implementation. All necessary operations are documented in detail within this subroutine.

## Contact Information

For questions regarding implementation details or requests for collaboration, please contact the principal investigator at [seojh@kfe.re.kr](mailto:seojh@kfe.re.kr).
