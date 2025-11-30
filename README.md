QuDrip is a modern C++ library for easy implementation of exact quantum many-body simulation.

It includes roughly three parts: 
  1) the core library (operator, index, conservation...) which helps create Hamiltonian from (local or not) operators and restrict the Hilbert space by imposing conservation law.
  2) Sparse matrix algorithms (Krylov-space based, Lanczos and Arnoldi algorithms...) for ground-state search (energy optimization) and real-time evolution.
  3) Various utilities for, e.g., (Monte Carlo) super jump method for (Markovian) Lindblad dynamics simulation, and lattice construction...

The library requires Eigen3.
