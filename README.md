# QuDrip

**QuDrip** is a modern **C++** library designed for easy and efficient implementation of **exact quantum many-body simulations**. It provides a flexible operator framework, Hilbert-space restrictions via conservation laws, sparse-matrix algorithms, and utilities for Hamiltonian and Lindblad-type dynamics.

---

## Features

### 1. Core Library
- Operator framework (local and non-local)
- Index and Hilbert-space representations
- Conservation-law machinery
- Tools for constructing Hamiltonians from gates/operators

### 2. Sparse-Matrix Algorithms
- Krylov-space methods
- **Lanczos** algorithm for ground-state search
- **Arnoldi** algorithm for dynamics and spectra
- Real-time evolution utilities

### 3. Utilities
- Monte Carlo **super-jump method** for Markovian Lindblad dynamics
- Lattice and geometry constructors
- Helper routines for simulation workflows

---

## Requirements

- **Eigen3**

---

## Basic Workflow

A typical workflow in QuDrip:

1. Create a Hilbert space (e.g., qubits or fermions/bosons) and impose a conserved quantity.
2. Define local or non-local operators (“gates”).
3. Construct a Hamiltonian either lazily (computational graph) or eagerly (explicit sparse matrices).
4. Use Krylov-based solvers such as Lanczos to obtain ground states or real-time evolution.

---

## Minimal Example (XY Spin Model)

Create qubit chain (the Hilbert space), here it is a space of $2^Nsite$ dimension organized as a one-dimensional spin chain.
```cpp
auto prechain = getQbits(Nsite);
```
Corresponding conserved quantity
```cpp
auto Nset = getNset(prechain);
```
Create boson (single-qbit) gates
```cpp
auto gc = getBoseGate(prechain);
auto ga = getBoseGate(prechain);
```

Constrain Hilbert space: Nsite/2 spin-up and spin-down
```cpp
auto chain = Constrain(prechain, Nset = Nsite / 2);
```

Feed spin matrices into gates. Here it is simple spin ladder operator
```cpp
Matrix Uc = 0.5 * (pauli_x + II * pauli_y);
Matrix Ua = 0.5 * (pauli_x - II * pauli_y);
gc << Uc;
ga << Ua;
```

Build Hamiltonian (the lazy way, where a computational graph is grenerated)
```cpp
auto ham = getEmptyOperator();
for (int i = 0; i < Nsite - 1; ++i) {
    ham = ham + gc(i + 1) * ga(i) + gc(i) * ga(i + 1);
}
```

Apply graph to Hilbert space → sparse matrix
```cpp
auto ham_mat = ham >> chain;
cout << "Time elapsed: " << c1.release() << std::endl;
```

Lanczos ground-state search
```cpp
int krylov_dim = 100;
auto psi = getState(chain, 1);
shuffle_state(psi(0));
auto solver = Arnoldi::getLanczosSolver(chain, krylov_dim);

solver.getGroundState(ham_mat, psi(0), 0);
```

Ground-state energy is then evaluated
```cpp
cout << "E = " << psi.eval(0, ham_mat) << endl;
```
