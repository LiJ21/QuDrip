// Regression tests for the QuDrip core machinery.
// Built as two translation units (see test_second_tu.cpp), so a successful
// link also guards against non-inline definitions in the headers.
#include "qudrip.hpp"
#include "lattice.hpp"

#include <unsupported/Eigen/KroneckerProduct>

#include <cstdlib>
#include <iostream>

using namespace qudrip;

static int n_checks = 0;
#define CHECK(cond)                                                       \
  do {                                                                    \
    ++n_checks;                                                           \
    if (!(cond)) {                                                        \
      std::cerr << "FAILED (line " << __LINE__ << "): " #cond << "\n";    \
      std::exit(1);                                                       \
    }                                                                     \
  } while (0)

// defined in test_second_tu.cpp
Matrix other_tu_gate_matrix();

static Matrix dense(const SpMatrix& m) { return Matrix(m); }

int main() {
  // --- a gate matrix must be applied exactly as fed, not transposed ---
  {
    auto chain = getQbits(1);
    auto g = getBoseGate(chain);
    Matrix U(2, 2);
    U << value_type(1, 2), value_type(3, -1),
         value_type(0, 5), value_type(-2, 0);
    g << U;
    CHECK(dense(g(0) >> chain).isApprox(U));

    auto g2 = getBoseGate(chain);
    g2 << pauli_y;
    CHECK(dense(g2(0) >> chain).isApprox(pauli_y));
  }

  // --- site placement: site 0 is the low bit of the index ---
  {
    auto chain = getQbits(2);
    auto g = getBoseGate(chain);
    g << pauli_y;
    Matrix I2 = Matrix::Identity(2, 2);
    CHECK(dense(g(0) >> chain).isApprox(Eigen::kroneckerProduct(I2, pauli_y).eval()));
    CHECK(dense(g(1) >> chain).isApprox(Eigen::kroneckerProduct(pauli_y, I2).eval()));
  }

  // --- fermionic gates: canonical anticommutation via the JW string ---
  {
    auto chain = getQbits(3);
    auto c = getFermiGate(chain);
    auto cd = getFermiGate(chain);
    Matrix a(2, 2), adag(2, 2);
    a << 0, 1, 0, 0;     // |1> -> |0>
    adag << 0, 0, 1, 0;  // |0> -> |1>
    c << a;
    cd << adag;
    Matrix Id8 = Matrix::Identity(8, 8);
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) {
        Matrix ci = dense(c(i) >> chain);
        Matrix cdj = dense(cd(j) >> chain);
        Matrix anti = ci * cdj + cdj * ci;
        if (i == j)
          CHECK(anti.isApprox(Id8));
        else
          CHECK(anti.norm() < 1e-12);
      }
  }

  // --- constrained XY chain: graph expansion, Lanczos, and lazy apply ---
  {
    const int N = 6;
    auto prechain = getQbits(N);
    auto Ns = getNset(prechain);
    auto gc = getBoseGate(prechain);
    auto ga = getBoseGate(prechain);
    auto chain = Constrain(prechain, Ns = N / 2);
    CHECK(chain.range() == 20);  // C(6,3)

    Matrix Uc = 0.5 * (pauli_x + II * pauli_y);
    Matrix Ua = 0.5 * (pauli_x - II * pauli_y);
    gc << Uc;
    ga << Ua;

    auto ham = getEmptyOperator();
    for (int i = 0; i < N - 1; ++i)
      ham = ham + gc(i + 1) * ga(i) + gc(i) * ga(i + 1);

    SpMatrix hmat = ham >> chain;
    Matrix hd = dense(hmat);
    CHECK(hd.isApprox(hd.adjoint().eval()));

    Eigen::SelfAdjointEigenSolver<Matrix> es(hd);
    double e0 = es.eigenvalues()(0);

    auto psi = getState(chain, 2);
    shuffle_state(psi(0));
    auto solver = Arnoldi::getLanczosSolver(chain, 19);
    solver.getGroundState(hmat, psi, 0);
    CHECK(std::abs(psi.eval(0, hmat).real() - e0) < 1e-8);

    // lazy apply: psi2(1) = ham * psi(0) must equal hmat * psi(0)
    // and must leave the source untouched
    auto psi2 = getState(chain, 2);
    Matrix before = psi.data().col(0);
    psi2(1) = ham * psi(0);
    Matrix expect = hmat * before;
    CHECK(psi2.data().col(1).isApprox(expect, 1e-8));
    CHECK(psi.data().col(0).isApprox(before));
  }

  // --- restrict(): keep in-subspace elements, drop the rest ---
  {
    auto pre = getQbits(2);
    auto Ns = getNset(pre);
    auto sub = Constrain(pre, Ns = 1);  // {|01>, |10>} -> dim 2
    CHECK(sub.range() == 2);

    SpMatrix full(4, 4);
    full.setIdentity();
    SpMatrix r = restrict(pre, sub, full);
    CHECK(r.rows() == 2 && r.cols() == 2);
    CHECK(dense(r).isApprox(Matrix::Identity(2, 2)));

    triplets_type ts;
    ts.emplace_back(0, 3, value_type(7, 0));  // |00><11|: outside the subspace
    ts.emplace_back(1, 2, value_type(0, 3));  // |01><10|: inside
    SpMatrix m2(4, 4);
    m2.setFromTriplets(ts.begin(), ts.end());
    SpMatrix r2 = restrict(pre, sub, m2);
    Matrix expect2 = Matrix::Zero(2, 2);
    expect2(0, 1) = value_type(0, 3);
    CHECK(dense(r2).isApprox(expect2));
  }

  // --- evolve: all three solver paths against the dense propagator ---
  {
    const int N = 6;
    const double h = 0.05;
    auto prechain = getQbits(N);
    auto Ns = getNset(prechain);
    auto gc = getBoseGate(prechain);
    auto ga = getBoseGate(prechain);
    auto chain = Constrain(prechain, Ns = N / 2);
    gc << Matrix(0.5 * (pauli_x + II * pauli_y));
    ga << Matrix(0.5 * (pauli_x - II * pauli_y));

    auto ham = getEmptyOperator();
    for (int i = 0; i < N - 1; ++i)
      ham = ham + gc(i + 1) * ga(i) + gc(i) * ga(i + 1);
    SpMatrix hmat = ham >> chain;
    Matrix prop = herm_exp(dense(hmat), -II * h);  // exp(-i H h), dense

    auto psi = getState(chain, 2);
    shuffle_state(psi(0));  // random normalized start
    Matrix v0 = psi.data().col(0);
    Matrix expect = prop * v0;

    // general Arnoldi (herm = 0)
    auto asolver = Arnoldi::getSolver(chain, 19);
    asolver.evolve(hmat, psi, 0, h);
    CHECK(psi.data().col(1).isApprox(expect, 1e-8));

    // General Arnoldi evolution must preserve the input amplitude.
    psi.data().col(0) = 2.0 * v0;
    psi.data().col(1).setZero();
    asolver.evolve(hmat, psi, 0, h);
    CHECK(psi.data().col(1).isApprox(Matrix(2.0 * expect), 1e-8));

    // three-vector Lanczos-mode Arnoldi (herm = 1)
    psi.data().col(0) = v0;
    psi.data().col(1).setZero();
    auto lsolver = Arnoldi::getLanczosSolver(chain, 19);
    lsolver.evolve(hmat, psi, 0, h);
    CHECK(psi.data().col(1).isApprox(expect, 1e-8));

    // Three-vector Lanczos-mode Arnoldi must preserve it as well.
    psi.data().col(0) = 2.0 * v0;
    psi.data().col(1).setZero();
    lsolver.evolve(hmat, psi, 0, h);
    CHECK(psi.data().col(1).isApprox(Matrix(2.0 * expect), 1e-8));

    // Hermitian LanczosSolver
    psi.data().col(0) = v0;
    psi.data().col(1).setZero();
    auto hsolver = Hermitian::getSolver(chain, 19);
    hsolver.evolve(hmat, psi, 0, h);
    CHECK(psi.data().col(1).isApprox(expect, 1e-8));

    // Hermitian evolve is linear: a scaled input scales the output
    psi.data().col(0) = 2.0 * v0;
    psi.data().col(1).setZero();
    hsolver.evolve(hmat, psi, 0, h);
    CHECK(psi.data().col(1).isApprox(Matrix(2.0 * expect), 1e-8));
  }

  // --- McManager::obsValue must respect the timestep ---
  {
    auto idx = getQbits(1);
    auto psi = getState(idx, 3);
    psi.data()(0, 0) = 1.0;  // |0> at tstp 0
    psi.data()(1, 1) = 2.0;  // 2|1> at tstp 1
    auto mc = getMC(3, psi);
    mc.addObs("E", 2);
    SpMatrix I2(2, 2);
    I2.setIdentity();
    mc.updateMatrix("E", 0, I2);
    mc.updateMatrix("E", 1, SpMatrix(I2 * value_type(2.0)));
    mc.incrAllObs(0, 0);  // data("E")(0, :) = (1, 2)
    mc.incrAllObs(1, 1);  // data("E")(1, :) = (4, 8)
    CHECK(std::abs(mc.obsValue("E", 0, 1) - value_type(2.0)) < 1e-12);
    CHECK(std::abs(mc.obsValue("E", 1, 0) - value_type(4.0)) < 1e-12);
    CHECK(std::abs(mc.obsValue("E", 1, 1) - value_type(8.0)) < 1e-12);
  }

  // --- lattice: periodic reduction, open lattice, and both DFT overloads ---
  {
    using Latt1 = Lattice<1>;
    auto C = [](int x) {
      Latt1::coord_type c;
      c(0) = x;
      return c;
    };

    Latt1 latt(Latt1::coords_type{C(0), C(1), C(2), C(3)},
               Latt1::coords_type{C(4)});
    CHECK(latt.vol() == 16);
    CHECK(latt(C(5)) == 1);  // periodic reduction 5 -> 1

    // single-argument DFT: delta at x = 2, k = 1
    // phase = inner(x, k * rpvec) / vol = 2 * 4 / 16 = 0.5
    LatticeFunction<value_type, 1> f(latt);
    f(2) = 1.0;
    value_type expect = value_type(std::cos(0.5), std::sin(0.5)) / 4.0;
    CHECK(std::abs(f.DFT(C(1)) - expect) < 1e-12);

    // two-argument DFT (previously OOB): delta at (x1, x2) = (1, 0)
    Latt1 latt2(Latt1::coords_type{C(0), C(1)}, Latt1::coords_type{C(2)});
    LatticeFunction<value_type, 1> g(latt2, 2);
    g(std::vector<int>{1, 0}) = 1.0;
    std::vector<Latt1::coord_type> ks{C(1), C(1)};
    value_type expect2 = value_type(std::cos(0.5), std::sin(0.5)) / 4.0;
    CHECK(std::abs(g.DFT(ks) - expect2) < 1e-12);

    // lattice without periodicity vectors (previously UB in the ctor)
    Latt1 open(Latt1::coords_type{C(0), C(1), C(2)});
    CHECK(open.vol() == 0);
    CHECK(open(C(2)) == 2);
    CHECK(open(C(5)) == open.N());  // outside: index() returns N
    bool dft_rejected = false;
    try {
      LatticeFunction<value_type, 1> open_f(open);
      (void)open_f.DFT(C(0));
    } catch (const std::logic_error&) {
      dft_rejected = true;
    }
    CHECK(dft_rejected);
  }

  // --- element access through a const State must not recurse ---
  {
    auto idx = getQbits(2);
    auto psi = getState(idx, 1);
    idx[2];
    psi[0] = value_type(0.5, -0.25);

    const State<QbitsIndex<>>& cref = psi;
    idx[2];
    CHECK(cref[0] == value_type(0.5, -0.25));
    idx[4];  // out of range -> null position
    CHECK(cref[0] == value_type(0.0, 0.0));
  }

  // --- header usable from a second translation unit ---
  CHECK(other_tu_gate_matrix().isApprox(pauli_x));

  std::cout << "all " << n_checks << " checks passed" << std::endl;
  return 0;
}
