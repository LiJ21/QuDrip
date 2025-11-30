#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <utility>

#include "find_param.hpp"
#include "qudrip.hpp"
#include "lattice.hpp"

#define TUNIF
using namespace qudrip;
using namespace std;
using cplx = complex<double>;

template <typename BIT2>
auto tlf(BIT2&& hopgate, int Nsite, int site1, int site2,
         const params_type<double>& params) {
  auto& thop = params("thop");
  auto& Jx = params("Jx");
  auto& Jz = params("Jz");
  int sm = params("sm");

  if (sm < 0) {
    auto& phi = params("phi");
    auto n1u = hopgate(site1, site1);
    auto n1d = hopgate(site1 + Nsite, site1 + Nsite);
    auto n2u = hopgate(site2, site2);
    auto n2d = hopgate(site2 + Nsite, site2 + Nsite);

    auto hop12u = hopgate(site1, site2);
    auto hop21u = hopgate(site2, site1);

    auto hop12d = hopgate(site1 + Nsite, site2 + Nsite);
    auto hop21d = hopgate(site2 + Nsite, site1 + Nsite);

    auto phase = cplx(cos(phi), sin(phi));
    auto phase2 = phase * phase;

    auto rstHop =
        -thop * (phase * (hop12d - n1u * hop12d - hop12d * n2u +
                          2. * n1u * hop12d * n2u + hop12u - n1d * hop12u -
                          hop12u * n2d + 2. * n1d * hop12u * n2d) +
                 std::conj(phase) * (hop21d - n2u * hop21d - hop21d * n1u +
                                     2. * n2u * hop21d * n1u

                                     + hop21u - n2d * hop21u - hop21u * n1d +
                                     2. * n2d * hop21u * n1d));

    auto tlInt =
        (Jx * (phase2 * hop12u * hop12d + std::conj(phase2) * hop21u * hop21d) -
         Jx * (hop12u * hop21d + hop21u * hop12d)) *
            value_type(.5) +
        ((Jx + Jz) * (n1u * n2u + n1d * n2d) +
         (-Jx + Jz) * (n1u * n2d + n1d * n2u) - Jz * (n1u + n2u + n1d + n2d)) *
            value_type(.25);

    return rstHop + tlInt;
  } else {
    auto& phi1 = params("phi1");  // s1 and sm
    auto& phi2 = params("phi2");  // sm and s2
    auto n1u = hopgate(site1, site1);
    auto n1d = hopgate(site1 + Nsite, site1 + Nsite);
    auto nmu = hopgate(sm, sm);
    auto nmd = hopgate(sm + Nsite, sm + Nsite);
    auto n2u = hopgate(site2, site2);
    auto n2d = hopgate(site2 + Nsite, site2 + Nsite);

    auto hop12u = hopgate(site1, site2);
    auto hop21u = hopgate(site2, site1);
    auto hop2mu = hopgate(site2, sm);
    auto hopm2u = hopgate(sm, site2);
    auto hop1mu = hopgate(site1, sm);
    auto hopm1u = hopgate(sm, site2);

    auto hop12d = hopgate(site1 + Nsite, site2 + Nsite);
    auto hop21d = hopgate(site2 + Nsite, site1 + Nsite);
    auto hop2md = hopgate(site2 + Nsite, sm + Nsite);
    auto hopm2d = hopgate(sm + Nsite, site2 + Nsite);
    auto hop1md = hopgate(site1 + Nsite, sm + Nsite);
    auto hopm1d = hopgate(sm + Nsite, site2 + Nsite);

    auto phase1 = cplx(cos(phi1), sin(phi1));
    auto phase2 = cplx(cos(phi2), sin(phi2));

    auto Int_h = -Jx * .25 * phase1 * phase2 *
                     ((nmd * hop12u - nmd * n2d * hop12u - nmd * n1d * hop12u +
                       nmd * n2d * n1d * hop12u) +
                      (nmu * hop12d - nmu * n2u * hop12d - nmu * n1u * hop12d +
                       nmu * n2u * n1u * hop12d)) -
                 Jx * .25 * std::conj(phase1 * phase2) *
                     ((nmd * hop21u - nmd * n2d * hop21u - nmd * n1d * hop21u +
                       nmd * n2d * n1d * hop21u) +
                      (nmu * hop21d - nmu * n2u * hop21d - nmu * n1u * hop21d +
                       nmu * n2u * n1u * hop21d)) -
                 Jx * .25 * phase1 * phase2 *
                     ((hop1md * hopm2u - n1u * hop1md * hopm2u -
                       hop1md * hopm2u * n2d + n1u * hop1md * hopm2u * n2d) +
                      (hop1mu * hopm2d - n1d * hop1mu * hopm2d -
                       hop1mu * hopm2d * n2u + n1d * hop1mu * hopm2d * n2u)) -
                 Jx * .25 * std::conj(phase1 * phase2) *
                     ((hopm1d * hop2mu - n1u * hopm1d * hop2mu -
                       hopm1d * hop2mu * n2d + n1u * hopm1d * hop2mu * n2d) +
                      (hopm1u * hop2md - n1d * hopm1u * hop2md -
                       hopm1u * hop2md * n2u + n1d * hopm1u * hop2md * n2u));

    auto Int_d =
        Jx * .25 * std::conj(phase1 * phase2) *
            ((hop21u * n2d * n1d - nmd * hop21u * n2d * n1d) +
             (hop21d * n2u * n1u - nmu * hop21d * n2u * n1u)) +
        Jx * .25 * phase1 * phase2 *
            ((hop12u * n2d * n1d - nmd * hop12u * n2d * n1d) +
             (hop12d * n2u * n1u - nmu * hop12d * n2u * n1u)) -
        Jx * .25 * std::conj(phase1 * phase2) *
            (hopm1d * n1u * n2d * hop2mu + hopm1u * n1d * n2u * hop2md) -
        Jx * .25 * phase1 * phase2 *
            (hop1md * n1u * n2d * hopm2u + hop1mu * n1d * n2u * hopm2d);

    auto Int_slide =
        Jx * .25 * std::conj(phase1) * phase2 *
            ((hopm2u * hopm1d * n1u - hopm2u * n2d * hopm1d * n1u) +
             (hopm2d * hopm1u * n1d - hopm2d * n2u * hopm1u * n1d)) +
        Jx * .25 * phase1 * std::conj(phase2) *
            ((hop2mu * hop1md * n1u - hop2mu * n2d * hop1md * n1u) +
             (hop2md * hop1mu * n1d - hop2md * n2u * hop1mu * n1d)) +
        Jx * .25 * phase1 * std::conj(phase2) *
            ((n2d * hop2mu * hop1md - n2d * hop2mu * n1u * hop1md) +
             (n2u * hop2md * hop1mu - n2u * hop2md * n1d * hop1mu)) +
        Jx * .25 * std::conj(phase1) * phase2 *
            ((n2d * hopm2u * hopm1d - n2d * hopm2u * n1u * hopm1d) +
             (n2u * hopm2d * hopm1u - n2u * hopm2d * n1d * hopm1u));

    return Int_h;
  }
}

template <typename BIT2, typename IDX>
auto tlfMatrix(BIT2&& hopgate, IDX&& idx, int Nsite, int site1, int site2,
               const params_type<double>& params) {
  return tlf(std::forward<BIT2>(hopgate), Nsite, site1, site2, params) >> idx;
}

int main(int argc, char** argv) {
  int branch = 0;
  int ave = 0;
  auto must = 1;
  int kry_dim, Nsite;
  double t0 = -1.0;
  double Jx, Jz;
  int store_correl = 0;

  int print_ham;
  double fr;
  string solver_type;

  int compute_chi = 0;
  param_finder pf(argv[1]);
  int docc = 0;
  int Niter = 1;
  int nconv = 0;
  double J1;
  double phi;

  params_type<double> param;
  param("Jx") = 1.0;
  param("Jz") = 1.0;
  param("sm") = -1;

  int Nx = 1, Ny = 1;
  int Nup = 0, Ndo = 0;

  try {
    pf.find_param("kry_dim", kry_dim, must);
    pf.find_param("Nx", Nx);
    pf.find_param("Ny", Ny);
    pf.find_param("t0", t0);
    pf.find_param("Nup", Nup);
    pf.find_param("Ndo", Ndo);
    pf.find_param("docc", docc);
    pf.find_param("Jx", Jx);
    pf.find_param("Jz", Jz);
    pf.find_param("phi", param("phi"));
    pf.find_param("solver_type", solver_type);
    pf.find_param("nconv", nconv);

    pf.find_param("Niter", Niter);
    pf.find_param("print", print_ham);
    pf.find_param("store_correl", store_correl);

    pf.find_param("branch", branch);
    pf.find_param("ave", ave);

  } catch (string& pp) {
    cout << "Missing " << pp << endl;
    abort();
  }

  vector<StaticIVector<2>> coords{}, pvecs{};
  pvecs.push_back(StaticIVector<2>{Nx, 0});
  pvecs.push_back(StaticIVector<2>{0, Ny});
  for (int i = 0; i < Nx; ++i) {
    for (int j = 0; j < Ny; ++j) {
      coords.push_back(StaticIVector<2>{i, j});
    }
  }

  Lattice<2> lat(coords, pvecs);

  Model<double, 2> m(lat);

  param("Jx") = Jx;
  param("Jz") = Jz;
  param("thop") = t0;
  param("phi") *= M_PI / 180;
  phi = param("phi");
  Nsite = lat.N();  // consider the 3*4 cluster

  auto fermi = getQbits(2 * Nsite);
  auto bitgate = qudrip::getFermiGate(fermi);
  auto hopgate = qudrip::getHopFermiGate(fermi);

  cout << fermi.range() << endl;
  cout << sizeof(fermi.range()) << endl;

  auto ntot = getNset(fermi);
  auto nup = getNup(fermi, Nsite);
  auto ndo = getNdo(fermi, Nsite);
  auto nd = getNd(fermi, Nsite);
  auto nh = getNh(fermi, Nsite);
  auto time1 = chrono::high_resolution_clock::now();
  int nup_v, ndo_v;
  nup_v = Nsite % 2 ? Nsite / 2 + 1 : Nsite / 2;
  ndo_v = Nsite / 2;
  auto hubbard = Constrain(fermi, nup = Nup, ndo = Ndo, nd = docc);
  auto time2 = chrono::high_resolution_clock::now();
  std::cout << hubbard.range() << endl;
  auto ndim = hubbard.range();
  time1 = chrono::high_resolution_clock::now();
  auto psi = getState(hubbard, 1);
  //-----------------------------------------

  if (print_ham)
    for (hubbard[0]; hubbard < hubbard.range(); ++hubbard)
      cout << bitset<8>(fermi) << endl;

  StaticIVector<2> a0{1, 0};
  StaticIVector<2> a1{0, 1};

  for (int i = 0; i < lat.N(); ++i) {
    auto& c = lat(i);
    m.addBond(i, lat(c + a0), param);
    m.addBond(i, lat(c + a1), param);
  }

  m.print();

  //...
  SpMatrix Ham(hubbard.range(), hubbard.range());
  if (branch) {
    auto Hamg =
        m.getHamiltonian([&](int s1, int s2, const params_type<double>& param) {
          return tlf(hopgate, Nsite, s1, s2, param);
        });
    Ham = Hamg >> hubbard;
    auto time1 = chrono::high_resolution_clock::now();
    Matrix tmp1 = Ham * psi(0).mat();
    auto time2 = chrono::high_resolution_clock::now();
    std::cout
        << "Matrix multi. Time elapsed: "
        << chrono::duration_cast<chrono::milliseconds>(time2 - time1).count() *
               1e-3
        << "s" << std::endl;
  } else {
    Ham = m.getHamiltonianMatrix(
        [&](int s1, int s2, const params_type<double>& param) {
          return tlfMatrix(hopgate, hubbard, Nsite, s1, s2, param);
        },
        ndim);
  }

  if (print_ham) {
    cout << Matrix(Ham) << endl;
  }
  time2 = chrono::high_resolution_clock::now();
  std::cout
      << "Time elapsed: "
      << chrono::duration_cast<chrono::milliseconds>(time2 - time1).count() *
             1e-3
      << "s" << std::endl;
  std::cout << "Done." << std::endl;

  time1 = chrono::high_resolution_clock::now();

  auto mc = getMC(1, psi);

  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> dice(-1.0, 1.0);

  for (auto i : range(ndim)) {
    hubbard[i];
    psi[0] = cplx(dice(mt), dice(mt));
  }
  psi(0) /= std::sqrt(Hermitian::norm(psi(0)));
  if (print_ham) cout << "randomized psi: " << endl << psi(0) << endl;
  cout << psi(0).mat().adjoint() * psi(0).mat() << endl;
  double Etot = psi.eval(0, Ham).real();
  value_type n = 0;

  // observables

  {
    int Nsq = Nsite * Nsite;
    mc.addObs("Etot", 1, MC_REAL);
    mc.addObs("phi_exp", 1, MC_REAL);
    mc.addObs("eta_exp", 1, MC_REAL);
    mc.addObs("tri_exp", 1, MC_REAL);

    mc.addObs("pair_correl", Nsite);
    mc.addObs("sxy_correl", Nsite);
    mc.addObs("cdw_correl", Nsite);
    mc.addObs("sz_correl", Nsite);

  }

  SpMatrix phi_mat(ndim, ndim);
  SpMatrix eta_mat(ndim, ndim);
  SpMatrix tri_mat(ndim, ndim);

  phi_mat.setZero();
  eta_mat.setZero();
  tri_mat.setZero();

  LatticeFunction<value_type, 2> sgn_eta(lat);
  for (int i = 0; i < lat.N(); ++i) {
    auto& c = lat(i);
    sgn_eta(i) = (c(0) % 2 ? 1 : -1) * (c(1) % 2 ? 1 : -1);
  }
  for (int s : range(lat.N())) {
    for (int sp : range(lat.N())) {
      SpMatrix pair_mat(ndim, ndim);
      pair_mat.setZero();

      if (s != sp) {
        int site1u = s;
        int site2u = sp;
        auto hop21u = hopgate(site2u, site1u);
        auto hop12u = hopgate(site1u, site2u);

        int site1d = s + Nsite;
        int site2d = sp + Nsite;
        auto hop21d = hopgate(site2d, site1d);
        auto hop12d = hopgate(site1d, site2d);

        pair_mat = hop21u * hop21d >> hubbard;

        {
          mc.incrMatrix("pair_correl", lat(lat(s) - lat(sp)),
                        hop21u * hop21d >> hubbard);
          mc.incrMatrix("sxy_correl", lat(lat(s) - lat(sp)),
                        -(hop21u * hop12d >> hubbard));

          auto ns = (hopgate(s, s) + hopgate(s + Nsite, s + Nsite));
          auto n0 = (hopgate(sp, sp) + hopgate(sp + Nsite, sp + Nsite));
          mc.incrMatrix("cdw_correl", lat(lat(s) - lat(sp)),
                        n0 * ns - n0 - ns >> hubbard);

          auto szs = hopgate(s, s) - hopgate(s + Nsite, s + Nsite);
          auto sz0 = (hopgate(sp, sp) - hopgate(sp + Nsite, sp + Nsite));
          mc.incrMatrix("sz_correl", lat(lat(s) - lat(sp)),
                        (szs * sz0 >> hubbard) / 4.0);
        }

      } else {
        pair_mat = hopgate(s, s) * hopgate(s + Nsite, s + Nsite) >> hubbard;

        {
          mc.incrMatrix(
              "pair_correl", 0,
              hopgate(s, s) * hopgate(s + Nsite, s + Nsite) >> hubbard);
          mc.incrMatrix(
              "sxy_correl", 0,
              hopgate(s, s) - hopgate(s, s) * hopgate(s + Nsite, s + Nsite) >>
                  hubbard);

          auto ns = (hopgate(s, s) + hopgate(s + Nsite, s + Nsite));
          auto n0 = (hopgate(s, s) + hopgate(s + Nsite, s + Nsite));
          mc.incrMatrix("cdw_correl", 0, n0 * ns - n0 - ns >> hubbard);

          auto szs = hopgate(s, s) - hopgate(s + Nsite, s + Nsite);
          auto sz0 = (hopgate(s, s) - hopgate(s + Nsite, s + Nsite));
          mc.incrMatrix("sz_correl", 0, (szs * sz0 >> hubbard) / 4.0);
        }
      }

      auto s1 = lat(s)(0);
      auto s2 = lat(s)(1);

      auto sp1 = lat(sp)(0);
      auto sp2 = lat(sp)(1);

      value_type tri_sgn;
      int ss = (s1 - sp1 - s2 + sp2) % 3;
      if (ss < 0) ss += 3;

      int sss = (s1 - sp1 - s2 + sp2) % 2;
      if (sss < 0) sss += 2;

      switch (ss) {
        case 0:
          tri_sgn = 1.;
          break;
        case 1:
          if (ave)
            tri_sgn = cos(2 * M_PI / 3.);
          else
            tri_sgn = value_type(cos(2 * M_PI / 3.), sin(2 * M_PI / 3.));
          break;
        case 2:
          if (ave)
            tri_sgn = cos(4. * M_PI / 3.);
          else
            tri_sgn = value_type(cos(4. * M_PI / 3.), sin(4. * M_PI / 3.));
      }

      phi_mat += pair_mat;
      eta_mat += sgn_eta(s) * sgn_eta(sp) * pair_mat;
      tri_mat += tri_sgn * pair_mat;
    }

  }

  for (int s = 0; s < Nsite; ++s) {
    mc.normMatrix("pair_correl", s, 1. / Nsite);
    mc.normMatrix("sxy_correl", s, 1. / Nsite);
    mc.normMatrix("cdw_correl", s, 1. / Nsite);
    mc.normMatrix("sz_correl", s, 1. / Nsite);
  }
  mc.updateMatrix("eta_exp", 0, eta_mat);
  mc.updateMatrix("tri120_exp", 0, eta_mat);
  mc.updateMatrix("tri_exp", 0, tri_mat);
  mc.updateMatrix("phi_exp", 0, phi_mat);

  std::cout << "initial energy:" << Etot << std::endl;
  double err;
  int i = 1;

  std::cout << "solver type : " << solver_type << std::endl;
  mc.updateMatrix("Etot", 0, Ham);
  if (solver_type == "full") {
    Eigen::SelfAdjointEigenSolver<Matrix> eigensolver(Ham);
    cout << "Full solution: start ground state" << endl;
    psi(0).mat() = eigensolver.eigenvectors().col(0);
    ofstream ef("e.out");
    for (auto i : range(eigensolver.eigenvalues().rows())) {
      ef << eigensolver.eigenvalues()(i, 0) << endl;
      if (print_ham)
        cout << i << "th vector: " << endl
             << eigensolver.eigenvectors().col(i) << endl;
    }
    ef.close();
  } else if (solver_type == "Arnoldi") {
    auto solver = Arnoldi::getSolver(hubbard, kry_dim);
    solver.convergeEV(10, Etot);
    while (solver.act_dim() < kry_dim) {
      cout << "Iter: " << solver.act_dim() << endl;

      cout << "Arnoldi: start to compute the ground state" << endl;
      solver.getGroundState(Ham, psi, 0, Niter);

      cout << "E = " << psi.eval(0, Ham) << endl;
      cout << "Error = " << solver.error() << endl;
      cout << "Dim = " << solver.act_dim() << endl;
      // if(solver.act_dim() > 30 && abs(solver.error()) < exact::eps)
      if (abs(solver.error()) < qudrip::eps) {
        break;
      }
    }

    solver.storeE("e.out");
  } else if (solver_type == "Lanczos") {
    auto solver = Arnoldi::getLanczosSolver(hubbard, kry_dim);
    solver.convergeEV(nconv, Etot);
    while (solver.act_dim() < kry_dim) {
      cout << "Iter: " << solver.act_dim() << endl;

      cout << "Hermitian: start to compute the ground state" << endl;
      solver.getGroundState(Ham, psi, 0, Niter);

      cout << "E = " << psi.eval(0, Ham) << endl;
      cout << "Error = " << solver.error() << endl;
      cout << "Dim = " << solver.act_dim() << endl;
      if (abs(solver.error()) < qudrip::eps) {
        break;
      }
    }

    solver.storeE("e.out");
  }

  cout << "ground state computed" << endl;
  mc.clearAll();
  mc.incrAllObs(0);

  cout << "norm = " << std::sqrt(Hermitian::norm(psi(0))) << endl;
  std::cout << "ground state energy: " << mc.obsValue("Etot", 0, 0) << endl;
  std::cout << "(check hermiticity) ground state energy: "
            << psi.eval(0, Ham.adjoint()) << endl;

  time2 = chrono::high_resolution_clock::now();
  std::cout
      << "Time elapsed: "
      << chrono::duration_cast<chrono::milliseconds>(time2 - time1).count() *
             1e-3
      << "s" << std::endl;

  if (print_ham) cout << psi(0) << endl;

#ifdef TUNIF
  cout << "docc distribution: ";
  for (auto i : range(Nsite)) {
    int site1u = i;
    int site1d = i + Nsite;
    auto nup = hopgate(site1u, site1u);
    auto ndo = hopgate(site1d, site1d);

    cout << "\t" << psi.eval(0, (nup * ndo) >> hubbard);
  }
  cout << endl;

  cout << "charge distribution: ";
  for (auto i : range(Nsite)) {
    int site1u = i;
    int site1d = i + Nsite;
    auto nup = hopgate(site1u, site1u);
    auto ndo = hopgate(site1d, site1d);

    cout << "\t" << psi.eval(0, (nup + ndo) >> hubbard);
  }
  cout << endl;

  cout << "spin distribution: ";
  for (auto i : range(Nsite)) {
    int site1u = i;
    int site1d = i + Nsite;
    auto nup = hopgate(site1u, site1u);
    auto ndo = hopgate(site1d, site1d);

    cout << "\t" << psi.eval(0, (nup - ndo) >> hubbard);
  }
  cout << endl;

#endif
  mc.normAll();
  std::cout << "Finished computation. Collecting data..." << std::endl;
  mc.collectAll();
  std::cout << "Finished collecting. Writing data to file..." << std::endl;
  mc.outFile(0, "obs.out",
             std::vector<std::string>{"Etot", "phi_exp", "eta_exp", "tri_exp"});
  mc.outFile(0, "pair_correl.out",
             std::vector<std::string>{"pair_correl", "cdw_correl", "sxy_correl",
                                      "sz_correl"});
  return 0;
}
