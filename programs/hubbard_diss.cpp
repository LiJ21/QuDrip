#include <chrono>
#include <cmath>
#include <iostream>
#include <random>
#include <string>

#include "find_param.hpp"
#include "qudrip.hpp"

#define TUNIF
using namespace qudrip;
using namespace std;
using cplx = complex<double>;

template <typename TG, typename TI>
SpMatrix lindblad(int i, int Nsite, TG& hopgate, TI& index) {
  return (hopgate(i, i) + hopgate(i + Nsite, i + Nsite) >> index);
}

int bsearch(std::vector<double>& values, double val);
void QJ_operator(SpMatrix& H, const SpMatrix& qj_tot, const SpMatrix& Id,
                 double wt);

template <typename T, typename Tg>
vector<double> lindblad_weight(int tstp, int Nsite, const State<T>& psi,
                               Tg& hopgate, T& hubbard);

template <typename GT, typename PT, typename IT>
void incr_correl(int Nsite, int i, const PT& psi, IT& cavity, GT& gate,
                 Matrix& pair_correl, Matrix& sxy_correl, Matrix& cdw_correl,
                 Matrix& sz_correl);

int main(int argc, char** argv) {
  auto must = 1;
  int kry_dim, Nsite, kry_dim_rt;
  int Nmc = 1;
  double U0, thop = -1.0;
  double h = 0.01;
  int store_correl = 0;
  int nt = 1;
  int compute_chi = 0;
  int OPEN = 0;  // periodic
  param_finder pf(argv[1]);

  double gamma = 0.0, nb = 0;
  int dh_quench = 0;
  int dh_sep = 1;
  int Nup = 0, Ndo = 0;
  string solver_type;

  try {
    pf.find_param("solver_type", solver_type);
    pf.find_param("kry_dim", kry_dim, must);
    pf.find_param("kry_dim_rt", kry_dim_rt, must);
    pf.find_param("Nsite", Nsite, must);
    pf.find_param("Nup", Nup, must);
    pf.find_param("Ndo", Ndo, must);

    pf.find_param("h", h, must);
    pf.find_param("U", U0, must);
    pf.find_param("nt", nt, must);

    pf.find_param("thop", thop);
    pf.find_param("store_correl", store_correl);
    pf.find_param("compute_chi", compute_chi);

    pf.find_param("boundary", OPEN);

    pf.find_param("gamma", gamma);
    pf.find_param("dh_quench", dh_quench);
    pf.find_param("dh_sep", dh_sep);
    pf.find_param("Nmc", Nmc);

  } catch (string& pp) {
    cout << "Missing " << pp << endl;
    abort();
  }

  vector<cplx> Tdr(nt + 1, 0);
  vector<double> Adr(nt + 1, 0);
  vector<double> hopdr(nt + 1, 1.);
  vector<double> hop2dr(nt + 1, 0);
  vector<double> ndr(nt + 1, 0);
  try {
    pf.find_param_tvec("Adr", Adr);
    pf.find_param_tvec("hopdr", hopdr);
    pf.find_param_tvec("hop2dr", hop2dr);
    pf.find_param_tvec("ndr", ndr);
  } catch (string& pp) {
    cout << "Missing " << pp << endl;
  }

  for (int i = 0; i < nt + 1; ++i) {
    Tdr[i] = cplx(cos(Adr[i]), sin(Adr[i])) * cplx(hopdr[i], hop2dr[i]);
  }

  auto fermi = getQbits(2 * Nsite);

  auto nfermi = getNset(fermi);
  auto nup = getNup(fermi, Nsite);
  auto ndo = getNdo(fermi, Nsite);
  auto hubbard = Constrain(fermi, nup = Nup, ndo = Ndo);
  std::cout << hubbard.range() << endl;
  auto gate = getHopFermiGate(fermi);
  auto ndim = hubbard.range();
  SpMatrix Tp(ndim, ndim);
  SpMatrix Tm(ndim, ndim);
  SpMatrix Htot(ndim, ndim);
  SpMatrix h0(ndim, ndim);
  SpMatrix Id(ndim, ndim);
  SpMatrix qj_tot(ndim, ndim);
  h0.setZero();
  Tp.setZero();
  Tm.setZero();
  Htot.setZero();
  qj_tot.setZero();
  Id.setIdentity();
  auto t1 = chrono::high_resolution_clock::now();

  SpMatrix docc(ndim, ndim);
  docc.setZero();

  SpMatrix Pc(ndim, ndim);
  Pc.setIdentity();
  for (auto site : range(Nsite))  // periodic boundary condition
  {
    int site1u = site;
    int site2u = (site + 1) % Nsite;
    auto hop21u = gate(site2u, site1u);
    auto hop12u = gate(site1u, site2u);

    if (!(site == Nsite - 1 && OPEN))  // open boundary condition for 2 sites
    {
      Tp += thop * hop21u >> hubbard;
      Tm += thop * hop12u >> hubbard;
    }

    int site1d = site + Nsite;
    int site2d = (site + 1) % Nsite + Nsite;
    auto hop21d = gate(site2d, site1d);
    auto hop12d = gate(site1d, site2d);

    if (!(site == Nsite - 1 && OPEN))  // open boundary condition for 2 sites
    {
      Tp += thop * hop21d >> hubbard;
      Tm += thop * hop12d >> hubbard;
    }

    auto nup = gate(site1u, site1u);
    auto ndo = gate(site1d, site1d);

    docc += nup * ndo >> hubbard;

    h0 += (nup + ndo) >> hubbard;
    auto ni = lindblad(site, Nsite, gate, hubbard);
    qj_tot += gamma * SpMatrix(ni.adjoint()) * ni;

    Pc = Pc * (Id - (nup + ndo >> hubbard) + 2. * (nup * ndo >> hubbard));
  }

#ifdef DEBUG_H
  for (auto i : range(cavity.range())) {
    hubbard[i];
    cout << hubbard << ": " << bits_type(fermi) << endl;
  }
  cout << Matrix(H) << endl;
#endif

  auto t2 = chrono::high_resolution_clock::now();
  std::cout << "Time elapsed: "
            << chrono::duration_cast<chrono::milliseconds>(t2 - t1).count() *
                   1e-3
            << "s" << std::endl;
  std::cout << "Done." << std::endl;

  t1 = chrono::high_resolution_clock::now();
  auto psi = getState(hubbard, 2);
#ifdef USE_MPI
  auto mc = getMC(nt + 1, psi, &argc, &argv);
#else
  auto mc = getMC(nt + 1, psi);
#endif

  cout << psi(0).mat().rows() << "," << psi(0).mat().cols() << ";" << endl;

  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> dice(0.0, 1.0);

  for (auto i : range(ndim)) {
    hubbard[i];
    psi[0] = dice(mt);
  }
  psi(0).mat() /= std::sqrt(Hermitian::norm(psi(0).mat()));
  cout << psi(0).mat().adjoint() * psi(0).mat() << endl;
  double Etot = psi.eval(0, Htot).real();
  value_type n = 0;

  // observables

  {
    mc.addObs("Etot", 1, MC_REAL);
    mc.addObs("d_exp", 1, MC_REAL);
    mc.addObs("J_exp", 1, MC_REAL);
    mc.addObs("phi_exp", 1, MC_REAL);
    mc.addObs("eta_exp", 1, MC_REAL);

    mc.addObs("pair_correl", Nsite);
    mc.addObs("sxy_correl", Nsite);
    mc.addObs("cdw_correl", Nsite);
    mc.addObs("sz_correl", Nsite);

    mc.updateMatrix("d_exp", 0, docc);
  }

  // construct RVB states
  vector<cplx> rvbt1(nt + 1), rvbs1(nt + 1), rvbt2(nt + 1), rvbs2(nt + 1);
  vector<double> pcmod(nt + 1);
  auto rvb_t1 = getState(hubbard, 1);
  auto rvb_s1 = getState(hubbard, 1);
  auto rvb_t2 = getState(hubbard, 1);
  auto rvb_s2 = getState(hubbard, 1);
  string r1_idx = "";
  string r2_idx = "";
  for (int i = 0; i < Nsite; ++i) {
    r1_idx += i % 2 == 0 ? "1" : "0";  // u, d, u, d,...
    r2_idx += i % 2 == 0 ? "0" : "1";  // d, u, d, u,...
  }

  {
    string strtmp = r2_idx;
    r2_idx += r1_idx;
    r1_idx += strtmp;
  }
  cout << r1_idx << "\t" << r2_idx << endl;
  fermi = r1_idx;
  cout << hubbard << endl;
  rvb_t1[0] = 1.0;
  rvb_s1[0] = 1.0;
  fermi = r2_idx;
  cout << hubbard << endl;
  rvb_t2[0] = 1.0;
  rvb_s2[0] = 1.0;

  cout << "constructing rvb" << endl;

  SpMatrix rvb_s_mat(ndim, ndim), rvb_t_mat(ndim, ndim);
  rvb_s_mat.setIdentity();
  rvb_t_mat.setIdentity();

  SpMatrix quench_mat(ndim, ndim);

  for (int s = 0; s < Nsite; s += 2) {
    auto hop12u = gate((s + 1) % Nsite, s);
    auto hop21d = gate(s + Nsite, (s + 1) % Nsite + Nsite);
    if (!OPEN || s != Nsite - 1) {
      rvb_s_mat = hop12u - hop21d >> hubbard;
      rvb_t_mat = hop12u + hop21d >> hubbard;

      rvb_s1(0).mat() = rvb_s_mat * rvb_s1(0).mat();
      rvb_t1(0).mat() = rvb_t_mat * rvb_t1(0).mat();
    }

    hop12u = gate((s + 2) % Nsite, (s + 1) % Nsite);
    hop21d = gate((s + 1) % Nsite + Nsite, (s + 2) % Nsite + Nsite);

    if (!OPEN || s != Nsite - 2) {
      rvb_s_mat = hop12u - hop21d >> hubbard;
      rvb_t_mat = hop12u + hop21d >> hubbard;

      rvb_s2(0).mat() = rvb_s_mat * rvb_s2(0).mat();
      rvb_t2(0).mat() = rvb_t_mat * rvb_t2(0).mat();
    }
  }

  rvb_s1(0).mat() /= std::sqrt(Hermitian::norm(rvb_s1(0).mat()));
  rvb_s2(0).mat() /= std::sqrt(Hermitian::norm(rvb_s2(0).mat()));
  rvb_t1(0).mat() /= std::sqrt(Hermitian::norm(rvb_t1(0).mat()));
  rvb_t2(0).mat() /= std::sqrt(Hermitian::norm(rvb_t2(0).mat()));

  cout << abs((rvb_s1(0).mat().adjoint() * rvb_s2(0).mat())(0, 0)) << "\t"
       << abs((rvb_t1(0).mat().adjoint() * rvb_t2(0).mat())(0, 0)) << endl;
  // finished RVB

  SpMatrix phi_mat(ndim, ndim);
  SpMatrix eta_mat(ndim, ndim);

  phi_mat.setZero();
  eta_mat.setZero();

  vector<cplx> rvbs_correl(Nsite), rvbt_correl(Nsite);
  Matrix rvb_s = rvb_s1(0).mat() + rvb_s2(0).mat();
  Matrix rvb_t = rvb_t1(0).mat() + rvb_t2(0).mat();
  rvb_s /= sqrt(Hermitian::norm(rvb_s));
  rvb_t /= sqrt(Hermitian::norm(rvb_t));

  for (auto s : range(Nsite)) {
    for (auto sp : range(Nsite)) {
      if (s != sp) {
        int site1u = s;
        int site2u = sp;
        auto hop21u = gate(site2u, site1u);
        auto hop12u = gate(site1u, site2u);

        int site1d = s + Nsite;
        int site2d = sp + Nsite;
        auto hop21d = gate(site2d, site1d);
        auto hop12d = gate(site1d, site2d);

        auto pair_mat = hop21u * hop21d >> hubbard;
        if (sp == 0) {
          mc.updateMatrix("pair_correl", s, pair_mat);
          mc.updateMatrix("sxy_correl", s, -(hop21u * hop12d >> hubbard));

          rvbs_correl[s] = (rvb_s.adjoint() * pair_mat * rvb_s)(0, 0);
          rvbt_correl[s] = (rvb_t.adjoint() * pair_mat * rvb_t)(0, 0);
        }

        phi_mat += pair_mat;
        eta_mat += ((s - sp) % 2 == 0 ? 1 : -1) * pair_mat;

        if (sp == 0 && s == sp + dh_sep) {
          rvb_s_mat = hop12u - hop21d >> hubbard;
          rvb_t_mat = hop12u + hop21d >> hubbard;

          if (dh_quench == 1)
            quench_mat = rvb_s_mat / sqrt(2.);
          else if (dh_quench == 2)
            quench_mat = rvb_t_mat / sqrt(2.);
          else if (dh_quench == 3)
            quench_mat = .5 * (rvb_s_mat + rvb_t_mat);
        }

      } else {
        auto pair_mat = gate(s, s) * gate(s + Nsite, s + Nsite) >> hubbard;
        phi_mat += pair_mat;
        eta_mat += pair_mat;

        if (sp == 0) {
          mc.updateMatrix("pair_correl", s, pair_mat);
          mc.updateMatrix(
              "sxy_correl", s,
              gate(s, s) - gate(s, s) * gate(s + Nsite, s + Nsite) >> hubbard);
          rvbs_correl[s] = (rvb_s.adjoint() * pair_mat * rvb_s)(0, 0);
          rvbt_correl[s] = (rvb_t.adjoint() * pair_mat * rvb_t)(0, 0);
        }
      }
    }
    auto ns = (gate(s, s) + gate(s + Nsite, s + Nsite));
    auto n0 = (gate(0, 0) + gate(Nsite, Nsite));
    mc.updateMatrix("cdw_correl", s, n0 * ns - n0 - ns >> hubbard);

    auto szs = gate(s, s) - gate(s + Nsite, s + Nsite);
    auto sz0 = (gate(0, 0) - gate(Nsite, Nsite));
    mc.updateMatrix("sz_correl", s, (szs * sz0 >> hubbard) / 4.0);
  }
  mc.updateMatrix("eta_exp", 0, eta_mat);
  mc.updateMatrix("phi_exp", 0, phi_mat);

  std::cout << "initial energy:" << Etot << std::endl;
  double err;
  int i = 1;
  if (solver_type == "full") {
    Htot = Tp + Tm + U0 * docc;
    mc.updateMatrix("J_exp", 0, II * thop * (Tp - Tm));

    mc.updateMatrix("Etot", 0, Htot);

    Eigen::SelfAdjointEigenSolver<Matrix> eigensolver(Htot);
    cout << "Full solution: start ground state" << endl;
    psi(0).mat() = eigensolver.eigenvectors().col(0);
    ofstream ef("e.out");
    for (auto i : range(eigensolver.eigenvalues().rows())) {
      ef << eigensolver.eigenvalues()(i, 0) << endl;
    }
    ef.close();
  } else {
    Htot = Tp + Tm + U0 * docc;
    mc.updateMatrix("J_exp", 0, II * thop * (Tp - Tm));

    mc.updateMatrix("Etot", 0, Htot);
    auto solver = Hermitian::getSolver(hubbard, kry_dim);
    cout << "start to compute the ground state" << endl;
    solver.getGroundState(Htot, psi, 0);
  }
  cout << "ground state computed" << endl;
  mc.clearAll();
  mc.incrAllObs(0);

  rvbs1[0] = (psi(0).mat().adjoint() * rvb_s1(0).mat())(0, 0);
  rvbs2[0] = (psi(0).mat().adjoint() * rvb_s2(0).mat())(0, 0);
  rvbt1[0] = (psi(0).mat().adjoint() * rvb_t1(0).mat())(0, 0);
  rvbt2[0] = (psi(0).mat().adjoint() * rvb_t2(0).mat())(0, 0);

  pcmod[0] = sqrt(Hermitian::norm(Pc * psi(0).mat()));

  std::cout << "ground state energy: " << Etot << endl;
  std::cout << "(check hermiticity) ground state energy: "
            << psi.eval(0, Htot.adjoint()) << endl;

  t2 = chrono::high_resolution_clock::now();
  std::cout << "Time elapsed: "
            << chrono::duration_cast<chrono::milliseconds>(t2 - t1).count() *
                   1e-3
            << "s" << std::endl;

#ifdef TUNIF
  cout << "docc distribution: ";
  for (auto i : range(Nsite)) {
    int site1u = i;
    int site1d = i + Nsite;
    auto nup = gate(site1u, site1u);
    auto ndo = gate(site1d, site1d);

    cout << "\t" << psi.eval(0, (nup * ndo) >> hubbard);
  }
  cout << endl;

  cout << "charge distribution: ";
  for (auto i : range(Nsite)) {
    int site1u = i;
    int site1d = i + Nsite;
    auto nup = gate(site1u, site1u);
    auto ndo = gate(site1d, site1d);

    cout << "\t" << psi.eval(0, (nup + ndo) >> hubbard);
  }
  cout << endl;

#endif
  if (dh_quench > 0) {
    psi(0).mat() = quench_mat * psi(0).mat();
    psi(0).mat() /= sqrt(Hermitian::norm(psi(0).mat()));
  }

  //

  // real-time
  auto l_solver = Arnoldi::getSolver(hubbard, kry_dim_rt);
  for (int imc = 0; imc < Nmc; ++imc) {
    cout << "Proc " << mc.report() << ", sample " << imc << endl;
    for (int i = 1; i <= nt; ++i) {
      Htot = Tdr[i] * Tp + std::conj(Tdr[i]) * Tm + (ndr[i] + U0) * docc;

      mc.updateMatrix("J_exp", 0,
                      II * thop * (Tdr[i] * Tp - std::conj(Tdr[i]) * Tm));

      mc.updateMatrix("Etot", 0, Htot);

      auto weights = lindblad_weight(0, Nsite, psi, gate, hubbard);
      vector<double> wt(weights.size() + 1);
      double wtot = 0.;
      wt[0] = 0.;
      for (int j = 1; j < Nsite + 1; ++j) wt[j] = wt[j - 1] + weights[j - 1];

      wtot = wt[wt.size() - 1];
      for (int j = 0; j < Nsite + 1; ++j) wt[j] /= wt[Nsite];
      if (dice(mt) < gamma * wtot * h)  // jump
      {
        auto dice_res = dice(mt);
        auto idx = bsearch(wt, dice_res);
        auto jump_op = gate(idx, idx) >> hubbard;
        psi(1).mat() = jump_op * (jump_op)*psi(0).mat();
        psi(1).mat() /= sqrt(Hermitian::norm(psi(1).mat()));
      } else  // evolve
      {
        QJ_operator(Htot, qj_tot, Id, gamma * wtot);

        l_solver.evolve(Htot, psi, 0, h);
        psi(1).mat() /= std::sqrt(Hermitian::norm(psi(1).mat()));
      }

      mc.incrAllObs(i, 1);
      rvbs1[i] = (psi(1).mat().adjoint() * rvb_s1(0).mat())(0, 0);
      rvbs2[i] = (psi(1).mat().adjoint() * rvb_s2(0).mat())(0, 0);
      rvbt1[i] = (psi(1).mat().adjoint() * rvb_t1(0).mat())(0, 0);
      rvbt2[i] = (psi(1).mat().adjoint() * rvb_t2(0).mat())(0, 0);

      pcmod[i] = sqrt(Hermitian::norm(Pc * psi(1).mat()));

      psi(0).mat() = psi(1).mat();
    }
  }
  mc.normAll();
  std::cout << "Finished computation. Collecting data..." << std::endl;
  mc.collectAll();
  std::cout << "Finished collecting. Writing data to file..." << std::endl;
  mc.outFile(
      h, "obs.out",
      std::vector<std::string>{"Etot", "d_exp", "phi_exp", "eta_exp", "J_exp"});
  mc.outFile(h, "pair_correl.out",
             std::vector<std::string>{"pair_correl", "cdw_correl", "sxy_correl",
                                      "sz_correl"});
  if (mc.report() == 0) {
    ofstream f("rvb.out");
    for (int i = 0; i <= nt; ++i) {
      f << i * h << "\t" << rvbs1[i].real() << "\t" << rvbs1[i].imag() << "\t"
        << rvbs2[i].real() << "\t" << rvbs2[i].imag() << "\t" << rvbt1[i].real()
        << "\t" << rvbt1[i].imag() << "\t" << rvbt2[i].real() << "\t"
        << rvbt2[i].imag() << "\t" << pcmod[i] << endl;
    }
    f.close();
    f.open("rvb_correl.out");
    for (int s = 0; s < Nsite; ++s)
      f << s << "\t" << rvbs_correl[s].real() << "\t" << rvbs_correl[s].imag()
        << "\t" << rvbt_correl[s].real() << "\t" << rvbt_correl[s].imag()
        << endl;
    f.close();
  }

  return 0;
}

void QJ_operator(SpMatrix& H, const SpMatrix& qj_tot, const SpMatrix& Id,
                 double wt) {
  H -= 0.5 * II * (qj_tot - wt * Id);
}

template <typename T, typename Tg>
vector<double> lindblad_weight(int tstp, int Nsite, const State<T>& psi,
                               Tg& hopgate, T& hubbard) {
  vector<double> res(Nsite);

  for (int i = 0; i < Nsite; ++i) {
    auto ni = lindblad(i, Nsite, hopgate, hubbard);
    res[i] = psi.eval(tstp, SpMatrix(ni.adjoint()) * ni).real();
  }
  return res;
}

template <typename GT, typename PT, typename IT>
void incr_correl(int Nsite, int i, const PT& psi, IT& cavity, GT& gate,
                 Matrix& pair_correl, Matrix& sxy_correl, Matrix& cdw_correl,
                 Matrix& sz_correl) {
  for (auto s : range(Nsite)) {
    if (s != 0) {
      int site1u = s;
      int site2u = 0;
      auto hop21u = gate(site2u, site1u);
      auto hop12u = gate(site1u, site2u);

      int site1d = s + Nsite;
      int site2d = Nsite;
      auto hop21d = gate(site2d, site1d);
      auto hop12d = gate(site1d, site2d);

      pair_correl(i, s) += psi.eval(i, hop21u * hop21d >> cavity);
      sxy_correl(i, s) += -psi.eval(i, hop21u * hop12d >> cavity);
    } else {
      pair_correl(i, s) +=
          psi.eval(i, gate(s, s) * gate(s + Nsite, s + Nsite) >> cavity);
      sxy_correl(i, s) += psi.eval(
          i, gate(s, s) - gate(s, s) * gate(s + Nsite, s + Nsite) >> cavity);
    }
    auto ns = (gate(s, s) + gate(s + Nsite, s + Nsite));
    auto n0 = (gate(0, 0) + gate(Nsite, Nsite));
    cdw_correl(i, s) += psi.eval(i, n0 * ns - n0 - ns >> cavity);

    auto szs = gate(s, s) - gate(s + Nsite, s + Nsite);
    auto sz0 = (gate(0, 0) - gate(Nsite, Nsite));
    sz_correl(i, s) += psi.eval(i, szs * sz0 >> cavity) / 4.0;
  }
}

int bsearch(std::vector<double>& values, double val) {
  if (val > values[values.size() - 1] || val < values[0]) return -1;
  int is = 0, ie = values.size();
  while (ie - is > 1) {
    int imid = is + (ie - is) / 2;
    if (val < values[imid])
      ie = imid;
    else
      is = imid;
  }

  return is;
}
