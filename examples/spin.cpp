#include "qudrip.hpp"
//#include "operator_abel.hpp"
#include<iostream> 
#include<fstream>
#include<chrono>
#include<array>
#include<random>

using namespace qudrip;
using namespace std;

int main(int argc, char** argv)
{
	auto Nsite = std::stoi(argv[1]);
	auto prechain = getQbits(Nsite);
	auto Nset = getNset(prechain);
	auto gc = getBoseGate(prechain);
	auto ga = getBoseGate(prechain);

	std::cout << "Constraining..." << std::endl;
	auto chain = Constrain(prechain, Nset = Nsite / 2);
	auto ndim = chain.range();
	std::cout << "dimension of Hil space = " << ndim << endl;

	Matrix Uc = 0.5 * (pauli_x + II * pauli_y);
	Matrix Ua = 0.5 * (pauli_x - II * pauli_y);
	gc << Uc;
	ga << Ua;

	cout << "Start to construct" << endl;

	auto c1 = get_clock();
	c1.note();
	auto ham = getEmptyOperator();
	for(int i = 0; i < Nsite - 1; ++i)
	{
		ham = ham + gc(i + 1) * ga(i) + gc(i) * ga(i + 1);
	}

	auto ham_mat = ham >> chain;
	cout << "Time elapsed: " << c1.release() << std::endl;

	//solving the ground state
	int krylov_dim = 100;
	auto psi = getState(chain, 1);
	shuffle_state(psi(0));
	auto solver = Arnoldi::getLanczosSolver(chain, krylov_dim);

	solver.getGroundState(ham_mat, psi(0), 50);

	cout << "E = " << psi.eval(ham_mat) << endl;
	return 0;
}
