#pragma once

#include <unordered_map>
#ifdef USE_MPI
	#include "mpi.h"
#endif
namespace qudrip
{
#define MC_REAL 0
#define MC_COMPLEX 1
template<typename TS>
class McManager
{
	//std::vector<Matrix> data;
	public:
		using functor_type = std::function<value_type(const value_type &)>;
#ifdef USE_MPI
		McManager(TS* p_psi, int* argc, char*** argv):
			nt_(p_psi -> nt()), p_psi_(p_psi),
			Nmc_(p_psi -> nt(), 0)
		{
			MPI_Init(argc, argv);
			MPI_Comm_size(MPI_COMM_WORLD, &Nproc_);
			MPI_Comm_rank(MPI_COMM_WORLD, &my_rank_);
		}

		McManager(int nt, TS* p_psi, int* argc, char*** argv):
			nt_(nt), p_psi_(p_psi),
			Nmc_(nt, 0)
		{
			MPI_Init(argc, argv);
			MPI_Comm_size(MPI_COMM_WORLD, &Nproc_);
			MPI_Comm_rank(MPI_COMM_WORLD, &my_rank_);
		}
		~McManager()
		{
			MPI_Finalize();
		}
#else 

		McManager(int nt, TS* p_psi):
			nt_(nt), p_psi_(p_psi), 
			my_rank_(0), Nproc_(1),
			Nmc_(nt, 0)
		{}
		McManager(TS* p_psi):
			nt_(p_psi -> nt()), p_psi_(p_psi), 
			my_rank_(0), Nproc_(1),
			Nmc_(p_psi -> nt(), 0)
		{}
#endif

		void addObs(const std::string & name, int dim, int type = MC_COMPLEX)
		{
			data_.insert(std::make_pair(name, Matrix(nt_, dim)));
			obs_.insert(std::make_pair(name, std::vector<SpMatrix>(dim)));
			type_.insert(std::make_pair(name, type));
			aux_.insert(std::make_pair(name, std::vector<std::string>{}));
		}

		void addAux(const std::string & name, const std::string & name2, const functor_type & func, int dim, int aux_type = MC_COMPLEX)
		{
			if(obs_.find(name2) == obs_.end())
			{
				throw std::string{"addAux: wrong dependent variable!"};
			}
			data_.insert(std::make_pair(name, Matrix(nt_, dim)));
			func_.insert(std::make_pair(name, func));
			type_.insert(std::make_pair(name, aux_type));
			//aux_.insert(std::make_pair(name2, name));
			aux_.find(name2) -> second.push_back(name);
		}

		template<typename MT>
		void updateMatrix(const std::string & name, int idx, MT && mat)
		{
			obs_.find(name) -> second[idx] = mat;
		}

		void incrAllObs(int tstp, int psi_tstp = -1)
		{
			if (psi_tstp == -1) psi_tstp = tstp;
			for(auto & obs : obs_)
			{
				//std::cout << obs.first << std::endl;
				auto & var_mat = data_.find(obs.first) -> second;
				const auto & dep_vars = aux_.find(obs.first) -> second;
				for(int d = 0; d < var_mat.cols(); ++d)
				{
					auto term = p_psi_ -> eval(psi_tstp, obs.second[d]);
					var_mat(tstp, d) += term;
					//aux vars
					for(const auto & n : dep_vars)
					{
						auto & func = func_.find(n) -> second;
						data_.find(n) -> second(tstp, d) += func(term);
					}
				}
			}
			//incrAllAux(tstp);
			++Nmc_[tstp];
		}

		void clearAll(int tstp)
		{
			for(auto & data : data_)
			{
				for(int d = 0; d < (data.second).cols(); ++d)
				{
					data.second(tstp, d) = 0.;
				}
			}
		}

		void clearAll()
		{
			for(auto & data : data_)
			{
				data.second.setZero();
			}
		}

		void normAll()
		{
			for(auto & data : data_)
			{
				for(int t = 0; t < nt_; ++t)
				{
					data.second.block(t, 0, 1, data.second.cols()) *= 1. / Nmc_[t];
				}
			}
		}

		void collectAll()
		{
#ifdef USE_MPI
			for(auto & data : data_)
			{
				int size = data.second.cols() * data.second.rows();
				if(my_rank_ == 0)
				{
					MPI_Reduce(MPI_IN_PLACE, (data.second).data(), size, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD);
					data.second *= 1. / Nproc_;
				}
				else
				{
					MPI_Reduce((data.second).data(), (data.second).data(), size, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD);
				}
			}

#endif 
		}

	void outFile(double h, const std::string & fname, const std::vector<std::string> & obs_names, int precision = 6)
	{
		if(my_rank_ != 0) return;
		std::ofstream f(fname);
		auto it_d = data_.find(obs_names[0]);
		if(it_d == data_.end()) throw std::string{"outFile: wrong name!!"};
		int dim = it_d -> second.cols();

		std::string label;
		if(dim > 1) label = "#time\tindex";
		else label = "#time";

		for(auto & n : obs_names)
		{
			auto it_t = type_.find(n);
			if(it_t == type_.end())
			{
				throw std::string{"outFile: wrong name!!"};
			}
			if(it_t -> second == MC_REAL) label += "\t" + n;
			else label += "\t" + n + "(real)" + "\t" + n + "(imag)";
			
			auto it_d = data_.find(n);
			if(it_d -> second.cols() != dim) 
			{
				std::cout << "mismatched dimensions" << std::endl;
				return;
			}
		}
		f << label << std::endl;
		for(int t = 0; t < nt_; ++t)
		{
			for(int idx = 0; idx < dim; ++idx)
			{
				if(dim > 1) label = std::to_string(t * h) + "\t" + std::to_string(idx);
				else label =  std::to_string(t * h);

				for(auto & n : obs_names)
				{
					auto it = data_.find(n);
					label += "\t" + std::to_string((it -> second)(t, idx).real());
					if(type_.find(n) -> second == MC_COMPLEX) label += "\t" + std::to_string((it -> second)(t, idx).imag());
				}
				f << label << std::endl;
			}
			if(dim > 1) f << std::endl;
		}
		f.close();
	}

	value_type obsValue(const std::string & name, int tstp, int idx)
	{
		return data_.find(name) -> second(idx);
	}

	int report()
	{
		return my_rank_;
	}

	private:
/*
	void incrAllAux(int tstp)
	{
		for(auto & func : func_)
		{
			auto it = data_.find(func.first);
			const auto & aux = data_.find(aux_.find(func.first) -> second) -> second;
			for(int d = 0; d < (it -> second).cols(); ++d)
			{
				(it -> second)(tstp, d) = func.second(aux(tstp, d));
			}
		}
	}
*/
	TS * p_psi_;
	int nt_, Nproc_;
	int my_rank_;
	std::unordered_map<std::string, Matrix> data_;
	std::unordered_map<std::string, std::vector<SpMatrix>> obs_;
	std::unordered_map<std::string, int> type_;
	std::unordered_map<std::string, functor_type> func_;
	std::unordered_map<std::string, std::vector<std::string>> aux_;
	std::vector<size_t> Nmc_;
};

#ifdef USE_MPI
template<typename TS>
McManager<TS> getMC(TS & psi, int * argc, char*** argv)
{
	return McManager<TS>(&psi, argc, argv);
}

template<typename TS>
McManager<TS> getMC(int nt, TS & psi, int * argc, char*** argv)
{
	return McManager<TS>(nt, &psi, argc, argv);
}
#else
template<typename TS>
McManager<TS> getMC(TS & psi)
{
	return McManager<TS>(&psi);
}
template<typename TS>
McManager<TS> getMC(int nt, TS & psi)
{
	return McManager<TS>(nt, &psi);
}
#endif
}
