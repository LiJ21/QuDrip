#include<fstream>
#include<string>
#include<iostream>
#include<unordered_map>
#include <algorithm>

#ifndef FIND_PARAM
#define FIND_PARAM
std::string delete_space(std::string s)
{
	s.erase(std::remove_if(s.begin(), s.end(),
				[](char &c) {
				return std::isspace<char>(c, std::locale::classic());
				}),
			s.end());

	return s;
}

template<class T>
auto convert(const std::string & str)
{
	T res;
	try
	{
		res = static_cast<T>(std::stod(str));
	}
	catch(std::out_of_range & err)
	{
		res = static_cast<T>(std::stold(str));
	}

	return res;
		
}

class param_finder
{
	std::ifstream f;
	std::unordered_map<std::string, std::string> params;
	public:
	param_finder(std::string fname) : f(fname)
	{
		std::string oline;
		while(std::getline(f, oline))
		{
			auto line = delete_space(oline);
			auto peq = line.find("=");
			if(line[0] == '#' || peq == std::string::npos) 
			{
				continue;
			}
			else
			{
				params.insert({line.substr(0, peq), line.substr(peq + 1)});
			}
		}
	}
	template<class T>
	void find_param(std::string param, T & value, int warning = 0);

	template<class T>
	void find_param_tvec(std::string param, std::vector<T> & value, int warning = 0);

};

template<class T>
void param_finder::find_param(std::string param, T & value, int warning) 
{
	auto pval = params.find(param);
	if(pval == params.end())
	{
		if(warning) throw param;
	}
	else
	{
		value = convert<T>(pval -> second);
	}
}

template<>
auto convert<int>(const std::string & str)
{
	return stoi(str);
}

template<>
auto convert<std::string>(const std::string & str)
{
	return str;
}

template<>
auto convert<long>(const std::string & str)
{
	return stoi(str);
}

template<typename T>
void param_finder::find_param_tvec(std::string param, std::vector<T> & value, int warning) 
{
	auto pval = params.find(param);

	if(pval == params.end())
	{
		if(warning) throw param;
	}
	else
	{
		if((pval -> second).substr(0, 2) == "--")
		{
			std::ifstream f((pval -> second).substr(2));
			std::string line;
			for(auto t = 0; t < value.size(); ++t)
			{
				while(std::getline(f, line) && line[0] == '#');
				value[t] = convert<T>(line);	
			}

		}
		else 
		{
			auto v = convert<T>(pval -> second);
			std::fill(value.begin(), value.end(), v);
		}
	}
}

template<>
void param_finder::find_param_tvec<std::complex<double>>(std::string param, std::vector<std::complex<double>> & value, int warning) 
{
	auto pval = params.find(param);
	auto II = std::complex<double>(0.0, 1.0);

	if(pval == params.end())
	{
		if(warning) throw param;
	}
	else
	{
		if((pval -> second).substr(0, 2) == "--")
		{
			std::ifstream f((pval -> second).substr(2));
			std::string line;
			for(auto t = 0; t < value.size(); ++t)
			{
				while(std::getline(f, line) && line[0] == '#');
				double re, im;
				std::stringstream iss(line);
				iss >> re >> im;
				value[t] = re + II * im;	
			}

		}
		else 
		{
			std::stringstream iss(pval -> second);
			auto v = std::complex<double>(0.0, 0.0);
			iss >> v;
			std::fill(value.begin(), value.end(), v);
		}
	}
}

#endif
