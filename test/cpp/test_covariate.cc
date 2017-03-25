#include <iostream>
#include "../../include/dslice/covariate.h"

int main(int argc, char** argv)
{

//	test integer type covariate
	std::vector<int> int_vec(4, -1);
	int_vec[0] = 0;
	int_vec[1] = 2;
	int_vec[2] = 1;
	int_vec[3] = 2;
	dslice::Covariate<int> int_cov(int_vec, false);

	std::vector<int> buf_s_1 = int_cov.GetSample();
	std::vector<dslice::ds_uint> buf_f_1 = int_cov.GetFactor();
	std::map<int, dslice::ds_uint> buf_d_1 = int_cov.GetDict();

	std::cout << "integer type covariate\n";
	std::cout << "sample: ";
	for (std::vector<int>::iterator it = buf_s_1.begin(); it != buf_s_1.end(); ++it) {
		std::cout << *it << " ";
	}
	std::cout << std::endl;
	std::cout << "factor: ";
	for (std::vector<dslice::ds_uint>::iterator it = buf_f_1.begin(); it != buf_f_1.end(); ++it) {
		std::cout << *it << " ";
	}
	std::cout << std::endl;
	std::cout << "size: " << int_cov.GetSize() << std::endl;
	std::cout << "level: " << int_cov.GetLevel() << std::endl;
	std::cout << "map:\n";
	for (std::map<int, dslice::ds_uint>::iterator it = buf_d_1.begin(); it != buf_d_1.end(); ++it) {
		std::cout << it->first << "=" << it->second << std::endl;
	}
	std::cout << std::endl;

//	test string type covariate
	std::vector<std::string> str_vec(4, "");
	str_vec[0] = "cat";
	str_vec[1] = "bird";
	str_vec[2] = "bird";
	str_vec[3] = "dog";
	dslice::Covariate<std::string> str_cov(str_vec);

	std::vector<std::string> buf_s_2 = str_cov.GetSample();
	std::vector<dslice::ds_uint> buf_f_2 = str_cov.GetFactor();
	std::map<std::string, dslice::ds_uint> buf_d_2 = str_cov.GetDict();

	std::cout << "string type covariate\n";
	std::cout << "sample: ";
	for (std::vector<std::string>::iterator it = buf_s_2.begin(); it != buf_s_2.end(); ++it) {
		std::cout << *it << " ";
	}
	std::cout << std::endl;
	std::cout << "factor: ";
	for (std::vector<dslice::ds_uint>::iterator it = buf_f_2.begin(); it != buf_f_2.end(); ++it) {
		std::cout << *it << " ";
	}
	std::cout << std::endl;
	std::cout << "size: " << str_cov.GetSize() << std::endl;
	std::cout << "level: " << str_cov.GetLevel() << std::endl;
	std::cout << "map:\n";
	for (std::map<std::string, dslice::ds_uint>::iterator it = buf_d_2.begin(); it != buf_d_2.end(); ++it) {
		std::cout << it->first << "=" << it->second << std::endl;
	}
	std::cout << std::endl;

	return 0;
}
