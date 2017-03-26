#include "../../include/dslice/darray.h"

int main(int argc, char** argv)
{

	int d_row = 4;
	int d_col = 2;
	std::cout << "test constructor with args: " << std::endl;
	dslice::IntDarray int_darray(d_row, d_col);
	for (int i = 0; i != d_row; ++i) {
		for (int j = 0; j != d_col; ++j) {
			int_darray.SetValue(i, j, i * 10 + j);
		}
	}
	for (int i = 0; i != d_row; ++i) {
		for (int j = 0; j != d_col; ++j) {
			std::cout << int_darray(i, j) << ' ';
		}
		std::cout << std::endl;
	}
	std::cout << "test constructor:" << std::endl;
	dslice::IntDarray test_array = int_darray;
	for (int i = 0; i != d_row; ++i) {
		for (int j = 0; j != d_col; ++j) {
			std::cout << test_array(i, j) << ' ';
		}
		std::cout << std::endl;
	}
	std::cout << "test sum:" << std::endl;
	std::vector<int> sum_vector = int_darray.Sum(2);
	for (int i = 0; i != sum_vector.size(); ++i) {
		std::cout << sum_vector[i] << ' ';
	}
	std::cout << std::endl;

	dslice::IntDarray test_array2;

	test_array2 = int_darray;
	std::cout << "test =" << std::endl;
	for (int i = 0; i != test_array2.Nrow(); ++i) {
		for (int j = 0; j != test_array2.Ncol(); ++j) {
			std::cout << test_array2(i, j) << ' ';
		}
		std::cout << std::endl;
	}

	std::cout << "test ==:" << std::endl;
	bool kk = int_darray == test_array;
	std::cout << kk << std::endl;

	std::string in_file = "test.txt";
	std::string delim = "\t";
	test_array = dslice::ReadDarray(in_file, delim, true, true);
	std::vector<std::string> rownames = test_array.RowName();
	std::vector<std::string> colnames = test_array.ColName();

	std::cout << "darray column name: " << std::endl;
	for (int j = 0; j != colnames.size(); ++j) {
		std::cout << colnames[j] << ' ';
	}
	std::cout << std::endl;

	std::cout << "darray row name: " << std::endl;
	for (int i = 0; i != rownames.size(); ++i) {
		std::cout << rownames[i] << ' ';
	}
	std::cout << std::endl;

	std::cout << "darray: " << std::endl;
	for (int i = 0; i != test_array.Nrow(); ++i) {
		for (int j = 0; j != test_array.Ncol(); ++j) {
			std::cout << test_array(i, j) << ' ';
		}
		std::cout << std::endl;
	}

	std::string out_file = "test.csv";
	delim = ",";
	dslice::WriteDarray(test_array, out_file, delim, true, true);

	return 0;
}
