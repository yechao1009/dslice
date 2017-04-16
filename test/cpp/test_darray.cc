#include "../../include/dslice/darray.h"

int main(int argc, char** argv)
{
	std::string delim = "\t";
	dslice::Darray<double> test_array1 = dslice::ReadDoubleDarray("test.txt", "\t", true, true);
	dslice::Darray<double> test_array2 = test_array1;
	dslice::Darray<double> test_array;
	test_array = test_array1 + test_array2;
	test_array1.Print();

	std::string k;
	std::vector<double> testv = test_array1(2, k);
	std::cout << "darray row 3 plus 4.9: " << std::endl;
	for (int j = 0; j != test_array1.Ncol(); ++j) {
		testv[j] += 4.9;
		std::cout << testv[j] << ' ';
	}
	std::cout << std::endl;

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
	std::string out_file = "test.csv";
	delim = ",";
	dslice::WriteDarray(test_array1, out_file, delim, true, true);
	return 0;
}
