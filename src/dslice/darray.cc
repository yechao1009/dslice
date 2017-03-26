#include <fstream>
#include <sstream>
#include "../../include/dslice/darray.h"

namespace dslice {

bool operator==(const IntDarray& lhs, const IntDarray& rhs) {
	if (lhs.Nrow() != rhs.Nrow() || lhs.Ncol() != rhs.Ncol()) {
		return false;
	}
	ds_uint nrow = lhs.Nrow();
	ds_uint ncol = lhs.Ncol();
	for (int i = 0; i != nrow; ++i) {
		for (int j = 0; j != ncol; ++j) {
			if (lhs[i][j] != rhs[i][j]) {
				return false;
			}
		}
	}
	return true;
}

bool operator!=(const IntDarray& lhs, const IntDarray& rhs) {
	return !(lhs == rhs);
}

const IntDarray operator+(const IntDarray& lhs, const IntDarray& rhs) {
	IntDarray int_darray;
	if (lhs.Nrow() != rhs.Nrow() || lhs.Ncol() != rhs.Ncol()) {
		return int_darray;
	}
	int_darray = lhs;
	int_darray += rhs;
	return int_darray;
}

const IntDarray operator-(const IntDarray& lhs, const IntDarray& rhs) {
	IntDarray int_darray;
	if (lhs.Nrow() != rhs.Nrow() || lhs.Ncol() != rhs.Ncol()) {
		return int_darray;
	}
	int_darray = lhs;
	int_darray -= rhs;
	return int_darray;
}

const IntDarray& IntDarray::operator+=(const IntDarray& int_darray) {
	if (Nrow() != int_darray.Nrow() || Ncol() != int_darray.Ncol()) {
		return *this;
	}
	ds_uint nrow = Nrow();
	ds_uint ncol = Ncol();
	for (int i = 0; i != nrow; ++i) {
		for (int j = 0; j != ncol; ++j) {
			darray_[i][j] += int_darray[i][j];
		}
	}
	return *this;
}

const IntDarray& IntDarray::operator-=(const IntDarray& int_darray) {
	if (Nrow() != int_darray.Nrow() || Ncol() != int_darray.Ncol()) {
		return *this;
	}
	ds_uint nrow = Nrow();
	ds_uint ncol = Ncol();
	for (int i = 0; i != nrow; ++i) {
		for (int j = 0; j != ncol; ++j) {
			darray_[i][j] -= int_darray[i][j];
		}
	}
	return *this;
}

char StringToChar(const std::string& in_str) {
	const char* cstr = in_str.c_str();
	return cstr[0];
}

int StringToInt(const std::string& in_str) {
	int out_value = 0;
	if (~in_str.empty()) {
		std::stringstream converter(in_str);
		converter >> out_value;
	}
	return out_value;
}

const IntDarray ReadDarray(std::string file, std::string delim, bool header, bool with_row_name) {
	std::ifstream file_in(file.c_str());
	IntDarray int_darray;
	if (!file_in) {
		std::cerr << "Error: fail to open file " << file << std::endl;
		return int_darray;
	}
	char sep = StringToChar(delim);
	int_darray = ReadDarray(file_in, sep, header, with_row_name);
	file_in.close();
	return int_darray;
}

const IntDarray ReadDarray(std::istream& in_stream, char delim, bool header, bool with_row_name) {
	IntDarray int_darray;
	std::string buffer, str;
	int value;
	std::vector<int> new_row;
	std::vector<std::string> col_names;
	std::vector<std::string> row_names;

	if (header) {
		std::getline(in_stream, buffer);
		std::istringstream istr_stream(buffer);
		if (with_row_name) {
			std::getline(istr_stream, str, delim);
		}
		while (std::getline(istr_stream, str, delim)) {
			col_names.push_back(str);
		}
	}
	while (std::getline(in_stream, buffer)) {
		std::istringstream istr_stream(buffer);
		new_row.resize(0);
		if (with_row_name) {
			std::getline(istr_stream, str, delim);
			row_names.push_back(str);
		}
		while (std::getline(istr_stream, str, delim)) {
			value = StringToInt(str);
			new_row.push_back(value);
		}
		if (new_row.size() == 0) {
			continue;
		}
		if (!int_darray.PushBack(new_row)) {
			int_darray.Resize(0, 0);
			return int_darray;
		}
	}
	if (header) {
		int_darray.SetDimName(col_names, 2);
	}
	if (with_row_name) {
		int_darray.SetDimName(row_names, 1);
	}
	return int_darray;
}

void WriteDarray(const IntDarray& int_darray, std::string file, std::string delim, bool col_name, bool row_name) {
	std::ofstream file_out(file.c_str());
	if (!file_out) {
		return;
	}
	WriteDarray(int_darray, file_out, delim, col_name, row_name);
	file_out.close();
}

void WriteDarray(const IntDarray& int_darray, std::ostream& out_stream, std::string delim, bool col_name, bool row_name) {
	if (int_darray.Empty()) {
		return;
	}
	ds_uint nrow = int_darray.Nrow();
	ds_uint ncol = int_darray.Ncol();
	std::vector<std::string> rownames = int_darray.RowName();
	std::vector<std::string> colnames = int_darray.ColName();
	row_name = row_name && (!rownames.empty());
	col_name = col_name && (!colnames.empty());

	if (row_name && col_name) {
		out_stream << "Darray" << delim;
	}
	if (col_name) {
		out_stream << colnames[0];
		for (ds_uint j = 1; j != ncol; ++j) {
			out_stream << delim << colnames[j];
		}
		out_stream << std::endl;
	}
	if (row_name) {
		for (ds_uint i = 0; i != nrow; ++i) {
			out_stream << rownames[i];
			for (ds_uint j = 0; j != ncol; ++j) {
				out_stream << delim << int_darray(i, j);
			}
			out_stream << std::endl;
		}
	} else {
		for (ds_uint i = 0; i != nrow; ++i) {
			out_stream << int_darray(i, 0);
			for (ds_uint j = 1; j != ncol; ++j) {
				out_stream << delim << int_darray(i, j);
			}
			out_stream << std::endl;
		}
	}
}

} // namespace

