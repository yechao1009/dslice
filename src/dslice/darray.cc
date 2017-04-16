#include "../../include/dslice/darray.h"

namespace dslice {

const Darray<ds_int> StringDarray::ToInt() {
	Darray<ds_int> int_darray;
	ds_uint nrow = Nrow();
	ds_uint ncol = Ncol();
	int_darray.Resize(nrow, ncol);
	std::string buffer;
	for (int i = 0; i != nrow; ++i) {
		for (int j = 0; j != ncol; ++j) {
			buffer = darray_[i][j];
			if (buffer.empty()) {
				int_darray.SetValue(i, j, 0);
			} else {
				int_darray.SetValue(i, j, stoi(buffer));
			}
		}
	}
	if (row_name_.size() != nrow) {
		row_name_.resize(nrow);
	}
	int_darray.SetDimName(row_name_, 1);
	if (col_name_.size() != ncol) {
		col_name_.resize(ncol);
	}
	int_darray.SetDimName(col_name_, 2);
	return int_darray;
}

const Darray<ds_uint> StringDarray::ToUInt() {
	Darray<ds_uint> int_darray;
	ds_uint nrow = Nrow();
	ds_uint ncol = Ncol();
	int_darray.Resize(nrow, ncol);
	std::string buffer;
	for (int i = 0; i != nrow; ++i) {
		for (int j = 0; j != ncol; ++j) {
			buffer = darray_[i][j];
			if (buffer.empty()) {
				int_darray.SetValue(i, j, 0);
			} else {
				int_darray.SetValue(i, j, stoul(buffer));
			}
		}
	}
	if (row_name_.size() != nrow) {
		row_name_.resize(nrow);
	}
	int_darray.SetDimName(row_name_, 1);
	if (col_name_.size() != ncol) {
		col_name_.resize(ncol);
	}
	int_darray.SetDimName(col_name_, 2);
	return int_darray;
}

const Darray<ds_double> StringDarray::ToDouble() {
	Darray<ds_double> double_darray;
	ds_uint nrow = Nrow();
	ds_uint ncol = Ncol();
	double_darray.Resize(nrow, ncol);
	std::string buffer;
	for (int i = 0; i != nrow; ++i) {
		for (int j = 0; j != ncol; ++j) {
			buffer = darray_[i][j];
			if (buffer.empty()) {
				double_darray.SetValue(i, j, 0);
			} else {
				double_darray.SetValue(i, j, stod(buffer));
			}
		}
	}
	if (row_name_.size() != nrow) {
		row_name_.resize(nrow);
	}
	double_darray.SetDimName(row_name_, 1);
	if (col_name_.size() != ncol) {
		col_name_.resize(ncol);
	}
	double_darray.SetDimName(col_name_, 2);
	return double_darray;
}

char StringToChar(const std::string& in_str) {
	const char* cstr = in_str.c_str();
	return cstr[0];
}

StringDarray ReadDarray(std::string file, std::string delim, bool header, bool with_row_name) {
	const char* in_file = file.c_str();
	char sep = StringToChar(delim);
	dslice::StringDarray darray = ReadDarray(in_file, sep, header, with_row_name);
	return darray;
}

StringDarray ReadDarray(const char* file, char delim, bool header, bool with_row_name) {
	StringDarray str_darray;
	std::ifstream file_in(file);
	if (!file_in) {
		std::cerr << "Error: fail to open file " << file << std::endl;
		return str_darray;
	}
	str_darray = ReadDarray(file_in, delim, header, with_row_name);
	file_in.close();
	return str_darray;
}

StringDarray ReadDarray(std::istream& in_stream, char delim, bool header, bool with_row_name) {
	StringDarray str_darray;
	std::string buffer, str;
	std::vector<std::string> new_row;
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
			new_row.push_back(str);
		}
		if (new_row.size() == 0) {
			continue;
		}
		if (!str_darray.PushBack(new_row)) {
			str_darray.Resize(0, 0);
			return str_darray;
		}
	}
	if (header) {
		str_darray.SetDimName(col_names, 2);
	}
	if (with_row_name) {
		str_darray.SetDimName(row_names, 1);
	}
	return str_darray;
}

Darray<ds_int> ReadIntDarray(std::string file, std::string delim, bool header, bool with_row_name) {
	const char* in_file = file.c_str();
	char sep = StringToChar(delim);
	dslice::Darray<ds_int> darray = ReadIntDarray(in_file, sep, header, with_row_name);
	return darray;
}

Darray<ds_int> ReadIntDarray(const char* file, char delim, bool header, bool with_row_name) {
	dslice::StringDarray buffer = ReadDarray(file, delim, header, with_row_name);
	dslice::Darray<ds_int> darray = buffer.ToInt();
	return darray;
}

Darray<ds_uint> ReadUIntDarray(std::string file, std::string delim, bool header, bool with_row_name) {
	const char* in_file = file.c_str();
	char sep = StringToChar(delim);
	dslice::Darray<ds_uint> darray = ReadUIntDarray(in_file, sep, header, with_row_name);
	return darray;
}

Darray<ds_uint> ReadUIntDarray(const char* file, char delim, bool header, bool with_row_name) {
	dslice::StringDarray buffer = ReadDarray(file, delim, header, with_row_name);
	dslice::Darray<ds_uint> darray = buffer.ToUInt();
	return darray;
}

Darray<ds_double> ReadDoubleDarray(std::string file, std::string delim, bool header, bool with_row_name) {
	const char* in_file = file.c_str();
	char sep = StringToChar(delim);
	dslice::Darray<ds_double> darray = ReadDoubleDarray(in_file, sep, header, with_row_name);
	return darray;
}

Darray<ds_double> ReadDoubleDarray(const char* file, char delim, bool header, bool with_row_name) {
	dslice::StringDarray buffer = ReadDarray(file, delim, header, with_row_name);
	dslice::Darray<ds_double> darray = buffer.ToDouble();
	return darray;
}

} // namespace
