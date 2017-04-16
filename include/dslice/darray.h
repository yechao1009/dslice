#ifndef DSLICE_DARRAY_H_
#define DSLICE_DARRAY_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include "base.h"

namespace dslice
{

// two dimensional data array
template <typename Object>
class Darray {

public:
	explicit Darray() : darray_(0) {}

	Darray(ds_uint rows, ds_uint cols) : darray_(rows) {
		for (ds_uint i = 0; i != rows; ++i) {
			darray_[i].resize(cols);
		}
		row_name_.resize(rows);
		col_name_.resize(cols);
	}
	Darray(const Darray<Object>& darray) { *this = darray; }

	ds_uint Nrow() const { return darray_.size(); }
	ds_uint Ncol() const { return Nrow() ? (darray_[0].size()) : 0; }
	bool Empty() const { return Nrow() == 0; }

	const std::vector< std::vector<Object> > GetValue() const { return darray_; };
	void SetValue(ds_uint row_id, ds_uint col_id, Object value);
	void SetValue(const std::vector< std::vector<Object> >& value);
	bool SetDimName(const std::vector<std::string>& dim_name, int type);

	std::vector<std::string> RowName() const { return row_name_; };
	std::vector<std::string> ColName() const { return col_name_; };
	std::vector<Object> Sum(int direction) const;
	std::vector< std::vector<Object> > CumSum(int direction) const;
	void Resize(ds_uint rows, ds_uint cols);
	bool PushBack(const std::vector<Object>& vkt);
	void Print();

	const std::vector<Object>& operator[](ds_uint row) const { return darray_[row]; };
	std::vector<Object>& operator[](ds_uint row) { return darray_[row]; };
	const Object& operator()(ds_uint row, ds_uint col) const { return darray_[row][col]; };
	Object& operator()(ds_uint row, ds_uint col) { return darray_[row][col]; };
	std::vector<Object> operator()(ds_uint row, std::string index);
	std::vector<Object> operator()(std::string index, ds_uint col);
	const Darray<Object>& operator+=(const Darray<Object>& darray);
	const Darray<Object>& operator-=(const Darray<Object>& darray);

protected:
	std::vector< std::vector<Object> > darray_;
	std::vector<std::string> row_name_;
	std::vector<std::string> col_name_;
};

template <typename Object>
void Darray<Object>::SetValue(ds_uint row_id, ds_uint col_id, Object value) {
	if (row_id < darray_.size() && col_id < darray_[0].size()) {
		darray_[row_id][col_id] = value;
	}
}

template <typename Object>
void Darray<Object>::SetValue(const std::vector< std::vector<Object> >& value) {
	ds_uint nrow = value.size();
	ds_uint ncol = value[0].size();
	Resize(nrow, ncol);
	for (size_t i = 0; i != nrow; ++i) {
		for (size_t j = 0; j != ncol; ++j) {
			darray_[i][j] = value[i][j];
		}
	}
}

template <typename Object>
bool Darray<Object>::SetDimName(const std::vector<std::string>& dim_name, int type) {
	if (type == 2 && Ncol() == dim_name.size()) {
		col_name_.assign(dim_name.begin(), dim_name.end());
	} else if (type == 1 && Nrow() == dim_name.size()) {
		row_name_.assign(dim_name.begin(), dim_name.end());
	} else {
		std::cout << "Input vector length does not unmatch darray." << std::endl;
		return false;
	}
	return true;
}

template <typename Object>
void Darray<Object>::Resize(ds_uint rows, ds_uint cols) {
	ds_uint nrow = this->Nrow();
	ds_uint ncol = this->Ncol();
	if (nrow == rows && ncol == cols) {
		return;
	} else if (nrow == rows && ncol != cols) {
		for (ds_uint i = 0; i != rows; ++i) {
			darray_[i].resize(cols);
		}
	} else if (nrow != rows && ncol == cols) {
		darray_.resize(rows);
		for (ds_uint i = nrow; i < rows; ++i) {
			darray_[i].resize(cols);
		}		
	} else {
		darray_.resize(rows);
		for (ds_uint i = 0; i != rows; ++i) {
			darray_[i].resize(cols);
		}
	}
}

template <typename Object>
bool Darray<Object>::PushBack(const std::vector<Object>& vkt) {
	if (Nrow() == 0 || Ncol() == (ds_uint)vkt.size()) {
		darray_.push_back(vkt);
	} else {
		return false;
	}
	return true;
}

template <typename Object>
std::vector<Object> Darray<Object>::Sum(int direction) const {
	ds_uint nrow = this->Nrow();
	ds_uint ncol = this->Ncol();
	std::vector<Object> margin_sum;
	if (direction == 1) {
		std::vector<Object> sum_by_row(nrow, 0);
		for (size_t i = 0; i != nrow; ++i) {
			for (size_t j = 0; j != ncol; ++j) {
				sum_by_row[i] += darray_[i][j];
			}
		}
		margin_sum.assign(sum_by_row.begin(), sum_by_row.end());
	} else if (direction == 2) {
		std::vector<Object> sum_by_col(ncol, 0);
		for (size_t j = 0; j != ncol; ++j) {
			for (size_t i = 0; i != nrow; ++i) {
				sum_by_col[j] += darray_[i][j];
			}
		}
		margin_sum.assign(sum_by_col.begin(), sum_by_col.end());
	} else {
		std::cout << "Invalid type of sum" << std::endl;
	}
	return margin_sum;
}

template <typename Object>
std::vector< std::vector<Object> > Darray<Object>::CumSum(int direction) const {
	ds_uint nrow = this->Nrow();
	ds_uint ncol = this->Ncol();
	std::vector< std::vector<Object> > cum_sum = darray_;
	if (direction == 2) {
		for (ds_uint i = 1; i != nrow; ++i) {
			for (ds_uint j = 0; j != ncol; ++j) {
				cum_sum[i][j] += cum_sum[i-1][j];
			}
		}
	} else if (direction == 1) {
		for (ds_uint j = 1; j != ncol; ++j) {
			for (ds_uint i = 0; i != nrow; ++i) {
				cum_sum[i][j] += cum_sum[i][j-1];
			}
		}
	} else {
		std::cout << "Invalid type of sum" << std::endl;
	}
	return cum_sum;
}

template <typename Object>
void Darray<Object>::Print() {
	std::cout << "Darray";
	for (ds_uint j = 0; j < Ncol(); ++j) {
		std::cout << '\t' << col_name_[j];
	}
	std::cout << std::endl;
	for (ds_uint i = 0; i != Nrow(); ++i) {
		std::cout << row_name_[i];
		for (ds_uint j = 0; j < Ncol(); ++j) {
			std::cout << '\t' << darray_[i][j];
		}
		std::cout << std::endl;
	}
}

template <typename Object>
void WriteDarray(const Darray<Object>& darray, std::ostream& out_stream, std::string delim, bool col_name, bool row_name) {
	if (darray.Empty()) {
		return;
	}
	ds_uint nrow = darray.Nrow();
	ds_uint ncol = darray.Ncol();
	std::vector<std::string> rownames = darray.RowName();
	std::vector<std::string> colnames = darray.ColName();
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
				out_stream << delim << darray(i, j);
			}
			out_stream << std::endl;
		}
	} else {
		for (ds_uint i = 0; i != nrow; ++i) {
			out_stream << darray(i, 0);
			for (ds_uint j = 1; j != ncol; ++j) {
				out_stream << delim << darray(i, j);
			}
			out_stream << std::endl;
		}
	}
}

template <typename Object>
void WriteDarray(const Darray<Object>& darray, std::string file, std::string delim, bool col_name, bool row_name) {
	std::ofstream file_out(file.c_str());
	if (!file_out) {
		return;
	}
	WriteDarray(darray, file_out, delim, col_name, row_name);
	file_out.close();
}


template <typename Object>
std::vector<Object> Darray<Object>::operator()(ds_uint row, std::string index) {
	std::vector<Object> row_vector;
	if (index.empty()) {
		return darray_[row];
	} else {
		return row_vector;
	}
}

template <typename Object>
std::vector<Object> Darray<Object>::operator()(std::string index, ds_uint col) {
	std::vector<Object> column_vector;
	if (index.empty()) {
		column_vector.resize(Nrow());
		for (size_t i = 0; i != Nrow(); ++i) {
			column_vector[i] = darray_[i][col];
		}
	}
	return column_vector;
}

template <typename Object>
bool operator==(const Darray<Object>& lhs, const Darray<Object>& rhs) {
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

template <typename Object>
bool operator!=(const Darray<Object>& lhs, const Darray<Object>& rhs) {
	return !(lhs == rhs);
}

template <typename Object>
const Darray<Object>& Darray<Object>::operator+=(const Darray<Object>& darray) {
	if (Empty()) {
		*this = darray;
		return *this;
	}
	ds_uint nrow = Nrow();
	ds_uint ncol = Ncol();
	if (nrow != darray.Nrow() || ncol != darray.Ncol()) {
		return *this;
	}
	for (int i = 0; i != nrow; ++i) {
		for (int j = 0; j != ncol; ++j) {
			darray_[i][j] += darray[i][j];
		}
	}
	return *this;
}

template <typename Object>
const Darray<Object>& Darray<Object>::operator-=(const Darray<Object>& darray) {
	if (Empty()) {
		*this = darray;
		return *this;
	}
	ds_uint nrow = Nrow();
	ds_uint ncol = Ncol();
	if (nrow != darray.Nrow() || ncol != darray.Ncol()) {
		return *this;
	}
	for (int i = 0; i != nrow; ++i) {
		for (int j = 0; j != ncol; ++j) {
			darray_[i][j] -= darray[i][j];
		}
	}
	return *this;
}

template <typename Object>
const Darray<Object> operator+(const Darray<Object>& lhs, const Darray<Object>& rhs) {
	Darray<Object> darray;
	if (lhs.Nrow() != rhs.Nrow() || lhs.Ncol() != rhs.Ncol()) {
		return darray;
	}
	darray = lhs;
	darray += rhs;
	return darray;
}

template <typename Object>
const Darray<Object> operator-(const Darray<Object>& lhs, const Darray<Object>& rhs) {
	Darray<Object> darray;
	if (lhs.Nrow() != rhs.Nrow() || lhs.Ncol() != rhs.Ncol()) {
		return darray;
	}
	darray = lhs;
	darray -= rhs;
	return darray;
}

class StringDarray : public Darray<std::string> {

public:
	StringDarray() : Darray<std::string>() {};
	StringDarray(ds_uint rows, ds_uint cols) : Darray<std::string>(rows, cols) {};
	StringDarray(const StringDarray& str_darray) { *this = str_darray; }

	const Darray<ds_int> ToInt();
	const Darray<ds_uint> ToUInt();
	const Darray<ds_double> ToDouble();
};

char StringToChar(const std::string& in_str);
StringDarray ReadDarray(std::istream& in_stream, char delim, bool header = true, bool with_row_name = false);
StringDarray ReadDarray(std::string file, std::string delim, bool header = true, bool with_row_name = false);
StringDarray ReadDarray(const char* file, char delim, bool header = true, bool with_row_name = false);
Darray<ds_int> ReadIntDarray(std::string file, std::string delim, bool header = true, bool with_row_name = false);
Darray<ds_int> ReadIntDarray(const char* file, char delim, bool header = true, bool with_row_name = false);
Darray<ds_uint> ReadUIntDarray(std::string file, std::string delim, bool header = true, bool with_row_name = false);
Darray<ds_uint> ReadUIntDarray(const char* file, char delim, bool header = true, bool with_row_name = false);
Darray<ds_double> ReadDoubleDarray(std::string file, std::string delim, bool header = true, bool with_row_name = false);
Darray<ds_double> ReadDoubleDarray(const char* file, char delim, bool header = true, bool with_row_name = false);

} // namespace

#endif // DSLICE_DARRAY_H_
