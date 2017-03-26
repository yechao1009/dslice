#ifndef DSLICE_DARRAY_H_
#define DSLICE_DARRAY_H_

#include <iostream>
#include <vector>
#include <string>
#include "base.h"

namespace dslice
{

// Double dimensional array
template <typename Object>
class Darray {

public:
	explicit Darray() : darray_(0) {}

	Darray(ds_uint rows, ds_uint cols) : darray_(rows) {
		for (ds_uint i = 0; i != rows; ++i) {
			darray_[i].resize(cols);
		}
	}
	Darray(const Darray<Object>& darray) { *this = darray; }
	ds_uint Nrow() const { return darray_.size(); }
	ds_uint Ncol() const { return Nrow() ? (darray_[0].size()) : 0; }
	bool Empty() const { return Nrow() == 0; }
	const std::vector<Object>& operator[](ds_uint row) const { return darray_[row]; };
	std::vector<Object>& operator[](ds_uint row) { return darray_[row]; };
	const Object& operator()(ds_uint row, ds_uint col) const { return darray_[row][col]; };
	Object& operator()(ds_uint row, ds_uint col) { return darray_[row][col]; };

	void SetValue(ds_uint row_id, ds_uint col_id, Object value) {
		if (row_id < darray_.size() && col_id < darray_[0].size()) {
			darray_[row_id][col_id] = value;
		}
	}

	bool SetDimName(const std::vector<std::string>& dim_name, int type) {
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

	std::vector<std::string> RowName() const { return row_name_; };
	std::vector<std::string> ColName() const { return col_name_; };

	void Resize(ds_uint rows, ds_uint cols);
	bool PushBack(const std::vector<Object>& vkt);
	std::vector<Object> Sum(int direction) const;


protected:
	std::vector< std::vector<Object> > darray_;
	std::vector<std::string> row_name_;
	std::vector<std::string> col_name_;
};

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
		for (ds_uint i = nrow; i != rows; ++i) {
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

class IntDarray : public Darray<ds_int> {

public:
	IntDarray() : Darray<ds_int>() {};
	IntDarray(ds_uint rows, ds_uint cols) : Darray<ds_int>(rows, cols) {};
	IntDarray(const IntDarray& int_darray) { *this = int_darray; }

	const IntDarray& operator+=(const IntDarray& int_darray);
	const IntDarray& operator-=(const IntDarray& int_darray);

};

bool operator==(const IntDarray& lhs, const IntDarray& rhs);
bool operator!=(const IntDarray& lhs, const IntDarray& rhs);
const IntDarray operator+(const IntDarray& lhs, const IntDarray& rhs);
const IntDarray operator-(const IntDarray& lhs, const IntDarray& rhs);
const IntDarray ReadDarray(std::string file, std::string delim, bool header = true, bool with_row_name = false);
const IntDarray ReadDarray(std::istream& in_stream, char delim, bool header = true, bool with_row_name = false);
void WriteDarray(const IntDarray& int_darray, std::string file, std::string delim, bool col_name = true, bool row_name = false);
void WriteDarray(const IntDarray& int_darray, std::ostream& out_stream, std::string delim, bool colname = true, bool row_name = false);

} // namespace

#endif // DSLICE_DARRAY_H_
