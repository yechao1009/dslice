#ifndef DSLICE_SLICER_SLICING_REPORT_H_
#define DSLICE_SLICER_SLICING_REPORT_H_

#include "darray.h"

namespace dslice {

class SlicingReport {

public:
	SlicingReport() {};
	SlicingReport(const SlicingReport& report) { *this = report; }

	// fetch statistic of dynamic slicing
	ds_double GetStatistic() const { return statistic_; };
	// fetch the number of slice
	ds_uint GetSliceCount() const { return slice_count_; };
	// fetch slicing incision
	std::vector<ds_uint> GetIncision() const { return incision_; };
	// fetch sample counts in each slice
	dslice::Darray<ds_uint> GetDetail() const { return slices_; };

	void SetReport(ds_double statistic, std::vector<ds_uint> incision, dslice::Darray<ds_uint>& slices) {
		statistic_ = statistic;
		slice_count_ = slices.Nrow();
		incision_.assign(incision.begin(), incision.end());
		slices_ = slices;
	}

	void PrintReport() {
		std::cout << "Dynamic slicing statistic value: " << statistic_ << std::endl;
		std::cout << "The number of slices: " << slice_count_ << std::endl;
		std::cout << "The incisions: ";
		for (ds_uint i = 0; i != incision_.size(); ++i) {
			std::cout << incision_[i] << " ";
		}
		std::cout << std::endl;
		std::vector<std::string> rownames = slices_.RowName();
		std::vector<std::string> colnames = slices_.ColName();
		std::cout << "Slice detail: " << std::endl;
		std::cout << "dslice";
		for (int j = 0; j != colnames.size(); ++j) {
			std::cout << '\t' << colnames[j];
		}
		std::cout << std::endl;
		for (ds_uint i = 0; i != slices_.Nrow(); ++i) {
			std::cout << rownames[i];
			for (ds_uint j = 0; j != slices_.Ncol(); ++j) {
				std::cout << '\t' << slices_(i, j);
			}
			std::cout << std::endl;
		}
	}

private:
	ds_double statistic_;
	ds_uint slice_count_;
	std::vector<ds_uint> incision_;
	dslice::Darray<ds_uint> slices_;
	
};
	
} // namespace dslice

#endif // DSLICE_SLICER_SLICING_REPORT_H_
