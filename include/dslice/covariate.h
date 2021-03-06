#ifndef DSLICE_COVARIATE_H_
#define DSLICE_COVARIATE_H_

#include <iostream>
#include <utility>
#include <map>
#include <algorithm>
#include "base.h"

namespace dslice
{

template <typename Object>
class Covariate {

public:
	Covariate(std::vector<Object> sample) {
		sample_ = sample;
		size_ = sample.size();
		GenerateFactor();
	}

	Covariate(std::vector<int> sample, bool as_factor = true) {
		sample_ = sample;
		size_ = sample.size();
		GenerateFactor(as_factor);
	}

	Covariate(std::vector<unsigned> sample, bool as_factor = true) {
		sample_ = sample;
		size_ = sample.size();
		GenerateFactor(as_factor);
	}

	Covariate(const Covariate<Object>& var) { *this = var; }

	void SetSize(ds_uint size) {
		size_ = size;
	}
	// fetch input sample
	std::vector<Object> GetSample() const {
		return sample_;
	}
	// fetch factor of input sample, begin with 0
	dslice::Factor GetFactor() const {
		return factor_;
	}
	// fetch sample to factor map
	std::map<Object, ds_uint> GetDict() const {
		return dict_;
	}
	// fetch sample size
	ds_uint GetSize() const {
		return size_;
	}
	// fecth factor level
	ds_uint GetLevel() const {
		return factor_.level;
	}
	// reorder covariate according to given index
	void Reorder(std::vector<ds_uint> index) {
		if (index.size() != size_) {
			return;
		} else {
			std::vector<Object> sv;
			std::vector<ds_uint> fv;
			for (ds_uint i = 0; i != size_; ++i) {
				sv[i] = sample_[index[i]];
				fv[i] = factor_.factor[index[i]];
			}
			sample_ = sv;
			factor_.factor = sv;
		}
	}

private:
	std::vector<Object> sample_;
	Factor factor_;
	std::map<Object, ds_uint> dict_;
	ds_uint size_;

	// generate factor for non-integer input sample (e.g., string)
	void GenerateFactor() {
		dict_.clear();
		factor_.factor.resize(sample_.size());
		for (ds_uint i = 0; i != size_; ++i) {
			dict_.insert(std::pair<Object, ds_uint>(sample_[i], dict_.size()));
			factor_.factor[i] = dict_[sample_[i]];
		}
		factor_.level = dict_.size();
	}
	// generate factor for integer sample
	void GenerateFactor(bool as_factor) {
		dict_.clear();
		factor_.factor.resize(sample_.size());
		if (as_factor) {
			for (ds_uint i = 0; i != size_; ++i) {
				dict_.insert(std::pair<Object, ds_uint>(sample_[i], dict_.size()));
				factor_.factor[i] = dict_[sample_[i]];
			}
			factor_.level = dict_.size();
		} else {
			for (ds_uint i = 0; i != size_; ++i) {
				factor_.factor[i] = (ds_uint)sample_[i];
			}
			factor_.level = *max_element(factor_.factor.begin(), factor_.factor.end()) + 1;
		}
	}
};


} // namespace

#endif // DSLICE_COVARIATE_H_
