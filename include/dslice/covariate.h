#ifndef DSLICE_COVARIATE_H_
#define DSLICE_COVARIATE_H_

#include <iostream>
#include <utility>
#include <vector>
#include <string>
#include <map>
#include "base.h"

namespace dslice
{

template <typename Object>
class Covariate
{
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

	Covariate (const Covariate<Object>& var) { *this = var; }

	void SetSize(ds_uint size) {
		size_ = size;
	}
	// fetch input sample
	std::vector<Object> GetSample() const {
		return sample_;
	}
	// fetch factor of input sample, begin with 0
	std::vector<ds_uint> GetFactor() const {
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
		return level_;
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
				fv[i] = factor_[index[i]];
			}
			sample_ = sv;
			factor_ = sv;
		}
	}

private:

	std::vector<Object> sample_;
	std::vector<ds_uint> factor_;
	std::map<Object, ds_uint> dict_;
	ds_uint size_;
	ds_uint level_;
	// generate factor for non-integer input sample (e.g., string)
	void GenerateFactor() {
		dict_.clear();
		factor_.resize(sample_.size());
		for (ds_uint i = 0; i != size_; ++i) {
			dict_.insert(std::pair<Object, ds_uint>(sample_[i], dict_.size()));
			factor_[i] = dict_[sample_[i]];
		}
		level_ = dict_.size();
	}
	// generate factor for integer sample
	void GenerateFactor(bool as_factor) {
		dict_.clear();
		factor_.resize(sample_.size());
		if (as_factor) {
			for (ds_uint i = 0; i != size_; ++i) {
				dict_.insert(std::pair<Object, ds_uint>(sample_[i], dict_.size()));
				factor_[i] = dict_[sample_[i]];
			}
			level_ = dict_.size();
		} else {
			for (ds_uint i = 0; i != size_; ++i) {
				factor_[i] = (ds_uint)sample_[i];
			}
			level_ = *max_element(factor_.begin(), factor_.end()) + 1;
		}
	}
};


} // namespace

#endif // DSLICE_COVARIATE_H_
