#ifndef DSLICE_SLICER_DYNAMIC_SLICER_H_
#define DSLICE_SLICER_DYNAMIC_SLICER_H_

#include <cmath>
#include <algorithm>
#include "slicing_report.h"

#define EPSILON 1e-6

namespace dslice {

class DynamicSlicer {
/*
 * For dynamic slicing
 * Parameter lambda in all slicing methods is the penalty for introducing an 
   additional slice, i.e., it is used to avoid making too many slices.
   lambda corresponds to the type I error under the scenario that the two 
   variables are independent.
   lambda should be greater than 0.
 * Parameter granularity is used in the scenario of large sample size.
   granularity specifies the minimum number of sample in each slice.
   granularity should be a positive integer.
 * Parameter alpha in OneSampleSlicing penalizes both the width and the number 
   of slices to avoid too many slices and degenerate slice (interval). 
   alpha should be greater than 1.
 * For Bayesian slicing
 * Parameter lambda corresponds to the probability that makes slice in each 
   possible position.
   lambda should be greater than 0.
 * Parameter alpha is hyper-parameter of the prior distribution of frequency in
   each slice. alpha should be greater than 0 and less equal than the factor 
   level of current investigated variable.
 */

public:
	/*
	 * \brief: dynamic slicing
	 * \param: factor, with value and level
	 * \param: penalty of introducing an additional slice
	 * \param: graunlarity of slicing
	 * \return: details of dynamic slicing
	 */
	dslice::SlicingReport Slicing(const dslice::Factor& variable, 
	                              ds_double lambda, 
	                              ds_uint granularity = 1);
	/*
	 * \brief: dynamic slicing
	 * \param: factor, with value and level
	 * \param: penalty of introducing an additional slice
	 * \param: rank of response
	 * \return: details of dynamic slicing
	 */
	dslice::SlicingReport Slicing(const dslice::Factor& variable, 
	                              ds_double lambda, 
	                              const std::vector<ds_uint>& rank);
	/*
	 * \brief: dynamic slicing for one-sample testing, equal bin version
	 * \param: quantile mapped to a given cdf
	 * \param: penalty of introducing an additional slice
	 * \param: rank of response, used if there are ties in ranking
	 * \return: statictics of dynamic slicing
	 */
	ds_double OneSampleSlicing(const std::vector<ds_double>& quantile, 
	                           ds_double lambda);
	/*
	 * \brief: dynamic slicing for one-sample testing
	 * \param: quantile mapped to a given cdf
	 * \param: penalty of introducing an additional slice
	 * \param: penalizes both the width and the number of slices
	 * \return: statictics of dynamic slicing
	 */
	ds_double OneSampleSlicing(const std::vector<ds_double>& quantile, 
	                           ds_double lambda, 
	                           ds_double alpha);
	/*
	 * \brief: Bayesian slicing
	 * \param: factor, with value and level
	 * \param: corresponds the probability that makes slice in each possible position
	 * \param: hyper-parameter of the prior distribution of frequency in each slice
	 * \param: graunlarity of slicing
	 * \return: statictics of Bayesian slicing
	 */
	ds_double BayesianSlicing(const dslice::Factor& variable, 
	                          ds_double lambda, 
	                          ds_double alpha, 
	                          ds_uint granularity = 1);
	/*
	 * \brief: Bayesian slicing
	 * \param: factor, with value and level
	 * \param: preselected factor, with value and level
	 * \param: corresponds the probability that makes slice in each possible position
	 * \param: hyper-parameter of the prior distribution of frequency in each slice
	 * \param: graunlarity of slicing
	 * \return: statictics of Bayesian slicing
	 */
	ds_double BayesianSlicing(const dslice::Factor& variable, 
	                          const dslice::Factor& selected, 
	                          ds_double lambda, 
	                          ds_double alpha, 
	                          ds_uint granularity = 1);
	
private:
	/*
	 * \brief: clustering adjacent samples with the sample value
	 * \param: factor, with value and level
	 * \return: cumsum count table of each value with level of variable
	 */
	dslice::Darray<ds_double> Clump(const dslice::Factor& variable);
	/*
	 * \brief: clustering adjacent samples with the sample value
	 * \param: factor, with value and level
	 * \param: graunlarity of clumping
	 * \return: cumsum count table of each value with level of variable
	 */
	dslice::Darray<ds_double> Clump(const dslice::Factor& variable, 
	                                ds_uint granularity);
	/*
	 * \brief: clustering adjacent samples with the sample value
	 * \param: factor, with value and level
	 * \param: rank information of continuous variable
	 * \return: cumsum count table of each value with level of variable
	 */
	dslice::Darray<ds_double> Clump(const dslice::Factor& variable, 
	                                const std::vector<ds_uint>& rank);
	/*
	 * \brief: wrapper of DynamicSlicing, return details of slicing
	 * \param: cumsum count table of each value
	 * \param: penalty of introducing an additional slice
	 * \return: details of dynamic slicing
	 */
	dslice::SlicingReport SlicingImplement(const dslice::Darray<ds_double>& count_table, 
	                                       ds_double penalty);
	/*
	 * \brief: return slicing strategy generate the max score
	 * \param: cumsum count table of each value
	 * \param: penalty of introducing an additional slice
	 * \return: dynamic programming strategy table and log likelihood of slicing
	 */
	dslice::DPTable DynamicSlicing(const dslice::Darray<ds_double>& count_table, 
	                               ds_double penalty);
	/*
	 * \brief: generate all possible slicing strategy via dynamic programming
	 * \param: dynamic programming strategy table
	 * \return: optimal slicing strategy
	 */
	std::vector<ds_uint> Backtracking(const std::vector<ds_uint>& global_incision);
	/*
	 * \brief: generate counts in each slice.
	 * \param: optimal slicing strategy
	 * \param: cumsum count table of each value
	 * \return: counts of each value in slice
	 */
	dslice::Darray<ds_uint> Summary(const std::vector<ds_uint>& optimal_incision, 
	                                const dslice::Darray<ds_double>& count_table);
	/*
	 * \brief: genreate cumsum count table of each value
	 * \param: factor, with value and level
	 * \return: cumsum count table of each value
	 */
	dslice::Darray<ds_uint> GenerateCountTable(const dslice::Factor& variable);
	/*
	 * \brief: genreate cumsum count table of each value
	 * \param: factor, with value and level
	 * \param: graunlarity of clumping
	 * \return: cumsum count table of each value
	 */
	dslice::Darray<ds_uint> GenerateCountTable(const dslice::Factor& variable, 
	                                           ds_uint granularity);
	/*
	 * \brief: genreate cumsum count table of each value
	 * \param: factor, with value and level
	 * \param: preselected factor, with value and level
	 * \return: cumsum count table of each value
	 */
	dslice::Darray<ds_uint> GenerateCountTable(const dslice::Factor& variable, 
	                                           const dslice::Factor& selected);
	/*
	 * \brief: genreate cumsum count table of each value
	 * \param: factor, with value and level
	 * \param: preselected factor, with value and level
	 * \param: graunlarity of clumping
	 * \return: cumsum count table of each value
	 */
	dslice::Darray<ds_uint> GenerateCountTable(const dslice::Factor& variable, 
	                                           const dslice::Factor& selected, 
	                                           ds_uint granularity);
	/*
	 * \brief: calculate Bayes factor for unconditional independent
	 * \param: cumsum count table of each value with level of current investigated variable
	 * \param: penalty of slicing
	 * \param: parameter of prior probability of slicing
	 * \return: Bayes factor
	 */
	ds_double CalculateBayesFactor(const dslice::Darray<ds_uint>& count_table, 
	                               ds_double lambda, 
	                               ds_double alpha);
	/*
	 * \brief: calculate Bayes factor for conditional independent
	 * \param: cumsum count table of each value with level = levlel_1 * levlel_2
	 * \param: penalty of slicing
	 * \param: parameter of prior probability of slicing
	 * \param: factor level of current investigated variable
	 * \return: Bayes factor
	 */
	ds_double CalculateBayesFactor(const dslice::Darray<ds_uint>& count_table, 
	                               ds_double lambda, 
	                               ds_double alpha, 
	                               ds_uint variable_level);

};

} // namespace dslice

#endif // DSLICE_SLICER_DYNAMIC_SLICER_H_
