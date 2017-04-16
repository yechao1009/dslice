#include "../../include/dslice/base.h"
#include "../../include/dslice/covariate.h"
#include "../../include/dslice/dynamic_slicer.h"

int main(int argc, char** argv)
{
	std::string in_file;
	std::string delim = "\t";
	bool with_head = true;
	bool with_rowname = false;
	ds_double lambda = 1;
	ds_double alpha = 1;
	ds_double stat;

	in_file = "covariate.txt";
	dslice::Darray<ds_uint> variable = dslice::ReadUIntDarray(in_file, delim);
	in_file = "rank.txt";
	dslice::Darray<ds_uint> rank = dslice::ReadUIntDarray(in_file, delim);
	in_file = "response.txt";
	dslice::Darray<ds_double> response = dslice::ReadDoubleDarray(in_file, delim);

	std::string row_index;
	std::vector<ds_uint> slice_variable = variable(row_index, 0);
	std::vector<ds_uint> slice_rank = rank(row_index, 0);
	dslice::Covariate<ds_uint> covariate(slice_variable, true);
	dslice::Factor factor = covariate.GetFactor();
	dslice::DynamicSlicer slicer;
	dslice::SlicingReport report;

	std::cout << "Test standard Slicing" << std::endl;
	report = slicer.Slicing(factor, lambda);
	std::cout << "Report of dynamic slicing: " << std::endl;
	report.PrintReport();
	report = slicer.Slicing(factor, lambda, 0);
	std::cout << "Report of sqrt-partition dynamic slicing: " << std::endl;
	report.PrintReport();

	in_file = "quantile.txt";
	dslice::Darray<ds_double> quantile = dslice::ReadDoubleDarray(in_file, delim);
	std::vector<ds_double> slice_quantile = quantile(row_index, 0);
	std::cout << "Test OneSampleSlicing" << std::endl;
	stat = slicer.OneSampleSlicing(slice_quantile, lambda);
	std::cout << "Statistics of equal partition one sample slicing is: " << stat << std::endl;
	stat = slicer.OneSampleSlicing(slice_quantile, lambda, alpha);
	std::cout << "Statistics of one sample slicing is: " << stat << std::endl;

	in_file = "bf_data.txt";
	dslice::Darray<ds_uint> bf_data = dslice::ReadUIntDarray(in_file, delim);
	dslice::Covariate<ds_uint> covariate_x(bf_data(row_index, 0), false);
	dslice::Factor factor_x = covariate_x.GetFactor();
	dslice::Covariate<ds_uint> covariate_z(bf_data(row_index, 1), false);
	dslice::Factor factor_z = covariate_z.GetFactor();
	std::vector<ds_uint> z0(bf_data.Nrow(), 0);
	dslice::Covariate<ds_uint> covariate_z0(z0, false);
	dslice::Factor factor_z0 = covariate_z0.GetFactor();
	std::cout << "Test Bayesian Slicing" << std::endl;
	stat = slicer.BayesianSlicing(factor_x, lambda, alpha);
	std::cout << "Statistics of unconditional Bayesian slicing is: " << stat << std::endl;
	stat = slicer.BayesianSlicing(factor_x, factor_z0, lambda, alpha);
	std::cout << "Statistics of unconditional Bayesian slicing via conditional version is: " << stat << std::endl;
	stat = slicer.BayesianSlicing(factor_x, factor_z, lambda, alpha);
	std::cout << "Statistics of conditional Bayesian slicing is: " << stat << std::endl;
	std::cout << "Test Bayesian Slicing with sqrt-partition" << std::endl;
	stat = slicer.BayesianSlicing(factor_x, lambda, alpha, 0);
	std::cout << "Statistics of sqrt-partition unconditional Bayesian slicing is: " << stat << std::endl;
	stat = slicer.BayesianSlicing(factor_x, factor_z0, lambda, alpha, 0);
	std::cout << "Statistics of sqrt-partition unconditional Bayesian slicing via conditional version is: " << stat << std::endl;
	stat = slicer.BayesianSlicing(factor_x, factor_z, lambda, alpha, 0);
	std::cout << "Statistics of sqrt-partition conditional Bayesian slicing is: " << stat << std::endl;

}