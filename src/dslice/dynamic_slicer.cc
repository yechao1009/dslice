#include "../../include/dslice/dynamic_slicer.h"

namespace dslice {

dslice::SlicingReport DynamicSlicer::Slicing(const dslice::Factor& variable, ds_double lambda, ds_uint granularity) {
	dslice::Darray<ds_double> count_table;
	if (granularity == 1) {
		count_table = Clump(variable);
	} else {
		count_table = Clump(variable, granularity);
	}
	ds_uint sample_size = variable.factor.size();
	ds_double penalty = -lambda * log((ds_double)sample_size);
	dslice::SlicingReport slicing_report = SlicingImplement(count_table, penalty);
	return slicing_report;
}

dslice::SlicingReport DynamicSlicer::Slicing(const dslice::Factor& variable, ds_double lambda, const std::vector<ds_uint>& rank) {
	dslice::Darray<ds_double> count_table = Clump(variable, rank);
	ds_uint sample_size = variable.factor.size();
	ds_double penalty = -lambda * log((ds_double)sample_size);
	dslice::SlicingReport slicing_report = SlicingImplement(count_table, penalty);
	return slicing_report;
}

dslice::Darray<ds_double> DynamicSlicer::Clump(const dslice::Factor& variable) {
	ds_uint sample_size = variable.factor.size();
	dslice::Darray<ds_double> clump_array(sample_size + 1, variable.level);
	ds_uint index = 1;
	ds_uint clump_count = 1;
	ds_uint clump_num = 1;
	while (index < sample_size) {
		if (variable.factor[index] - variable.factor[index - 1] != 0) {
			clump_array.SetValue(clump_num, variable.factor[index - 1], clump_count);
			clump_count = 1;
			clump_num++;
		} else {
			clump_count++;
		}
		index++;
	}
	clump_array.SetValue(clump_num, variable.factor[sample_size - 1], clump_count);
	clump_array.Resize(clump_num + 1, variable.level);
	clump_array.SetValue(clump_array.CumSum(2));
	return clump_array;
}

dslice::Darray<ds_double> DynamicSlicer::Clump(const dslice::Factor& variable, ds_uint granularity) {
	ds_uint sample_size = variable.factor.size();
	if (granularity == 0) {
		granularity = (ds_uint)sqrt((ds_double)sample_size);
	}
	std::vector<ds_uint> candidate(sample_size / granularity + 2, 0);
	ds_uint clump_num = 0;
	ds_uint index = granularity;
	while (index < sample_size) {
		if (variable.factor[index] - variable.factor[index - 1] != 0) {
			candidate[++clump_num] = index;
			index += granularity;
		}else{
			index++;
		}
	}
	candidate[++clump_num] = sample_size;
	dslice::Darray<ds_double> clump_array(clump_num + 1, variable.level);
	for (ds_uint i = 1; i != clump_num + 1; ++i) {
		for (ds_uint j = 0; j != variable.level; ++j) {
			clump_array.SetValue(i, j, clump_array(i - 1, j));
		}
		for (ds_uint j = candidate[i - 1]; j < candidate[i]; ++j) {
			clump_array.SetValue(i, variable.factor[j], clump_array(i, variable.factor[j]) + 1);
		}
	}
	clump_array.Resize(clump_num+1, variable.level);
	return clump_array;
}

dslice::Darray<ds_double> DynamicSlicer::Clump(const dslice::Factor& variable, const std::vector<ds_uint>& rank) {
	ds_uint sample_size = variable.factor.size();
	dslice::Darray<ds_double> clump_array(sample_size + 1, variable.level);
	ds_uint index = 1;
	ds_uint clump_count = 1;
	ds_uint clump_num = 1;
	ds_double temp_value;
	while (index < sample_size) {
		if ((variable.factor[index] - variable.factor[index - 1] != 0) && (rank[index] - rank[index - 1] != 0)) {
			temp_value = clump_array(clump_num, variable.factor[index - 1]) + 1;
			clump_array.SetValue(clump_num, variable.factor[index - 1], temp_value);
			clump_count = 1;
			clump_num++;
		} else {
			temp_value = clump_array(clump_num, variable.factor[index - 1]) + 1;
			clump_array.SetValue(clump_num, variable.factor[index - 1], temp_value);
		}
		index++;
	}
	temp_value = clump_array(clump_num, variable.factor[index - 1]) + 1;
	clump_array.SetValue(clump_num, variable.factor[index - 1], temp_value);
	clump_array.Resize(clump_num+1, variable.level);
	clump_array.SetValue(clump_array.CumSum(2));
	clump_array.Print();
	return clump_array;
}

dslice::SlicingReport DynamicSlicer::SlicingImplement(const dslice::Darray<ds_double>& count_table, ds_double penalty) {
	DPTable programming_result = DynamicSlicing(count_table, penalty);
	std::vector<ds_uint> optimal_incision = Backtracking(programming_result.strategy);
	dslice::Darray<ds_uint> slices = Summary(optimal_incision, count_table);
	dslice::SlicingReport slicing_report;
	slicing_report.SetReport(programming_result.score, optimal_incision, slices);
	return slicing_report;
}

dslice::DPTable DynamicSlicer::DynamicSlicing(const dslice::Darray<ds_double>& count_table, ds_double penalty) {
	ds_uint clump_num = count_table.Nrow() - 1;
	ds_uint level = count_table.Ncol();
	DPTable dp_table;
	dp_table.strategy.resize(clump_num + 1);
	std::vector<ds_double> score(clump_num + 1, 0);
	std::vector<ds_uint> counts(level);
	ds_uint incision;
	ds_double local_sum, candidate_score, slice_score;
	// dynamic programming
	for (ds_uint i = 1; i != clump_num + 1; ++i) {
		// j = 0
		slice_score = penalty + score[0];
		local_sum = 0;
		for (ds_uint u = 0; u != level; ++u) {
			counts[u] = count_table(i, u);
			local_sum += counts[u];
			if(counts[u] > EPSILON){
				slice_score += counts[u] * log(counts[u]);
			}
		}
		slice_score -= local_sum * log(local_sum);
		incision = 0;
		// j = 0 end
		for (ds_uint j = 1; j != i; ++j) {
			candidate_score = penalty + score[j];
			local_sum = 0;
			for (ds_uint u = 0; u != level; ++u) {
				counts[u] = count_table(i, u) - count_table(j, u);
				local_sum += counts[u];
				if(counts[u] > EPSILON){
					candidate_score += counts[u] * log(counts[u]);
				}
			}
			candidate_score -= local_sum * log(local_sum);
			if (slice_score < candidate_score) {
				slice_score = candidate_score;
				incision = j;
			}
		}
		score[i] = slice_score;
		dp_table.strategy[i] = incision;
	}
	// substract null log-likelihood (assume one slice, i.e., no cut)
	ds_double max_loglikelihood = score[clump_num] - penalty;
	ds_uint sample_size = 0;
	for(ds_uint i = 0; i != level; ++i) {
		if (count_table(clump_num, i) > EPSILON){
			max_loglikelihood -= count_table(clump_num, i) * log(count_table(clump_num, i));
			sample_size += count_table(clump_num, i);
		}
	}
	max_loglikelihood += sample_size * log(sample_size);
	dp_table.score = max_loglikelihood;
	return dp_table;
}

std::vector<ds_uint> DynamicSlicer::Backtracking(const std::vector<ds_uint>& global_incision) {
	std::vector<ds_uint> optimal_incision;
	ds_uint index = global_incision.size() - 1;
	while (index > 0) {
		optimal_incision.push_back(index);
		index = global_incision[index];
	}
	optimal_incision.push_back(0);
	reverse(optimal_incision.begin(), optimal_incision.end());
	return optimal_incision;
}

dslice::Darray<ds_uint> DynamicSlicer::Summary(const std::vector<ds_uint>& optimal_incision, 
	const dslice::Darray<ds_double>& count_table) {
	ds_uint level = count_table.Ncol();
	ds_uint slice_num = optimal_incision.size() - 1;
	dslice::Darray<ds_uint> slices(slice_num, level + 1);
	std::vector<std::string> rownames(slice_num);
	std::vector<std::string> colnames(level + 1);
	std::string prefix = "slice";
	for (ds_uint i = 0; i != slice_num; ++i) {
		rownames[i] = prefix + std::to_string(i + 1);
	}
	for (ds_uint j = 0; j != level; ++j) {
		colnames[j] = std::to_string(j);
	}
	colnames[level] = "total";
	slices.SetDimName(rownames, 1);
	slices.SetDimName(colnames, 2);

	ds_double category_count = 0;
	for (ds_uint i = 0; i != slice_num; ++i) {
		for (ds_uint j = 0; j != level; ++j) {
			category_count = count_table(optimal_incision[i + 1], j) - count_table(optimal_incision[i], j);
			slices.SetValue(i, j, category_count);
			slices.SetValue(i, level, category_count + slices(i, level));
		}
	}
	return slices;
}

ds_double DynamicSlicer::BayesianSlicing(const dslice::Factor& variable, ds_double lambda, ds_double alpha, ds_uint granularity) {
	dslice::Darray<ds_uint> count_table;
	if (granularity == 1) {
		count_table = GenerateCountTable(variable);
	} else {
		count_table = GenerateCountTable(variable, granularity);
	}
	ds_double dsbf = CalculateBayesFactor(count_table, lambda, alpha);
	return dsbf;
}

ds_double DynamicSlicer::BayesianSlicing(const dslice::Factor& variable, const dslice::Factor& selected, ds_double lambda, ds_double alpha, ds_uint granularity) {
	dslice::Darray<ds_uint> count_table;
	if (granularity == 1) {
		count_table = GenerateCountTable(variable, selected);
	} else {
		count_table = GenerateCountTable(variable, selected, granularity);
	}
	ds_double dsbf = CalculateBayesFactor(count_table, lambda, alpha, variable.level);
	return dsbf;
}

dslice::Darray<ds_uint> DynamicSlicer::GenerateCountTable(const dslice::Factor& variable) {
	ds_uint sample_size = variable.factor.size();
	ds_double temp_value;
	dslice::Darray<ds_uint> count_table(sample_size + 1, variable.level);
	for (ds_uint i = 1; i != sample_size + 1; ++i) {
		for (ds_uint j = 0; j != variable.level; ++j) {
			count_table.SetValue(i, j, count_table(i - 1, j));
		}
		temp_value = count_table(i, variable.factor[i - 1]) + 1;
		count_table.SetValue(i, variable.factor[i - 1], temp_value);
	}
	return count_table;
}

dslice::Darray<ds_uint> DynamicSlicer::GenerateCountTable(const dslice::Factor& variable, ds_uint granularity) {
	ds_uint sample_size = variable.factor.size();
	if (granularity == 0) {
		granularity = (ds_uint)sqrt((ds_double)sample_size);
	}
	ds_uint group_num = sample_size / granularity;
	ds_uint tail = sample_size % granularity;
	std::vector<ds_uint> candidate(group_num + 1, 0);
	for (ds_uint i = 1; i != group_num + 1; ++i) {
		candidate[i] = candidate[i - 1] + granularity;
	}
	if (tail != 0) {
		for (ds_int i = 1; i != tail + 1; ++i){
			candidate[i] += i;
		}
		for (ds_int i = tail + 1; i != group_num + 1; ++i) {
			candidate[i] += tail;
		}
	}
	ds_double temp_value;
	dslice::Darray<ds_uint> clump_array(group_num + 1, variable.level);
	for (ds_uint i = 1; i != group_num + 1; ++i) {
		for (ds_uint j = 0; j != variable.level; ++j) {
			clump_array.SetValue(i, j, clump_array(i - 1, j));
		}
		for (ds_uint j = candidate[i - 1]; j != candidate[i]; ++j) {
			temp_value = clump_array(i, variable.factor[j]) + 1;
			clump_array.SetValue(i, variable.factor[j], temp_value);
		}
	}
	return clump_array;
}

dslice::Darray<ds_uint> DynamicSlicer::GenerateCountTable(const dslice::Factor& variable, const dslice::Factor& selected) {
	ds_uint sample_size = variable.factor.size();
	ds_uint level = variable.level * selected.level;
	ds_int col_index;
	ds_double temp_value;
	dslice::Darray<ds_uint> count_table(sample_size + 1, level);
	for (ds_uint i = 1; i != sample_size + 1; ++i) {
		for (ds_uint j = 0; j != level; ++j) {
			count_table.SetValue(i, j, count_table(i - 1, j));
		}
		col_index = selected.factor[i - 1] * variable.level + variable.factor[i - 1];
		temp_value = count_table(i, col_index) + 1;
		count_table.SetValue(i, col_index, temp_value);
	}
	return count_table;
}

dslice::Darray<ds_uint> DynamicSlicer::GenerateCountTable(const dslice::Factor& variable, const dslice::Factor& selected, ds_uint granularity) {
	ds_uint sample_size = variable.factor.size();
	if (granularity == 0) {
		granularity = (ds_uint)sqrt((ds_double)sample_size);
	}
	ds_uint group_num = sample_size / granularity;
	ds_uint tail = sample_size % granularity;
	std::vector<ds_uint> candidate(group_num + 1, 0);
	for (ds_uint i = 1; i != group_num + 1; ++i) {
		candidate[i] = candidate[i - 1] + granularity;
	}
	if (tail != 0) {
		for (ds_uint i = 1; i != tail + 1; ++i){
			candidate[i] += i;
		}
		for (ds_uint i = tail + 1; i != group_num + 1; ++i) {
			candidate[i] += tail;
		}
	}
	ds_uint level = variable.level * selected.level;
	ds_uint col_index;
	ds_double temp_value;
	dslice::Darray<ds_uint> clump_array(group_num + 1, level);
	for (ds_uint i = 1; i != group_num + 1; ++i) {
		for (ds_uint j = 0; j != level; ++j) {
			clump_array.SetValue(i, j, clump_array(i - 1, j));
		}
		for (ds_uint j = candidate[i - 1]; j != candidate[i]; ++j) {
			col_index = selected.factor[j] * variable.level + variable.factor[j];
			temp_value = clump_array(i, col_index) + 1;
			clump_array.SetValue(i, col_index, temp_value);
		}
	}
	return clump_array;
}

ds_double DynamicSlicer::CalculateBayesFactor(const dslice::Darray<ds_uint>& count_table, ds_double lambda, ds_double alpha) {
	ds_uint length = count_table.Nrow() - 1;
	ds_uint level = count_table.Ncol();
	ds_double beta = alpha / level;
	ds_uint sample_size = 0;
	for (ds_uint j = 0; j < count_table.Ncol(); ++j) {
		sample_size += count_table(length, j);
	}
	std::vector<ds_double> alpha_table(sample_size + 1, 0);
	std::vector<ds_double> beta_table(sample_size + 1, 0);
	for (ds_uint i = 1; i != sample_size + 1; ++i) {
		alpha_table[i] = alpha_table[i - 1] + log(alpha + i - 1);
		beta_table[i] = beta_table[i - 1] + log(beta + i - 1);
	}
	// log Probability that null hypothesis holds
	ds_uint counts;
	std::vector<ds_double> log_pr_null(length + 1, 0);
	for (ds_uint i = 1; i != length + 1; ++i) {
		counts = 0;
		for (ds_uint v = 0; v != level; ++v){
			counts += count_table(i, v);
			log_pr_null[i] += beta_table[count_table(i, v)];
		}
		log_pr_null[i] -= alpha_table[counts];
	}
	// calculate Bayes factor
	ds_double log_delta_ji, temp_value;
	ds_uint count_diff;
	std::vector<ds_double> bayes_factor(length + 1, 0);
	ds_double prior = 1.0 / (1.0 + exp(lambda * log((ds_double)sample_size)));
	bayes_factor[0] = 1 / prior;
	for (ds_uint i = 1; i != length + 1; ++i) {
		for (ds_uint j = 0; j != i; ++j) {
			log_delta_ji = 0;
			counts = 0;
			for (ds_int v = 0; v != level; ++v) {
				count_diff = count_table(i, v) - count_table(j, v);
				counts += count_diff;
				log_delta_ji += beta_table[count_diff];
			}
			log_delta_ji += log(prior) - log(1 - prior) - alpha_table[counts];
			temp_value = (i - j) * log(1 - prior) + log_pr_null[j] - log_pr_null[i] + log_delta_ji;
			bayes_factor[i] += exp(temp_value) * bayes_factor[j];
		}
	}
	ds_double dsbf = bayes_factor[length];
	return dsbf;
}

ds_double DynamicSlicer::CalculateBayesFactor(const dslice::Darray<ds_uint>& count_table, ds_double lambda, ds_double alpha, ds_uint variable_level) {
	ds_uint length = count_table.Nrow() - 1;
	ds_uint level = count_table.Ncol();
	ds_uint selected_level = level / variable_level;
	ds_double beta = alpha / variable_level;
	ds_uint sample_size = 0;
	for (ds_uint j = 0; j < count_table.Ncol(); ++j) {
		sample_size += count_table(length, j);
	}
	std::vector<ds_double> alpha_table(sample_size + 1, 0);
	std::vector<ds_double> beta_table(sample_size + 1, 0);
	for (ds_uint i = 1; i != sample_size + 1; ++i) {
		alpha_table[i] = alpha_table[i - 1] + log(alpha + i - 1);
		beta_table[i] = beta_table[i - 1] + log(beta + i - 1);
	}
	ds_uint col_index;
	// log Probability that null hypothesis holds
	std::vector<ds_uint> selected_count(selected_level, 0);
	std::vector<ds_double> log_pr_null(length + 1, 0);
	for (ds_uint i = 1; i != length + 1; ++i) {
		for (ds_uint u = 0; u != selected_level; ++u) {
			selected_count[u] = 0;
			for (ds_uint v = 0; v != variable_level; ++v){
				col_index = u * variable_level + v;
				selected_count[u] += count_table(i, col_index);
				log_pr_null[i] += beta_table[count_table(i, col_index)];
			}
			log_pr_null[i] -= alpha_table[selected_count[u]];
		}
	}
	// calculate Bayes factor
	ds_double log_delta_ji, temp_value;
	ds_uint count_diff;
	std::vector<ds_double> bayes_factor(length + 1, 0);
	ds_double prior = 1.0 / (1.0 + exp(lambda * log((ds_double)sample_size)));
	bayes_factor[0] = 1 / prior;
	for (ds_uint i = 1; i != length + 1; ++i) {
		for (ds_uint j = 0; j != i; ++j) {
			log_delta_ji = 0;
			for (ds_int u = 0; u != selected_level; ++u) {
				selected_count[u] = 0;
				for (ds_int v = 0; v != variable_level; ++v) {
					col_index = u * variable_level + v;
					count_diff = count_table(i, col_index) - count_table(j, col_index);
					selected_count[u] += count_diff;
					log_delta_ji += beta_table[count_diff];
				}
				log_delta_ji -= alpha_table[selected_count[u]];
			}
			log_delta_ji += log(prior) - log(1 - prior);
			temp_value = (i - j) * log(1 - prior) + log_pr_null[j] - log_pr_null[i] + log_delta_ji;
			bayes_factor[i] += exp(temp_value) * bayes_factor[j];
		}
	}
	ds_double dsbf = bayes_factor[length];
	return dsbf;
}

ds_double DynamicSlicer::OneSampleSlicing(const std::vector<ds_double>& quantile, ds_double lambda, ds_double alpha) {
	ds_uint length = quantile.size();
	ds_double penalty = -lambda * log((ds_double)length);
	ds_double unit = 1.0 / length;
	std::vector<ds_double> left_score(length + 1, 0);
	std::vector<ds_double> right_score(length + 1, 0);
	std::vector<ds_uint> index(2 * length + 2, -1);

	ds_uint incision;
	ds_double candidate_score, slice_score;
	ds_double order_diff, quantile_diff;
	ds_double left_cut_value, right_cut_value;
	// i = 0, log likelihood = 0 for the left side
	left_cut_value = alpha * log(quantile[0]);
	right_cut_value = log(unit) + (alpha - 1) * log(quantile[0]); // log(unit / quantile[0]) + alpha * log(quantile[0])
	left_score[0] = left_cut_value + penalty;
	right_score[0] = right_cut_value + penalty;
	index[1] = 0;
	index[2] = 0;
	// i = 0 end
	for(ds_uint i = 1; i != length; ++i) {
		// on left side of i
		left_cut_value = i * log(i * unit / quantile[i]);
		slice_score = left_cut_value + penalty + alpha * log(quantile[i]);
		incision = 0;
		for (ds_uint j = 0; j != i - 1; ++j) {
			order_diff = i - j;
			quantile_diff = quantile[i] - quantile[j];
			left_cut_value = left_score[j] + order_diff * log(order_diff * unit / quantile_diff);
			right_cut_value = right_score[j] + (order_diff - 1) * log((order_diff - 1) * unit / quantile_diff);
			candidate_score = std::max(left_cut_value, right_cut_value) + penalty + alpha * log(quantile_diff);
			if (slice_score < candidate_score) {
				slice_score = candidate_score;
				incision = j;
			}
		}
		// j = i - 1
		order_diff = 1;
		quantile_diff = quantile[i] - quantile[i - 1];
		left_cut_value = left_score[i - 1] + log(unit / quantile_diff);
		right_cut_value = right_score[i - 1];
		candidate_score = std::max(left_cut_value, right_cut_value) + penalty + alpha * log(quantile_diff);
		if (slice_score < candidate_score) {
			slice_score = candidate_score;
			incision = i - 1;
		}
		// j = i - 1 end
		left_score[i] = slice_score;
		index[2 * i + 1] = incision;
		// on right side of i
		right_cut_value = (i + 1) * log((i + 1) * unit / quantile[i]);
		slice_score = right_cut_value + penalty + alpha * log(quantile[i]);
		incision = 0;
		for (ds_uint j = 0; j != i; ++j) {
			order_diff = i - j;
			quantile_diff = quantile[i] - quantile[j];
			left_cut_value = left_score[j] + (order_diff + 1) * log((order_diff + 1) * unit / quantile_diff);
			right_cut_value = right_score[j] + order_diff * log(order_diff * unit / quantile_diff);
			candidate_score = std::max(left_cut_value, right_cut_value) + penalty + alpha * log(quantile_diff);
			if(slice_score < candidate_score){
				slice_score = candidate_score;
				incision = j;
			}
		}
		right_score[i] = slice_score;
		index[2 * i + 2] = incision;
	}
	// i = length, i.e. for quantile_{n+1} = 1.0
	// on left side of len, no right side
	slice_score = penalty; // n*log(n/n /1) + penalty + alpha * log(1)
	incision = 0;
	for (ds_uint j = 0; j != length - 1; ++j) {
		order_diff = length - j;
		quantile_diff = 1.0 - quantile[j];
		left_cut_value = left_score[j] + order_diff * log(order_diff * unit / quantile_diff);
		right_cut_value = right_score[j] + (order_diff - 1) * log((order_diff - 1) * unit / quantile_diff);
		candidate_score = std::max(left_cut_value, right_cut_value) + penalty + alpha * log(quantile_diff);
		if (slice_score < candidate_score) {
			slice_score = candidate_score;
			incision = j;
		}
	}
	// j = length - 1
	order_diff = 1;
	quantile_diff = 1.0 - quantile[length - 1];
	left_cut_value = left_score[length - 1] + log(unit / quantile_diff);
	right_cut_value = right_score[length - 1];
	candidate_score = std::max(left_cut_value, right_cut_value) + penalty + alpha * log(quantile_diff);
	if (slice_score < candidate_score) {
		slice_score = candidate_score;
		incision = length - 1;
	}
	// j = length - 1 end
	left_score[length] = slice_score;
	index[2 * length + 1] = incision;
	// i = length end
	double max_likelihood = left_score[length] - penalty;
	return max_likelihood;
}

ds_double DynamicSlicer::OneSampleSlicing(const std::vector<ds_double>& quantile, ds_double lambda) {
	ds_uint length = quantile.size();
	ds_double penalty = -lambda * log((ds_double)length);
	ds_double unit = 1.0 / length;
	std::vector<ds_double> cum_uint(length + 1, 0); // cumulated unit counts
	ds_uint current_index = 1;
	ds_double current_cum = unit;
	for (ds_uint i = 0; i < length; ++i) {
		while (quantile[i] > current_cum) {
			cum_uint[current_index++] = i;
			current_cum += unit;
		}
	}
	for(; current_index < length + 1; ++current_index) {
		cum_uint[current_index] = length;
	}
	std::vector<ds_double> score(length + 1, 0);
	std::vector<ds_double> index(length + 1, -1);
	ds_uint incision;
	ds_double candidate_score, slice_score;
	ds_double order_diff, quantile_diff;
	for (ds_uint i = 1; i != length + 1; ++i) {
		slice_score = penalty + score[0];
		quantile_diff = i;
		order_diff = cum_uint[i] - cum_uint[0];
		if (order_diff > EPSILON) {
			slice_score += order_diff * log(order_diff / quantile_diff);
		}
		incision = 0;
		for (ds_uint j = 1; j != i; ++j) {
			candidate_score = penalty + score[j];
			quantile_diff = i - j;
			order_diff = cum_uint[i] - cum_uint[j];
			if (order_diff > EPSILON) {
				candidate_score += order_diff * log(order_diff / quantile_diff);
			}
			if(slice_score < candidate_score) {
				slice_score = candidate_score;
				incision = j;
			}
		}
		score[i] = slice_score;
		index[i] = incision;
	}
	ds_double max_likelihood = score[length] - penalty;
	quantile_diff = length;
	order_diff = cum_uint[length];
	max_likelihood -= order_diff * log(order_diff / quantile_diff);
	return max_likelihood; // maximized log-likelihood
}

} // namespace dslice