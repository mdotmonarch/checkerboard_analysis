# -*- coding: utf-8 -*-
# @author: max (max@mdotmonar.ch)

using HDF5
using Plots
using StatsPlots
using DSP
using HypothesisTests

include("processing.jl")
include("entropy.jl")

groups = ["A", "B", "C", "D", "E", "F", "G", "H"]

group_labels = Dict()
group_labels["A"] = "WT young"
group_labels["B"] = "WT adult"
group_labels["C"] = "5xFAD young"
group_labels["D"] = "5xFAD adult"
group_labels["E"] = "XBP1s young"
group_labels["F"] = "XBP1s adult"
group_labels["G"] = "Double young"
group_labels["H"] = "Double adult"
group_labels_plot = ["WT young" "WT adult" "5xFAD young" "5xFAD adult" "XBP1s young" "XBP1s adult" "Double young" "Double adult"]

grouped_datasets = Dict()
grouped_datasets["A"] = [
	# WT 3m male
	"MR-0474",
	"MR-0311",
	"MR-0309",
	"MR-0306",
	"MR-0299-t2",
	"MR-0299-t1",
	"MR-0298-t2",
	"MR-0298-t1",
	"MR-0296-t2",
	"MR-0296-t1",
	# WT 3m female
	"MR-0491",
	"MR-0303",
	"MR-0300-t2",
	"MR-0300-t1",
]
grouped_datasets["B"] = [
	# WT 6m male
	"MR-0282",
	"MR-0276",
	"MR-0273",
	"MR-0270",
	# WT 6m female
	"MR-0294",
	"MR-0289",
	"MR-0288-t2",
	"MR-0288-t1",
	"MR-0284",
	"MR-0283-t2",
	"MR-0283-t1",
]
grouped_datasets["C"] = [
	# 5xFAD 3m male
	"MR-0310",
	"MR-0305",
	# 5xFAD 3m female
	"MR-0313",
	"MR-0312",
	"MR-0307-t2",
	"MR-0307-t1",
	"MR-0304-t2",
	"MR-0304-t1",
	"MR-0302-t2",
	"MR-0302-t1",
	"MR-0301-t2",
	"MR-0301-t1",
	"MR-0297",
]
grouped_datasets["D"] = [
	# 5xFAD 6m male
	"MR-0448",
	"MR-0447",
	"MR-0446",
	"MR-0293",
	"MR-0292-t2",
	"MR-0292-t1",
	"MR-0280-t2",
	"MR-0280-t1",
	"MR-0278",
	"MR-0275",
	# 5xFAD 6m female
	"MR-0291",
	"MR-0290",
	"MR-0287",
	"MR-0285-t2",
	"MR-0285-t1",
	"MR-0274",
]
grouped_datasets["E"] = [
	# XBP1s 3m male
	"MR-0592_nd4",
	"MR-0591_nd4",
	"MR-0483",
	"MR-0460",
	"MR-0456",
	# XBP1s 3m female
	"MR-0621_nd4",
	"MR-0620_nd4",
	"MR-0593_nd4",
	"MR-0573_nd4",
]
grouped_datasets["F"] = [
	# XBP1s 6m male
	"MR-0599_nd4",
	"MR-0597_nd4",
	"MR-0596_nd4",
	"MR-0569_nd4",
	"MR-0554",
	"MR-0465-t2",
	"MR-0465-t1",
	# XBP1s 6m female
	"MR-0625_nd4",
	"MR-0624_nd4",
	"MR-0623_nd4",
	"MR-0622_nd4",
]
grouped_datasets["G"] = [
	# Double 3m male
	"MR-0575_nd4",
	"MR-0583_nd4",
	# Double 3m female
	"MR-0577_nd4",
	"MR-0579_nd4",
	"MR-0585_nd4",
	"MR-0586_nd4",
]
grouped_datasets["H"] = [
	# Double 6m male
	"MR-0630_nd4",
	"MR-0629-t2_nd4",
	"MR-0629-t1_nd4",
	"MR-0552",
	# Double 6m female
	"MR-0530",
	"MR-0548-t1",
	"MR-0548-t2",
	"MR-0568-t1_nd4",
	"MR-0568-t2_nd4",
	"MR-0582_nd4",
	"MR-0587_nd4",
	"MR-0588_nd4",
]

t_list = ["RCMSE"]
r_list = ["0.1", "0.2", "0.3", "0.4", "0.5"]

# load entropy data
entropy_data = Dict()

for d in datasets
	entropy_data[d] = Dict()
	for t in t_list
		entropy_data[d][t] = Dict()
		for r in r_list
			entropy_data[d][t][r] = Dict()
		end
	end
end

for dataset in datasets
	file = h5open("./entropy_data/$(dataset)_checkerboard_entropy.h5", "r")

	for type in t_list
		for r in r_list
			entropy = read(file["/$(type)/$(r)/mean_signal"])
			entropy_data[dataset][type][r]["nAUC"] = entropy["nAUC"]
			entropy_data[dataset][type][r]["LRS"] = entropy["LRS"]
		end
	end

	close(file)
end

# tests
for type in t_list
	for r in r_list
		local grouped_lrs = Dict()
		local grouped_nauc = Dict()
		for group in groups
			grouped_lrs[group] = Float64[]
			grouped_nauc[group] = Float64[]
		end
		for g in groups
			for d in grouped_datasets[g]
				push!(grouped_lrs[g], entropy_data[d][type][r]["LRS"])
				push!(grouped_nauc[g], entropy_data[d][type][r]["nAUC"])
			end
		end

		# kruskal-wallis test for LRS
		f = open("./analysis/$(type)_$(r)_LRS.txt", "w")
		println("$(type) $(r) LRS: Analysis #####################################################")
		write(f, "$(type) $(r) LRS: Analysis #####################################################\n")
		kruskal_lrs = KruskalWallisTest(grouped_lrs["A"], grouped_lrs["B"], grouped_lrs["C"], grouped_lrs["D"], grouped_lrs["E"], grouped_lrs["F"], grouped_lrs["G"], grouped_lrs["H"])
		println("Kruskal-Wallis test p-value: $(pvalue(kruskal_lrs))")
		write(f, "Kruskal-Wallis test p-value: $(pvalue(kruskal_lrs))\n\n")
		println("")
		if pvalue(kruskal_lrs) < 0.05
			# perform mann-whitney u test for each pair
			println("Significant differences (Mann-Whitney U test):")
			write(f, "Significant differences (Mann-Whitney U test):\n")
			for i in 1:8
				for j in i+1:8
					mannwhitney = MannWhitneyUTest(grouped_lrs[groups[i]], grouped_lrs[groups[j]])
					if pvalue(mannwhitney) < 0.05
						println("$(group_labels[groups[i]]) vs $(group_labels[groups[j]]):\t\t\t$(pvalue(mannwhitney))")
						write(f, "$(group_labels[groups[i]]) vs $(group_labels[groups[j]]):\t\t\t$(pvalue(mannwhitney))\n")
					end
				end
			end
		end
		println("#############################################################################")
		write(f, "#############################################################################\n")
		println("")
		close(f)

		# krustal-wallis test for nAUC
		f = open("./analysis/$(type)_$(r)_nAUC.txt", "w")
		println("$(type) $(r) nAUC: Analysis #####################################################")
		write(f, "$(type) $(r) nAUC: Analysis #####################################################\n")
		kruskal_nauc = KruskalWallisTest(grouped_nauc["A"], grouped_nauc["B"], grouped_nauc["C"], grouped_nauc["D"], grouped_nauc["E"], grouped_nauc["F"], grouped_nauc["G"], grouped_nauc["H"])
		println("Kruskal-Wallis test p-value: $(pvalue(kruskal_nauc))")
		write(f, "Kruskal-Wallis test p-value: $(pvalue(kruskal_nauc))\n\n")
		println("")
		if pvalue(kruskal_nauc) < 0.05
			#perform mann-whitney u test for each pair
			println("Significant differences (Mann-Whitney U test):")
			write(f, "Significant differences (Mann-Whitney U test):\n")
			for i in 1:8
				for j in i+1:8
					mannwhitney = MannWhitneyUTest(grouped_nauc[groups[i]], grouped_nauc[groups[j]])
					if pvalue(mannwhitney) < 0.05
						println("$(group_labels[groups[i]]) vs $(group_labels[groups[j]]):\t\t\t$(pvalue(mannwhitney))")
						write(f, "$(group_labels[groups[i]]) vs $(group_labels[groups[j]]):\t\t\t$(pvalue(mannwhitney))\n")
					end
				end
			end
		end
		println("##############################################################################")
		write(f, "##############################################################################\n")
		println("")
		close(f)
	end
end

println("Analysis performed.")