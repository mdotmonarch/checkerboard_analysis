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

datasets = []
for group in groups
	global datasets = [datasets; grouped_datasets[group]]
end

v_color = [:skyblue :skyblue :tomato :tomato :seagreen :seagreen :goldenrod :goldenrod]
v_color_index = Dict()
v_color_index["A"] = :skyblue 
v_color_index["B"] = :skyblue
v_color_index["C"] = :tomato
v_color_index["D"] = :tomato
v_color_index["E"] = :seagreen
v_color_index["F"] = :seagreen
v_color_index["G"] = :goldenrod
v_color_index["H"] = :goldenrod
v_fill = [0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5]
v_ls = [:dash :dashdot :dash :dashdot :dash :dashdot :dash :dashdot]
v_ls_index = Dict()
v_ls_index["A"] = :dash
v_ls_index["B"] = :dashdot
v_ls_index["C"] = :dash
v_ls_index["D"] = :dashdot
v_ls_index["E"] = :dash
v_ls_index["F"] = :dashdot
v_ls_index["G"] = :dash
v_ls_index["H"] = :dashdot

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
			entropy_data[dataset][type][r]["curve"] = entropy["curve"]
			entropy_data[dataset][type][r]["nAUC"] = entropy["nAUC"]
			entropy_data[dataset][type][r]["LRS"] = entropy["LRS"]
		end
	end

	close(file)
end

# plots
if !isdir("./plots")
	mkdir("./plots")
end

for type in t_list
	for r in r_list
		# check if directory exists
		if !isdir("./plots/$(type)_$(r)")
			mkdir("./plots/$(type)_$(r)")
		end

		# entropy curves per group, with average curve
		for g in groups
			plot(xlims=(1, 45), ylims=(0, 2.25), size=(800, 600), legend=:topright)
			plot!(xlabel="Scale", ylabel="Sample Entropy")
			plot!(title=" $(type) $(r)")

			avg_entropy_curve = zeros(45)
			for d in grouped_datasets[g]
				avg_entropy_curve += entropy_data[d][type][r]["curve"]
				plot!(1:45, entropy_data[d][type][r]["curve"], color=:black, lw=1, alpha=0.2, label=:none)
			end
			avg_entropy_curve /= length(grouped_datasets[g])
			plot!(1:45, avg_entropy_curve, color=v_color_index[g], lw=2, ls=v_ls_index[g], label=group_labels[g])
			savefig("./plots/$(type)_$(r)/entropy_curves_$(replace(group_labels[g], " " => "_")).png")
		end

		# average entropy curves of all groups
		plot(xlims=(1, 45), ylims=(0, 2.25), size=(800, 600), legend=:topright)
		plot!(xlabel="Scale", ylabel="Sample Entropy")
		plot!(title=" $(type) $(r)")
		for g in groups
			avg_entropy_curve = zeros(45)
			for d in grouped_datasets[g]
				avg_entropy_curve += entropy_data[d][type][r]["curve"]
			end
			avg_entropy_curve /= length(grouped_datasets[g])
			plot!(1:45, avg_entropy_curve, color=v_color_index[g], lw=2, ls=v_ls_index[g], label=group_labels[g])
		end
		savefig("./plots/$(type)_$(r)/average_entropy_curves.png")

		# LRS distribution
		local grouped_lrs = Dict()
		for g in groups
			grouped_lrs[g] = Float64[]
		end
		for g in groups
			for d in grouped_datasets[g]
				push!(grouped_lrs[g], entropy_data[d][type][r]["LRS"])
			end
		end
		plot(ylims=(-0.025, 0.025), size=(800, 600), legend=:none)
		a_data = [grouped_lrs[g] for g in groups]

		violin!(group_labels_plot, a_data, label=group_labels_plot, color = v_color, fill = v_fill, ls=v_ls)
		dotplot!(group_labels_plot, a_data, label=false, line = 0, marker=:black, side=:left, mode=:none)
		plot!(xlabel="Group", ylabel="LRS")

		savefig("./plots/$(type)_$(r)/lrs_distribution.png")

		# nAUC distribution
		local grouped_nauc = Dict()
		for g in groups
			grouped_nauc[g] = Float64[]
		end
		for g in groups
			for d in grouped_datasets[g]
				push!(grouped_nauc[g], entropy_data[d][type][r]["nAUC"])
			end
		end
		plot(ylims=(0, 2.0), size=(800, 600), legend=:none)
		a_data = [grouped_nauc[g] for g in groups]

		violin!(group_labels_plot, a_data, label=group_labels_plot, color = v_color, fill = v_fill, ls=v_ls)
		dotplot!(group_labels_plot, a_data, label=false, line = 0, marker=:black, side=:left, mode=:none)
		plot!(xlabel="Group", ylabel="nAUC")

		savefig("./plots/$(type)_$(r)/nauc_distribution.png")
	end
end

println("Plots saved.")