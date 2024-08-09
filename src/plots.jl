# -*- coding: utf-8 -*-
# @author: max (max@mdotmonar.ch)

using HDF5
using Plots
using StatsPlots
using DSP

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
r_list = ["0.2"]
e_f_list = [
	"none",
	"snr_0",
	"snr_1",
	"snr_2",
	"snr_3",
	"snr_4",
	"snr_5",
	"snr_6",
	"snr_7"
]

# Signal and spectrogram
function plot_signal(signal, fs=250)
	plot!([i/fs for i in 1:length(signal)], signal, xlabel="Time (s)", ylabel="Amplitude")
end

function plot_spectrogram(signal, fs=250)
	nw = length(signal) ÷ 30
	spec = spectrogram(signal, nw, nw÷2, fs=fs, window=hanning)
	heatmap!(spec.time, spec.freq, 10*log.(spec.power), c=:viridis, xlabel="Time (s)", ylabel="Frequency (Hz)")
end

for dataset in datasets
	if !isdir("./plots")
		mkdir("./plots")
	end
	if !isdir("./plots/$(dataset)")
		mkdir("./plots/$(dataset)")
	end
	if !isdir("./plots/$(dataset)/signals")
		mkdir("./plots/$(dataset)/signals")
	end
	if !isdir("./plots/$(dataset)/signals/electrodes")
		mkdir("./plots/$(dataset)/signals/electrodes")
	end

	### open processed file
	processed_file = h5open("./processed_data/$(dataset)_checkerboard_processed.h5", "r")

	### plot signal and spectrogram for each electrode
	#=
	for n in 1:252
		signal = read(processed_file, "electrode_$(n-1)/normalized/data")

		p1 = plot(title="Normalized signal and spectrogram of $(dataset), electrode $(n-1)", legend=:none)
		plot_signal(signal)

		p2 = plot()
		plot_spectrogram(signal)

		plot(p1, p2, layout=grid(2, 1, heights=[0.2 ,0.8]), size=(800, 900), dpi=300)
		savefig("./plots/$(dataset)/signals/electrodes/electrode_$(n-1).png")
	end
	=#

	### plot signal and spectrogram for electrode mean
	e_f_computed = keys(processed_file["electrode_mean"])

	for e_f in e_f_computed
		signal = read(processed_file, "electrode_mean/$(e_f)/data")

		p1 = plot(title="Normalized signal and spectrogram of $(dataset), \nelectrode mean w/ electrode filter: $(e_f)", legend=:none)
		plot_signal(signal)

		p2 = plot()
		plot_spectrogram(signal)

		plot(p1, p2, layout=grid(2, 1, heights=[0.2 ,0.8]), size=(800, 900), dpi=300)
		savefig("./plots/$(dataset)/signals/electrode_mean_$(e_f).png")
	end

	### close processed file
	close(processed_file)
end

# load entropy data
entropy_data = Dict()

for dataset in datasets
	entropy_data[dataset] = Dict()
	for t in t_list
		entropy_data[dataset][t] = Dict()
		for r in r_list
			entropy_data[dataset][t][r] = Dict()
		end
	end
end

for dataset in datasets
	file = h5open("./entropy_data/$(dataset)_checkerboard_entropy.h5", "r")

	for type in t_list
		for r in r_list
			entropy = read(file["/$(type)/$(r)/electrode_mean"])
			for e_f in keys(entropy)
				entropy_data[dataset][type][r][e_f] = Dict()
				entropy_data[dataset][type][r][e_f]["curve"] = entropy[e_f]["curve"]
				entropy_data[dataset][type][r][e_f]["LRS_all"] = compute_LRS(entropy[e_f]["curve"], [i for i in 1:45])
				entropy_data[dataset][type][r][e_f]["LRS_15"] = compute_LRS(entropy[e_f]["curve"][1:15], [i for i in 1:15])
				entropy_data[dataset][type][r][e_f]["LRS_30"] = compute_LRS(entropy[e_f]["curve"][16:30], [i for i in 16:30])
				entropy_data[dataset][type][r][e_f]["LRS_45"] = compute_LRS(entropy[e_f]["curve"][31:45], [i for i in 31:45])
				entropy_data[dataset][type][r][e_f]["nAUC_all"] = compute_nAUC(entropy[e_f]["curve"])
				entropy_data[dataset][type][r][e_f]["nAUC_15"] = compute_nAUC(entropy[e_f]["curve"][1:15])
				entropy_data[dataset][type][r][e_f]["nAUC_30"] = compute_nAUC(entropy[e_f]["curve"][16:30])
				entropy_data[dataset][type][r][e_f]["nAUC_45"] = compute_nAUC(entropy[e_f]["curve"][31:45])
			end
		end
	end

	close(file)
end

### plot SNR comparison
for dataset in datasets
	for t in t_list
		for r in r_list
			if !isdir("./plots/$(dataset)/entropy")
				mkdir("./plots/$(dataset)/entropy")
			end
			global l = @layout [a
				[grid(1,4)]
				[grid(1,4)]
			]
			p1 = plot(xlims=(1, 45))
			plot!(xlabel="Scale", ylabel="SampEn")
			plot!(title="RCMSE curves of $(dataset)\ncomparison of electrode filters")
			vspan!([1, 15], color=:black, alpha=0.1, label=:none)
			vspan!([16, 30], color=:black, alpha=0.1, label=:none)
			vspan!([31, 45], color=:black, alpha=0.1, label=:none)
			for e_f in ["none", "snr_3", "snr_7"]
				if !haskey(entropy_data[dataset][t][r], e_f)
					continue
				end
				entropy_curve = entropy_data[dataset][t][r][e_f]["curve"]
				scales = [i for i in 1:45]
				plot!(scales, entropy_curve, label=e_f)
			end
			# LRS
			p2 = []
			for (segment, lim) in zip(["all", "15", "30", "45"], ["1:45", "1:15", "16:30", "31:45"])
				p = plot(framestyle = :origin, xlims=(-1, 1), grid=false, xaxis=false, title="LRS "*lim)
				for e_f in ["none", "snr_3", "snr_7"]
					if !haskey(entropy_data[dataset][t][r], e_f)
						continue
					end
					lrs = entropy_data[dataset][t][r][e_f]["LRS_"*segment]
					scatter!([0], [lrs], label=e_f, ms=5, markerstrokewidth=0)
				end
				p2 = [p2; p]
			end
			# nAUC
			p3 = []
			for (segment, lim) in zip(["all", "15", "30", "45"], ["1:45", "1:15", "16:30", "31:45"])
				p = plot(framestyle = :origin, xlims=(-1, 1), grid=false, xaxis=false, title="nAUC "*lim)
				for e_f in ["none", "snr_3", "snr_7"]
					if !haskey(entropy_data[dataset][t][r], e_f)
						continue
					end
					nauc = entropy_data[dataset][t][r][e_f]["nAUC_"*segment]
					scatter!([0], [nauc], label=e_f, ms=5, markerstrokewidth=0)
				end
				p3 = [p3; p]
			end
			plot(p1, p2..., p3..., layout=l, size=(900, 900), dpi=300)
			savefig("./plots/$(dataset)/entropy/$(t)_$(r)_comparison.png")
		end
	end
end

for type in t_list
	for r in r_list
		if !isdir("./plots/$(type)_$(r)")
			mkdir("./plots/$(type)_$(r)")
		end
		for (segment, lim) in zip(["all", "15", "30", "45"], ["1-45", "1-15", "16-30", "31-45"])
			for e_f in ["none", "snr_3", "snr_7"]
				# LRS distribution
				local grouped_lrs = Dict()
				for group in groups
					grouped_lrs[group] = Float64[]
				end
				for group in groups
					for dataset in grouped_datasets[group]
						if !haskey(entropy_data[dataset][type][r], e_f)
							continue
						end
						push!(grouped_lrs[group], entropy_data[dataset][type][r][e_f]["LRS_"*segment])
					end
				end
				plot(size=(1000, 600), dpi=300, legend=:none)
				a_data = [grouped_lrs[group] for group in groups]

				violin_labels = [group_labels[group]*" ($(length(grouped_lrs[group])))" for group in groups]
				violin_labels = reshape(violin_labels, 1, length(violin_labels))

				violin!(violin_labels, a_data, label=violin_labels, color = v_color, fill = v_fill, ls=v_ls)
				dotplot!(violin_labels, a_data, label=false, line = 0, marker=:black, side=:left, mode=:none, alpha=0.3)
				plot!(xlabel="Group", ylabel="LRS", title="LRS distribution of $(type) $(r) $(lim), electrode filter: $(e_f)")

				println("subject count:", length([(a_data...)...]) )

				savefig("./plots/$(type)_$(r)/LRS_distribution_$(e_f)_$(lim).png")
			end
		end
	end
end

println("Plots saved.")