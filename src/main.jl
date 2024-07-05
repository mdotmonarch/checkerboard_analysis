# -*- coding: utf-8 -*-
# @author: max (max@mdotmonar.ch)

using HDF5
using Plots
using StatsPlots

include("preprocessing.jl")
include("processing.jl")

# parameters
filter = Dict(
	"type" => "BPF",
	"cutoff_low" => 0.1, # Hz
	"cutoff_high" => 40 # Hz
)
resampling_rate = 250 # Hz

datasets = [
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
	# Double 3m male
	"MR-0575_nd4",
	"MR-0583_nd4",
	# Double 3m female
	"MR-0577_nd4",
	"MR-0579_nd4",
	"MR-0585_nd4",
	"MR-0586_nd4",
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

println("Dataset length: ", length(datasets))

for dataset in datasets
	println("Processing dataset: ", dataset)

	# preprocessing
	# filter_and_resample(dataset, filter, resampling_rate)

	# processing
	# get_segments(dataset, 35)
	# normalize_signals(dataset)
	# get_mean_signal(dataset, 35)

	# entropy and complexity
	compute_entropy_curve(dataset, "RCMSE", 2, 0.1, [i for i in 1:45])
	compute_entropy_curve(dataset, "RCMSE", 2, 0.2, [i for i in 1:45])
	compute_entropy_curve(dataset, "RCMSE", 2, 0.3, [i for i in 1:45])
	compute_entropy_curve(dataset, "RCMSE", 2, 0.4, [i for i in 1:45])
	compute_entropy_curve(dataset, "RCMSE", 2, 0.5, [i for i in 1:45])
end

# plots
include("plots.jl")

# analysis
include("analysis.jl")

println("Pipeline complete.")