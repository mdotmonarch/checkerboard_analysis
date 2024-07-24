# -*- coding: utf-8 -*-
# @author: max (max@mdotmonar.ch)

using HDF5

include("preprocessing.jl")
include("processing.jl")

# parameters
filter = Dict(
	"type" => "BPF",
	"cutoff_low" => 0.1, # Hz
	"cutoff_high" => 40 # Hz
)
resampling_rate = 250 # Hz

r_list = [0.2]
e_f_list = ["none", "snr_0", "snr_1", "snr_2", "snr_3", "snr_4", "snr_5", "snr_6", "snr_7"]

# divide the datasets into N chunks
N = 8
i = parse(Int, ARGS[1])

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

for electrode_filter in e_f_list

	println("Processing datasets with electrode filter: ", electrode_filter)

	f_datasets = datasets

	if electrode_filter != "none"
		f_datasets = []
		for dataset in datasets
			#check if snr file exists
			if !isfile("../SNR/"*dataset*"_SNR.h5")
				println("SNR file not found for dataset: ", dataset)
				continue
			end
			snr_file = h5open("../SNR/"*dataset*"_SNR.h5", "r")
			threshold = parse(Float64, electrode_filter[5:end])
			at_least_one = false
			for i in 1:252
				if read(snr_file, "electrode_"*string(i-1)*"/SNR") >= threshold
					push!(f_datasets, dataset)
					break
				end
			end
		end
	else 
		f_datasets = datasets
	end

	println("Dataset length: ", length(f_datasets))

	# split according to N cores
	f_datasets = f_datasets[i:N:end]

	#=
	## PREPROCESSING
	for dataset in f_datasets
		println("Preprocessing dataset: ", dataset)
		filter_and_resample(dataset, filter, resampling_rate)
	end
	=#

	## PROCESSING
	for dataset in f_datasets
		println("Processing dataset: ", dataset)
		get_segments(dataset, 35)
		normalize_signals(dataset)
		get_electrode_mean(dataset, 35, electrode_filter)
	end

	## ENTROPY
	for dataset in f_datasets
		println("Computing entropy and complexity for dataset: ", dataset)
		for r in r_list
			compute_entropy_curve(dataset, electrode_filter, "RCMSE", 2, r, [i for i in 1:45])
		end
	end

end

println("Pipeline complete.")