# -*- coding: utf-8 -*-
# @author: max (max@mdotmonar.ch)

include("entropy.jl")

using HDF5
using DSP
using JSON
using Statistics
using Plots

# parameters
sampling_rate = 20000
resampling_rate = 250
nyquist_frequency = 0.5 * sampling_rate

function group_check(file, list, i=0)
	if i == length(list)
		return true
	end
	path = "/"*join(list[1:i], "/")
	if list[i+1] in keys(file[path])
		return group_check(file, list, i+1)
	end
	return false
end

function get_segments(dataset, time)
	# open processed file
	h5open("./processed_data/"*dataset*"_checkerboard_processed.h5", "cw") do processed_file
		# check
		if group_check(processed_file, ["electrode_0", "data"])
			println("Skipping getting signal segment..")
			return
		end

		println("Getting signal segment... ")

		h5open("./preprocessed_data/"*dataset*"_checkerboard_preprocessed.h5", "r") do preprocessed_file
			# get 35s segment
			for i in 1:252
				signal = read(preprocessed_file, "electrode_"*string(i-1)*"/data")
				processed_file["electrode_"*string(i-1)*"/data"] = signal[1:resampling_rate*time]
			end
		end

		println("Done.")
	end
end

function normalize_signals(dataset)
	# open processed file
	h5open("./processed_data/"*dataset*"_checkerboard_processed.h5", "cw") do processed_file
		# check
		if group_check(processed_file, ["electrode_0", "normalized"])
			println("Skipping signal segment normalization...")
			return
		end

		println("Normalize signal segment... ")

		for i in 1:252
			signal = read(processed_file, "electrode_"*string(i-1)*"/data")
			# normalize
			u = mean(signal)
			s = std(signal)
			normalized_signal = (signal .- u) ./ s
			processed_file["electrode_"*string(i-1)*"/normalized/data"] = normalized_signal
		end

		println("Done.")
	end
end

function get_mean_signal(dataset, time)
	# open processed file
	h5open("./processed_data/"*dataset*"_checkerboard_processed.h5", "cw") do processed_file
		# check
		if group_check(processed_file, ["mean_signal"])
			println("Skipping getting mean signal...")
			return
		end

		println("Getting mean signal... ")

		mean_signal = zeros(resampling_rate*time)
		for i in 1:252
			mean_signal += read(processed_file, "electrode_"*string(i-1)*"/normalized/data")
		end
		mean_signal = mean_signal ./ 252
		processed_file["mean_signal/data"] = mean_signal

		println("Done.")
	end
end

function compute_entropy_curve(dataset, type, m, r, scales)
	# open entropy file
	h5open("./entropy_data/"*dataset*"_checkerboard_entropy.h5", "cw") do entropy_file
		# check
		if group_check(entropy_file, [type, string(r)])
			println("Skipping computing "*type*" curve with r = "*string(r)*"...")
			return
		end

		println("Computing "*type*" curve with r = "*string(r)*"...")

		h5open("./processed_data/"*dataset*"_checkerboard_processed.h5", "r") do processed_file

			#=
			for i in 1:252
				println("Processing electrode_"*string(i-1)*"...")
				signal = read(processed_file, "electrode_"*string(i-1)*"/data")
				# compute entropy curve
				if type == "MSE"
					entropy_curve = multiscale_entropy(signal, m, r*std(signal), "sample", scales)
				elseif type == "RCMSE"
					entropy_curve = refined_composite_multiscale_entropy(signal, m, r*std(signal), "sample", scales)
				elseif type == "FMSE"
					entropy_curve = multiscale_entropy(signal, m, r*std(signal), "fuzzy", scales)
				elseif type == "FRCMSE"
					entropy_curve = refined_composite_multiscale_entropy(signal, m, r*std(signal), "fuzzy", scales)
				end
				entropy_file[type*"/"*string(r)*"/electrode_"*string(i-1)*"/curve"] = entropy_curve
				entropy_file[type*"/"*string(r)*"/electrode_"*string(i-1)*"/nAUC"] = compute_nAUC(entropy_curve)
				entropy_file[type*"/"*string(r)*"/electrode_"*string(i-1)*"/LRS"] = compute_LRS(entropy_curve, scales)
			end
			=#

			println("Processing mean signal...")
			signal = read(processed_file, "mean_signal/data")
			# compute entropy curve
			if type == "MSE"
				entropy_curve = multiscale_entropy(signal, m, r*std(signal), "sample", scales)
			elseif type == "RCMSE"
				entropy_curve = refined_composite_multiscale_entropy(signal, m, r*std(signal), "sample", scales)
			elseif type == "FMSE"
				entropy_curve = multiscale_entropy(signal, m, r*std(signal), "fuzzy", scales)
			elseif type == "FRCMSE"
				entropy_curve = refined_composite_multiscale_entropy(signal, m, r*std(signal), "fuzzy", scales)
			end
			entropy_file[type*"/"*string(r)*"/mean_signal/curve"] = entropy_curve
			entropy_file[type*"/"*string(r)*"/mean_signal/nAUC"] = compute_nAUC(entropy_curve)
			entropy_file[type*"/"*string(r)*"/mean_signal/LRS"] = compute_LRS(entropy_curve, scales)
		end

		println("Done.")
	end
end