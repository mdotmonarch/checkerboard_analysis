# -*- coding: utf-8 -*-
# @author: max (max@mdotmonar.ch)

include("entropy.jl")

using HDF5
using DSP
using Statistics

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
	end
	
	println("Done.")
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
	end

	println("Done.")
end

function get_electrode_mean(dataset, time, electrode_filter="none") #filters: none, snr_n
	# get electrode mean
	h5open("./processed_data/"*dataset*"_checkerboard_processed.h5", "cw") do processed_file
		# check
		if group_check(processed_file, ["electrode_mean", electrode_filter])
			println("Skipping getting electrode mean...")
			return
		end

		println("Getting electrode mean... ")

		mean_signal = zeros(resampling_rate*time)
		if electrode_filter == "none"
			for i in 1:252
				mean_signal += read(processed_file, "electrode_"*string(i-1)*"/normalized/data")
			end
			mean_signal = mean_signal ./ 252
			processed_file["electrode_mean/none/data"] = mean_signal
		else
			threshold = parse(Float64, electrode_filter[5:end])
			# add SNR filter based on SNR h5 file
			snr_file = h5open("../SNR/"*dataset*"_SNR.h5", "r")

			electrode_number = 0

			for i in 1:252
				if read(snr_file, "electrode_"*string(i-1)*"/SNR") >= threshold
					electrode_number += 1
					mean_signal += read(processed_file, "electrode_"*string(i-1)*"/normalized/data")
				end
			end

			mean_signal = mean_signal ./ electrode_number
			processed_file["electrode_mean/"*electrode_filter*"/data"] = mean_signal
		end
	end

	println("Done.")
end

function compute_entropy_curve(dataset, e_f, type, m, r, scales)

	h5open("./entropy_data/"*dataset*"_checkerboard_entropy.h5", "cw") do entropy_file
		# check
		if group_check(entropy_file, [type, string(r), "electrode_0"])
			println("Skipping computing "*type*" curve with r = "*string(r)*" for all electrodes...")
			return
		end

		println("Computing "*type*" curve with r = "*string(r)*" for all electrodes...")

		h5open("./processed_data/"*dataset*"_checkerboard_processed.h5", "r") do processed_file

			# compute entropy curve for all events, for all electrodes
			for i in 0:251
				signal = read(processed_file, "electrode_"*string(i)*"/normalized/data")
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
				entropy_file[type*"/"*string(r)*"/electrode_"*string(i)*"/curve"] = entropy_curve
				entropy_file[type*"/"*string(r)*"/electrode_"*string(i)*"/nAUC"] = compute_nAUC(entropy_curve)
				entropy_file[type*"/"*string(r)*"/electrode_"*string(i)*"/LRS"] = compute_LRS(entropy_curve, scales)
			end
		end
	end

	h5open("./entropy_data/"*dataset*"_checkerboard_entropy.h5", "cw") do entropy_file
		# check
		if group_check(entropy_file, [type, string(r), "electrode_mean", e_f])
			println("Skipping computing "*type*" curve with r = "*string(r)*" for electrode mean with electrode filter = "*e_f*"...")
			return
		end

		println("Computing "*type*" curve with r = "*string(r)*" for electrode mean with electrode filter = "*e_f*"...")

		h5open("./processed_data/"*dataset*"_checkerboard_processed.h5", "r") do processed_file

			# compute entropy curve for all events, for all electrodes
			signal = read(processed_file, "electrode_mean/"*e_f*"/data")
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
			entropy_file[type*"/"*string(r)*"/electrode_mean/"*e_f*"/curve"] = entropy_curve
			entropy_file[type*"/"*string(r)*"/electrode_mean/"*e_f*"/nAUC"] = compute_nAUC(entropy_curve)
			entropy_file[type*"/"*string(r)*"/electrode_mean/"*e_f*"/LRS"] = compute_LRS(entropy_curve, scales)
		end
	end

	println("Done.")
end