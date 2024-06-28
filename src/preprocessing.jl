# -*- coding: utf-8 -*-
# @author: max (max@mdotmonar.ch)

using HDF5
using DSP
using Statistics
using Plots
using CurveFit

# parameters
sampling_rate = 20000
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

function signal_reconstruct(signal, info)
	return [(s-info[:ADZero])*(info[:ConversionFactor]*(10.0^info[:Exponent])) for s in signal]
end

function high_pass_filter(signal, nco, order)
	filter = digitalfilter(Highpass(nco), Butterworth(order))
	return filtfilt(filter, signal)
end

function band_pass_filter(signal, nco_low, nco_high, order)
	filter = digitalfilter(Bandpass(nco_low, nco_high), Butterworth(order))
	return filtfilt(filter, signal)
end

function filter_and_resample(dataset, filter, resampling_rate)
	h5open("./preprocessed_data/"*dataset*"_checkerboard_preprocessed.h5", "cw") do preprocessed_file
		# check if average_signal group exists
		if group_check(preprocessed_file, ["electrode_0"])
			println("Skipping preprocessing...")
			return
		end

		# open raw file
		file = h5open("./raw_data/"*dataset*"_checkerboard.h5", "r")
		stream = read(file, "Data/Recording_0/AnalogStream/Stream_0")
		close(file)

		# open event list file
		csv_file = readlines("./raw_data/event_list_"*dataset*"_checkerboard_.csv")
		event_start = [split(line, ",") for line in csv_file][2][2]
		event_start = parse(Int, event_start)
		event_end = [split(line, ",") for line in csv_file][2][3]
		event_end = parse(Int, event_end)

		signal_length = trunc(Int, length(stream["ChannelData"])/252)	# 252 channels

		for i in 1:252
			println("Processing electrode_"*string(i-1)*"...")

			info = stream["InfoChannel"][i]
			signal_raw = stream["ChannelData"][signal_length*(i-1)+1:signal_length*i]

			signal = signal_reconstruct(signal_raw, info)
			signal = [s*1000 for s in signal] # Convert from V to mV

			signal = signal[event_start:event_end]

			# filter signal
			if filter["type"] == "HPF"
				filtered_signal = high_pass_filter(signal, filter["cutoff"]/nyquist_frequency, 5)
			elseif filter["type"] == "BPF"
				filtered_signal = band_pass_filter(signal, filter["cutoff_low"]/nyquist_frequency, filter["cutoff_high"]/nyquist_frequency, 5)
			end

			# resample signal
			step = trunc(Int, sampling_rate / resampling_rate)
			resampled_signal = filtered_signal[1:step:end]

			preprocessed_file["electrode_"*string(i-1)*"/data"] = resampled_signal
		end
		println("Done.")
	end
end