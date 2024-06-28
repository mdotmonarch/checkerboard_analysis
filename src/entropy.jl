# -*- coding: utf-8 -*-
# @author: max (max@mdotmonar.ch)

using Statistics
using Combinatorics
using Trapz
using CurveFit

function chebyshev_distance(x, y)
	return maximum(abs.(x - y))
end

function generate_windows(signal, m)
	N = length(signal)
	return [signal[i:i + m - 1] for i in 1:N - m + 1]
end

function sample_entropy_matches(signal, m, r)
	# generate windows from signal
	m_vector = generate_windows(signal, m)
	m1_vector = generate_windows(signal, m + 1)

	# compute the number of matches
	A = sum([(chebyshev_distance(i, j) <= r) ? 1 : 0 for (i, j) in combinations(m1_vector, 2)])
	B = sum([(chebyshev_distance(i, j) <= r) ? 1 : 0 for (i, j) in combinations(m_vector, 2)])

	return A, B
end

function sample_entropy(signal, m, r)
	A, B = sample_entropy_matches(signal, m, r)
	return -log(A/B)
end

function fuzzy_entropy_matches(signal, m, r)
	# generate windows from signal
	m_vector = generate_windows(signal, m)
	m1_vector = generate_windows(signal, m + 1)

	# compute the number of matches
	A = sum([exp( -log(2)*( (chebyshev_distance(i, j)/r)^2 ) ) for (i, j) in combinations(m1_vector, 2)])
	B = sum([exp( -log(2)*( (chebyshev_distance(i, j)/r)^2 ) ) for (i, j) in combinations(m_vector, 2)])

	return A, B
end

function fuzzy_entropy(signal, m, r)
	A, B = fuzzy_entropy_matches(signal, m, r)
	return -log(A/B)
end

function multiscale_entropy(signal, m, r, e, scales = [i for i in 1:trunc(Int, length(signal)/(m+10))])
	N = length(signal)

	en_list = Float64[]

	for scale in scales
		# coarse-graining
		ratio = trunc(Int, N / scale)
		coarse_signal = zeros(ratio)
		for i in 1:ratio
			coarse_signal[i] = sum(signal[(i - 1) * scale + 1:i * scale]) / scale
		end

		# entropy
		if e == "sample"
			push!(en_list, sample_entropy(coarse_signal, m, r))
		elseif e == "fuzzy"
			push!(en_list, fuzzy_entropy(coarse_signal, m, r))
		end
	end

	return en_list
end

function composite_multiscale_entropy(signal, m, r, e, scales = [i for i in 1:trunc(Int, length(signal)/(m+10))])
	N = length(signal)

	en_list = Float64[]

	for scale in scales
		# multiple coarse-graining
		cumulative_en = 0
		tau = length(signal)%scale
		for t in 0:tau
			ratio = trunc(Int, N / scale)
			coarse_signal = zeros(ratio)
			for i in 1:ratio
				coarse_signal[i] = sum(signal[((i - 1) * scale)+1+t:(i * scale)+t]) / scale
			end

			# entropy
			if e == "sample"
				cumulative_en += sample_entropy(coarse_signal, m, r)
			elseif e == "fuzzy"
				cumulative_en += fuzzy_entropy(coarse_signal, m, r)
			end
		end

		# sample entropy
		push!(en_list, cumulative_en/(tau+1))
	end

	return en_list
end

function refined_composite_multiscale_entropy(signal, m, r, e, scales = [i for i in 1:trunc(Int, length(signal)/(m+10))])
	N = length(signal)

	en_list = Float64[]

	for scale in scales
		# cumulative matches
		A = 0
		B = 0

		# multiple coarse-graining
		tau = length(signal)%scale
		for t in 0:tau
			ratio = trunc(Int, N / scale)
			coarse_signal = zeros(ratio)
			for i in 1:ratio
				coarse_signal[i] = sum(signal[((i - 1) * scale)+1+t:(i * scale)+t]) / scale
			end

			# entropy
			if e == "sample"
				A_m, B_m = sample_entropy_matches(coarse_signal, m, r)
			elseif e == "fuzzy"
				A_m, B_m = fuzzy_entropy_matches(coarse_signal, m, r)
			end
			A += A_m
			B += B_m
		end

		en = -log(A/B)

		# sample entropy
		push!(en_list, en)
	end

	return en_list
end

function compute_nAUC(curve)
	#filter out NaN values
	curve = curve[.!isnan.(curve)]
	#filter out Inf values
	curve = curve[.!isinf.(curve)]

	return trapz([i for i in 1:length(curve)], curve)/length(curve)
end

function compute_LRS(curve, scales)
	a, b = linear_fit(scales, curve)
	return b
end