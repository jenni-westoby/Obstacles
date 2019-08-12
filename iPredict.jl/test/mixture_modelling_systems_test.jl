include("iPredict.jl/src/MixtureModelling.jl")
using DelimitedFiles
using Distributions
using StatsPlots
using DataFrames

test_data_1 = rand(LogNormal(1, 0.5), 20)
test_data_2 = rand(LogNormal(3, 0.5), 60)
append!(test_data_1, test_data_2)

relationships = Array{Float64}(undef, length(test_data_1), 2)
mu_arr = [1, 3]
sigma_arr = [0.5, 0.5]
k_arr = [0.5, 0.5]

expectation!(relationships, length(test_data_1), 2, mu_arr, sigma_arr, k_arr, test_data_1)
maximisation_k!(relationships, k_arr)

append!(test_data_1, rand(LogNormal(2, 0.7), 80))

relationships = Array{Float64}(undef, length(test_data_1), 3)
mu_arr = [1, 3, 2]
sigma_arr = [0.5, 0.5, 0.7]
k_arr = [1/3, 1/3, 1/3]

expectation!(relationships, length(test_data_1), 3, mu_arr, sigma_arr, k_arr, test_data_1)
maximisation_k!(relationships, k_arr)
