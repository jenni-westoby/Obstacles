include("../src/pValues.jl")
using Test

test_dict = Dict(1 => Float64[1,2,3,4])

#These tests need updating
@test FindpValues(test_dict, 2.5, 1) == Dict(1 => 0.5)
@test FindpValues(test_dict, 6, 1) == Dict(1 => 0)
@test FindpValues(test_dict, 0, 1) == Dict(1 => 1)

test_dict = Dict(1 => Float64[1,2,3,4], 2 => Float64[5,6,7,8])
@test FindpValues(test_dict, 2.5, 1) == Dict(1 => 0.5, 2 => 1)

#This test fails but output is correct
@test FindpValues(test_dict, 2.5, 40) == Dict(1 => missing, 2 => missing)

meanArr = [1,2,3,4]
mean_val = 1
threshold = 1

@test GlobalpValues(meanArr, mean_val, 1.5) == 0.75
