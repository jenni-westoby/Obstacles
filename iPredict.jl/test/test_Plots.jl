using DataFrames
using Test
include("../src/Plots.jl")

dict = Dict("I1" => [1,2,3,4], "I2" => [3,3,3,3])
expected_out = DataFrame(I1 = [1,2,3,4], I2 = [3,3,3,3])

@test DictToDf(dict) == expected_out

input = DataFrame(Genes = ["G1", "G1"], Isoforms = ["I1", "I3"],
Cell1 = [5,0],Cell2 = [1,1])

@test MeanIsoPerGenePerCell(input) == 1.5

genes_dict = Dict("G1" => dict)
expected_out = [2]

@test MeanNumIsoformsPerGenePerCellReal(genes_dict)[1] == expected_out

genes_dict = Dict("G1" => dict, "G2" => Dict("I1" => [0,0,1,1]))

@test MeanNumIsoformsPerGenePerCellReal(genes_dict)[1] == [2, 0.5]
@test MeanNumIsoformsPerGenePerCellReal(genes_dict)[2] == Dict("G1" => 2, "G2" => 0.5)
