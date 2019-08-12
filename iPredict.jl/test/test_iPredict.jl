include("../src/iPredict.jl")
#include("src/iPredict.jl")
using .iPredict
using Test

#An extremely general systems test - make sure the output is a dictionary!
result = globalPredict("test/data/small_counts.txt",
"test/data/small_gene_isoform_relationships.txt", 1, true)

@test isa(result, Dict) == true

#Check an error is thrown when the gene names don't match
@test_throws ArgumentError globalPredict("test/data/test_counts_matrix.txt",
"test/data/all_isoforms.txt", 1, false)

#Do some general systems checks for oneGenePrediction too
result = oneGenePrediction("test/data/small_counts.txt",
"test/data/small_gene_isoform_relationships.txt", "G1", 10, true)

@test isa(result, Dict) == true
@test length(result) == 1

result = oneGenePrediction("test/data/alt_small_counts.txt",
"test/data/small_gene_isoform_relationships.txt", "G1", 10, true)

@test isa(result, Dict) == true
@test length(result) == 2
