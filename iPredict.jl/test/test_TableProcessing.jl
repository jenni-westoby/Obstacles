using Test
using CSV
include("../src/TableProcessing.jl")
#include("src/TableProcessing.jl")

#Check that reading in test data works
df = CSV.read("test/data/test_counts_matrix.txt",
header = true, delim = '\t')
@test ReadCountsMatrix("test/data/test_counts_matrix.txt") == df

#Check that an error is thrown when no filename is passed
@test_throws MethodError ReadCountsMatrix()

#Same tests for gene-isoform table
genes = CSV.read("test/data/all_isoforms.txt", header = false, delim = ' ')
@test ReadGeneIsoformRelationships("test/data/all_isoforms.txt") == genes
@test_throws MethodError ReadGeneIsoformRelationships()

#Check function returns true when genes match isoforms
df = DataFrame(Column1 = ["I1", "I2", "I3", "I4"], Cell1 = [0, 0, 0, 0])
genes = DataFrame(Column1 = ["G1", "G1", "G2", "G3"],
Column2 = ["I2", "I3", "I1", "I4"])
@test CheckGenesMatchIsoforms(df, genes) == true

#Check function returns false when genes don't match isoforms
@test CheckGenesMatchIsoforms(df, genes[1:2, :]) == false

#Check we get expected output when adding gene information
df_test = DataFrame(Genes = ["G1", "G1", "G2", "G3"],
Isoforms = ["I2", "I3", "I1", "I4"], Cell1 = [0, 0, 0, 0])
@test AddGeneInfo(df, genes) == df_test

#Check that we correctly convert from dataframe to dict
df_test = DataFrame(Genes = ["G1", "G1", "G2", "G3"],
Isoforms = ["I2", "I3", "I1", "I4"], Cell1 = [5, 5, 0, 5], Cell2 = [5,5,5,5])

dict_test_two = Dict("G1" => Dict("I2" => [5,5], "I3" => [5,5]))
dict_test_one = Dict("G3" => Dict("I4" => [5,5]))

#Check that we keep expected genes when N = 1 to N = 4
@test KeepGenesWithNIsoforms(df_test, 1) == dict_test_one
@test KeepGenesWithNIsoforms(df_test, 2) == dict_test_two

#Check that an error is thrown when N<1
@test_throws ArgumentError KeepGenesWithNIsoforms(df_test, 0)
@test_throws ArgumentError KeepGenesWithNIsoforms(df_test, -10)

#Check that it works when two of three isoforms are expressed
df_test = DataFrame(Genes = ["G1", "G1", "G1", "G3"],
Isoforms = ["I2", "I3", "I1", "I4"], Cell1 = [5, 5, 0, 5], Cell2 = [5,5,5,5])
@test KeepGenesWithNIsoforms(df_test, 2) == dict_test_two

df_test = DataFrame(Genes = ["G1", "G1", "G1"],
Isoforms = ["I2", "I3", "I1"], Cell1 = [5, 5, 0], Cell2 = [5,5,5])
dict_test_two = Dict("G1" => Dict("I2" => [5,5], "I3" => [5,5]))
@test ConvertExpressedIsosToDict(df_test) == dict_test_two
