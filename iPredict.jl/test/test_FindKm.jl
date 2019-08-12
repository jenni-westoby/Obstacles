using Test
using CSV
include("../src/FindKm.jl")

################################################################################
# Tests for functions that normalise data like M3Drop expects
################################################################################

#Check can accurately find no. detected genes per cell
df = DataFrame(Column1 = ["I1", "I2", "I3", "I4"], Cell1 = Float64[0, 0, 0, 0],
Cell2 = Float64[0,1,1,1])
expected_output_true = DataFrame(Cells = [Symbol("Cell1"), Symbol("Cell2")],
DetectedGenes = Float64[0,3])
@test FindNumDetectedGenesPerCell(df, true) == expected_output_true

#Check can accurately find no. undetected genes per cell
expected_output_false = DataFrame(Cells = [Symbol("Cell1"), Symbol("Cell2")],
DetectedGenes = Float64[4,1])
@test FindNumDetectedGenesPerCell(df, false) == expected_output_false

#Check can find mean DetectedGenes
@test FindMean(expected_output_true) == 1.5
@test FindMean(expected_output_false) == 2.5

#Check can find sd DetectedGenes
@test FindStd(expected_output_true) ≈ 2.12 atol=0.01
@test FindStd(expected_output_false) ≈ 2.12 atol=0.01

#Check removing low quality cells works
expected_output = deepcopy(expected_output_false)
@test RemoveLowQuality(expected_output_false, 2.5, 2.12) == expected_output
df = DataFrame(Cells = [Symbol("Cell1")], DetectedGenes = Float64[8025])
@test RemoveLowQuality(df, 0, 1 ) == DataFrame(Cells = Symbol[],
DetectedGenes = Float64[])

#Check can remove cells correctly
orig_df = DataFrame(Isoforms = ["I1", "I2", "I3", "I4"], Cell1 = Float64[0, 0, 0, 0],
Cell2 = Float64[0,1,1,1])
@test KeepFiltered(df, orig_df) == DataFrame(Isoforms = ["I1", "I2", "I3", "I4"],
 Cell1 = Float64[0, 0, 0, 0])
@test KeepFiltered(DataFrame(Cells = []), orig_df) ==
DataFrame(Isoforms = ["I1", "I2", "I3", "I4"])

#Check can correctly remove undetected genes
test_data = DataFrame(Column1 = ["I1", "I2", "I3"], Cell1 = Float64[0, 0, 1],
Cell2 = Float64[0,1,1], Cell3 = Float64[1,1,1], Cell4 = [1,1,1])
expected_output = DataFrame(Column1 = ["I3"], Cell1 = Float64[1],
Cell2 = Float64[1], Cell3 = Float64[1], Cell4 = Float64[1])
@test RemoveUndetectedGenes(test_data) == expected_output

#Check cpm calculation is correct
input = DataFrame(Isoforms = ["I1", "I2", "I3", "I4"], Cell1 = Float64[0, 0, 0, 1],
Cell2 = Float64[0,0,1,1])
expected_output = DataFrame(Isoforms = ["I1", "I2", "I3", "I4"], Cell1 = Float64[0, 0, 0, 1000000],
Cell2 = Float64[0,0,500000,500000])
@test ConvertToCpm(input) == expected_output

#Check can correctly remove lowly expressed genes
expected_output = DataFrame(Isoforms = ["I3", "I4"], Cell1 = Float64[ 0, 1],
Cell2 = Float64[1,1])
@test RemoveLowlyExpressedGenes(input) == expected_output

################################################################################
# Tests for functions that calculate Km
################################################################################

input = DataFrame(Isoforms = ["I1", "I2"], Cell1 = [0, 1], Cell2 = [1, 1])
expected_output = DataFrame(Isoforms = ["I1", "I2"], p = [0.5, 0])
@test Findp(input) == expected_output

s_out = deepcopy(expected_output)
s_out[:s] = [0.5, 1]
@test Finds(expected_output, input) == s_out

#Check LogLikelihood function gives correct output
LogLikelihoodFn = LogLikelihood(Float64[1,2,3], Float64[1,2,3])
@test LogLikelihoodFn([1.0,1.0]) ≈ 8.05 atol=0.01

#Check maximum likelihood estimation is working
df = CSV.read("test/data/p.txt", header = true, delim = ',')
p=df.x
df = CSV.read("test/data/s.txt", header = true, delim = ',')
s=df.x
parameters = MaxLikelihoodEst(p, s)
@test parameters[1] ≈ 8.5 atol=0.1
@test parameters[2] ≈ 0.1 atol=0.1
