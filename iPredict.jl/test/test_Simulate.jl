using Test
include("../src/Simulate.jl")

#Check we can correctly find P(dropout)
test_dict = Dict("I1" => [1,1,1], "I2" => [5,0,0], "I3" => [1,1,0])
@test FindPDropout(test_dict, 1.0) == [0.5,(1-15/24), (1-6/15)]

#Check we can correctly rank isoforms
@test FindRank(test_dict) == [2,3,1]

#Check we can correctly find p(Obs)
@test FindPObs(test_dict) == [1, 1/3, 2/3]

#Check we can correctly find the harmonic mean
@test FindHarmonicMean(2) == exp(-2.25) + (0.5 * exp(-4))

#Check we can correctly find the median frequencies
@test FindMedianFrequencies(2, 1) == [exp(-(1.5 * 1.5)), (exp(-4))/2]

#Check we can correctly scale frequencies
@test ScaleFrequencies([0.15, 0.1]) == [0.6, 0.4]
@test ScaleFrequencies([1.2, 0.8]) == [0.6, 0.4]

#Check we can correctly find cumulative frequencies
@test FindCumulativeFrequencies([0.5, 0.1, 0.2, 0.2]) ≈ [0.5, 0.7,0.9,1.0]

#Check can return rank of picked isoform correctly
@test PickRank([0.5,0.75,1.0], 0.8) == 3
@test PickRank([0.5,0.75,1.0], 0.1) == 1

@test RescaleFreqs([0.1,0.4,1.0], 3) == [0.8, 1.0]

@test RescaleRankings([4,1,3,2], 1, 2) == [3,2,1]
@test RescaleRankings([4,1,3,2], 4, 1) == [1,3,2]

Random.seed!(1234)
out = RandomIsoformChoiceProbabilities(2,2)

@test out[1,1] ≈ 0.5106336456707476 atol=0.001
@test out[1,2] ≈ 0.4893663543292524 atol=0.001
@test out[2,1] ≈ 0.624996535934181 atol=0.001
@test out[2,2] ≈ 0.3750034640658189 atol=0.001

#Tests for QuantErrorsForExpressed
pFP = 0.05
pFN = 0.1
@test QuantErrorsForExpressed(true, 0.05) == false
@test QuantErrorsForExpressed(true, 0.1) == true
@test QuantErrorsForExpressed(false, 0.01) == true
@test QuantErrorsForExpressed(false, 0.1) == false

Random.seed!(1234)
@test RandomQuantErrorsForExpressed(true, 0.5) == true
@test RandomQuantErrorsForExpressed(true, 0.2) == false
@test RandomQuantErrorsForExpressed(false, 0.2) == true
@test RandomQuantErrorsForExpressed(false, 0.3) == false

Random.seed!(1234)
@test RandomQuantErrorsForUnexpressed(0.3) == false
@test RandomQuantErrorsForUnexpressed(0.3) == true

Random.seed!(1234)
unexprRandArr = rand(1,1,1,2)
NumSimulations = 1
genesList = ["G1"]
numCells = 1
unexpr_dict = Dict("G1" => 2)
filteringThreshold = 4
NumIsoformsToSimulate =4
out = rand(1,1,1)
out[1,1,1] = 0

@test UnexpressedErrorsIsoformDependence(unexprRandArr, NumSimulations,
genesList, numCells, unexpr_dict, NumIsoformsToSimulate, filteringThreshold) == out

unexprRandArr[1,1,1,1] = 0.05
out[1,1,1] = 1

@test UnexpressedErrorsIsoformDependence(unexprRandArr, NumSimulations,
genesList, numCells, unexpr_dict, NumIsoformsToSimulate, filteringThreshold) == out

readCapturedIsoforms = Array{Bool}(undef, 1,1,1,4)
readCapturedIsoforms[1,1,1,1] = true
readCapturedIsoforms[1,1,1,2] = true
readCapturedIsoforms[1,1,1,3] = true
readCapturedIsoforms[1,1,1,4] = false

Random.seed!(1234)
randIsoChoiceArr = rand(1,1,1,4)

@test ExpressedErrorsIsoformDependence(readCapturedIsoforms, randIsoChoiceArr,
NumSimulations, genesList, numCells, unexpr_dict, NumIsoformsToSimulate,
filteringThreshold) == readCapturedIsoforms

randIsoChoiceArr[1,1,1,1] = 0.03
output = Array{Bool}(undef, 1,1,1,4)
output[1,1,1,1] = false
output[1,1,1,2] = true
output[1,1,1,3] = true
output[1,1,1,4] = false

@test ExpressedErrorsIsoformDependence(readCapturedIsoforms, randIsoChoiceArr,
NumSimulations, genesList, numCells, unexpr_dict, NumIsoformsToSimulate,
filteringThreshold) == output

randIsoChoiceArr[1,1,1,4] = 0.05
output[1,1,1,4] = true

@test ExpressedErrorsIsoformDependence(readCapturedIsoforms, randIsoChoiceArr,
NumSimulations, genesList, numCells, unexpr_dict, NumIsoformsToSimulate,
filteringThreshold) == output

###############################################################################
#Test that we can correctly pick isoforms

Random.seed!(1234)
randIsoChoiceArr = rand(1,1,1,1) #0.59
IsoChoiceProbabilities = [0.1, 0.2]
pObsDict = Dict("G1" => IsoChoiceProbabilities)
NumSimulations = 1
genesList = ["G1"]
rankingDict = Dict("G1" => [2,1])
pDropoutDict = Dict("G1" => [0.1, 0.05])
numCells = 1
NumIsoformsToSimulate = 1

Random.seed!(1234) #Seed so sampling is reproducible

testArr = pickIsoforms(IsoChoiceProbabilities, NumSimulations,
    genesList, rankingDict, pDropoutDict, numCells,
    NumIsoformsToSimulate)

@test testArr[1,1,1,1] == 0.1

IsoChoiceProbabilities = rand(1, 2)
IsoChoiceProbabilities[1,1] = 0.1
IsoChoiceProbabilities[1,2] = 0.2

Random.seed!(1234) #Seed so sampling is reproducible

testArr = pickRandIsoforms(IsoChoiceProbabilities, NumSimulations,
    genesList, pDropoutDict, numCells,
    NumIsoformsToSimulate)

@test testArr[1,1,1,1] == 0.05

Random.seed!(1234)
testArr[1,1,1,1] = 0.0

testArr = PickUnifObsIsoforms(NumIsoformsToSimulate, genesList, pObsDict,
pDropoutDict, numCells, NumSimulations)

@test testArr[1,1,1,1] == 0.05

Random.seed!(1234)
testArr[1,1,1,1] = 0.0

testArr = PickCellVarIsoforms(NumIsoformsToSimulate, genesList, pObsDict,
pDropoutDict, numCells, NumSimulations)

@test testArr[1,1,1,1] == 0.05

################################################################################
#Test can correctly simulate software FPs
NumIsoformsToSimulate = 4
filteringThreshold =4
Random.seed!(1234)
MaxUnexpressedDetectedIsoforms = map(x -> x > 0.6, rand(1,1,1,4)) #[0.59, 0.76, 0.57, 0.46]
unexpr_dict = Dict("G1" => 1)

#Test can do negatives
testArr = pickUnexpressedErrors(MaxUnexpressedDetectedIsoforms, NumSimulations,
genesList, numCells, unexpr_dict, NumIsoformsToSimulate, filteringThreshold)
@test testArr[1,1,1,1] == 0

#Test can do negatives
testArr = pickUnexpressedErrors(MaxUnexpressedDetectedIsoforms, NumSimulations,
genesList, numCells, unexpr_dict, 3, filteringThreshold)
@test testArr[1,1,1,1] == 1


#Test can do postives
unexpr_dict = Dict("G1" => 2)
testArr = pickUnexpressedErrors(MaxUnexpressedDetectedIsoforms, NumSimulations,
genesList, numCells, unexpr_dict, NumIsoformsToSimulate, filteringThreshold)
@test testArr[1,1,1,1] == 1

#Test can do random isoform choice correctly


###############################################################################
#Test globalArraySimulation function

Random.seed!(1234)
#randIsoChoiceArr = 0.59
#IsoChoiceArr = 0.1
#randIsoChoiceArr =0.85
#unexprRandArr = all false
unexpr_dict = Dict("G1" => 1)
filteringThreshold = 2
MaxNumIsos = 6

testArr = globalArraySimulation(NumSimulations, genesList, rankingDict,
    pDropoutDict, numCells, NumIsoformsToSimulate, filteringThreshold,
    unexpr_dict, MaxNumIsos)

@test testArr == [1]

################################################################################
#Test big simulation function
cumulativeFreqs = [0.2, 1.0]
rankings = [1,2]
pDropouts = [0.9, 0.1]
numCells = 1
IsoformsPerGenePerCell = 1
test = true
pFP = 0.17
pFN = 0.007

#Test the simplest case
@test Simulate(cumulativeFreqs,rankings, pDropouts, numCells,
IsoformsPerGenePerCell, 2, test) == 1

#Test when two isoforms are expressed
@test Simulate(cumulativeFreqs,rankings, pDropouts, numCells,
2, 2, test) == 1

#Test when we swap cumulativeFreqs
@test Simulate([0.8, 1.0], [2,1], pDropouts, numCells,
IsoformsPerGenePerCell, 2, test) == 1

#Test when we swap cumulativeFreqs
@test Simulate([0.8, 1.0], [2,1], pDropouts, numCells,
2, 2, test) == 1

#Test the simplest case when there is more than 1 cell
@test Simulate(cumulativeFreqs,rankings, [0.9, 0.8], 2,
IsoformsPerGenePerCell, 2, test) == 0.5

#Test when we have 3 isoforms
cumulativeFreqs = [0.4, 0.8, 1.0]
rankings = [1,3,2]
pDropouts = [0.9, 0.1, 0.5]

#Test for 1 isoform out of 3
@test Simulate(cumulativeFreqs,rankings, pDropouts, numCells,
IsoformsPerGenePerCell, 3, test) == 1

#Test for 2 isoforms out of 3
@test Simulate(cumulativeFreqs,rankings, pDropouts, numCells,
2, 3, test) == 1

#Test for 3 isoforms out of 3
@test Simulate(cumulativeFreqs,rankings, pDropouts, numCells,
3, 3, test) == 2

#Test that FP and FN logic works
cumulativeFreqs = [0.2, 1.0]
rankings = [1,2]
pDropouts = [0.9, 0.1]
numCells = 1
IsoformsPerGenePerCell = 1
test = true
pFP = 0.6
pFN = 0.007

@test Simulate(cumulativeFreqs,rankings, pDropouts, numCells,
IsoformsPerGenePerCell, 2, test) == 2

pFP = 0.0
pFN = 0.8

@test Simulate(cumulativeFreqs,rankings, pDropouts, numCells,
IsoformsPerGenePerCell, 2, test) == 0

#Test for when ObsExIndex != index
cumulativeFreqs = [0.25,0.5,0.75, 1.0]
rankings = [1,2,3,4]
pDropouts = [0.0,0.0,0.0,0.0]
IsoformsPerGenePerCell = 3
pFP = 0.0
pFN = 0.0

@test Simulate(cumulativeFreqs,rankings, pDropouts, numCells,
IsoformsPerGenePerCell, 4, test) == 3

#Test that number of isoforms logic works
pFP = 1.0
pFN = 0.0

@test Simulate(cumulativeFreqs,rankings, pDropouts, numCells,
IsoformsPerGenePerCell, 10, test) == 10

pFP = 0.0
pFN = 1.0

@test Simulate(cumulativeFreqs,rankings, pDropouts, numCells,
IsoformsPerGenePerCell, 10, test) == 0
