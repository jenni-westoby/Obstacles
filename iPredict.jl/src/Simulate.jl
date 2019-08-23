using StatsBase
using Random
using Distributions

#Simulation process
#Big function that takes dict of isoforms => [counts] for each gene as input
#TODO

################################################################################
#1.Rank Isoforms. See
#https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005761#sec011
################################################################################

function FindRank(InputDict)

    arr = []

    #find sums of values in InputDict
    for (key, value) in InputDict
        push!(arr, sum(value))
    end

    return ordinalrank(arr)
end

function FindPObs(InputDict)

    arr = []

    for (key, value) in InputDict
        pObs = sum(x -> x>0, value) / length(value)
        push!(arr, pObs)
    end

    return arr
end

function FindHarmonicMean(numIsoforms)

    if numIsoforms ==1
        throw(MethodError("It doesn't make sense to find the Harmonic Mean for
        genes with only 1 isoform!"))
    end

    harmonicMean = 0

    for i in 1:numIsoforms
        harmonicMean += (1/i) * exp(-(1 + i/numIsoforms)^2)
    end

    return harmonicMean
end

function FindMedianFrequencies(numIsoforms, harmonicMean)

    arr = Float64[]

    for i in 1:numIsoforms
        top = exp(-(1 + i/numIsoforms)^2)
        bottom = i * harmonicMean
        push!(arr, top/bottom)
    end

    return arr
end

function ScaleFrequencies(frequencies)

    total = sum(frequencies)
    frequencies = frequencies ./ total
end

function FindCumulativeFrequencies(scaledFreqs)

    arr = Float64[]
    scaledFreqs = reverse(sort(scaledFreqs))

    for i in 1:length(scaledFreqs)
        push!(arr, sum(scaledFreqs[1:i]))
    end

    return arr
end

#Function to pick isoforms based on cumulative frequencies and a random number
function PickRank(cumulativeFreqs, random)

    for i in 1:length(cumulativeFreqs)
        if random < cumulativeFreqs[i]
            return i
            break
        end
    end
end

#Function to rescale frequencies after one isoform has been picked
function RescaleFreqs(cumulativeFreqs, rank)

    if rank!=length(cumulativeFreqs)
        for i in (rank + 1):(length(cumulativeFreqs))
            cumulativeFreqs[i] -= cumulativeFreqs[rank]
        end
    end

    deleteat!(cumulativeFreqs, rank)
    cumulativeFreqs = ScaleFrequencies(cumulativeFreqs)
    cumulativeFreqs = FindCumulativeFrequencies(cumulativeFreqs)
end

#Function to rescale rankings after one isoform has been picked
function RescaleRankings(rankings, rank, index)


    deleteat!(rankings, index)

    #small internal rescaling function
    function Rescale(element)
        if element > rank
            return element - 1
        end
        return element
    end

    return Rescale.(rankings)
end

#Wrapper function - only have to do this once for each numIsoforms
function GetIsoformChoiceProbabilities(numIsoforms)

    if numIsoforms > 1
        harmonicMean = FindHarmonicMean(numIsoforms)
        frequencies = FindMedianFrequencies(numIsoforms, harmonicMean)
        frequencies = ScaleFrequencies(frequencies)
        #frequencies = FindCumulativeFrequencies(frequencies)
    else
        frequencies = [1.0]
    end
    return frequencies
end

function RandomIsoformChoiceProbabilities(numIsoforms, numGenes)

    rand_probs = rand(numGenes, numIsoforms)
    output = Array{Float64}(undef, numGenes, numIsoforms)

    for i in 1:numGenes
        rollingSum = sum(rand_probs[i, :])

        for j in 1:numIsoforms
            output[i,j] = rand_probs[i,j] / rollingSum
        end
    end
    return output
end



################################################################################
#2.Find P(dropout) for each isoform
################################################################################
function FindPDropout(InputDict, Km)

    arr = []

    #find PDropout for each value in InputDict
    for (key, value) in InputDict
        S = mean(value)
        PDropout = 1 - (S / (Km + S))
        push!(arr, PDropout)
    end

    return arr
end


function FindMeanExpr(InputDict)
    arr = []

    for (key, value) in InputDict
        S = mean(value)
        push!(arr, S)
    end

    return arr
end

#Introduce quantification errors for expressed isoforms using global pFN and pFP
function QuantErrorsForExpressed(BooleanValue::Bool, randomNumber::Float64)

    if BooleanValue == true
        if randomNumber < pFN
            return false
        end
    else
        if randomNumber < pFP
            return true
        end
    end

    return BooleanValue
end

function RandomQuantErrorsForExpressed(BooleanValue::Bool, randomNumber::Float64)

    if BooleanValue == true
        if randomNumber < rand(Uniform(0,0.5))
            return false
        end
    else
        if randomNumber < rand(Uniform(0,0.5))
            return true
        end
    end

    return BooleanValue
end

function RandomQuantErrorsForUnexpressed(randomNumber::Float64)

    if randomNumber < rand(Uniform(0,0.5))
        return true
    end
    return false
end

################################################################################
# Simulation function
################################################################################

function Simulate(cumulativeFreqs::Array{Float64, 1}, rankings::Array{Int64, 1},
    pDropouts::Array{Float64, 1}, numCells::Int64, IsoformsPerGenePerCell::Int64,
    NumIsoformsForGene::Int64, test::Bool)

    #Initialise some temporary variables
    tmpFreqs = deepcopy(cumulativeFreqs)
    tmpRankings = deepcopy(rankings)
    tmpDropouts = deepcopy(pDropouts)
    staticRankings = deepcopy(rankings)
    ObservedIsoformsPerCell = Float64[]

    #Use a seed when testing so results are reproducible and thus testable
    if test == true
        Random.seed!(1234)
    end

    #Iterate over cells
    for i in 1:numCells

        #Re-initialise some variables
        ObservedExpression = fill(false, NumIsoformsForGene)

        #Necessary to reset temporary variables if IsoformsPerGenePerCell > 1
        if tmpFreqs != cumulativeFreqs
            tmpFreqs = deepcopy(cumulativeFreqs)
            tmpRankings = deepcopy(rankings)
            tmpDropouts = deepcopy(pDropouts)
            staticRankings = deepcopy(rankings)
        end


        #Iterate over the no. expressed isoforms per cell
        for isoform in 1:IsoformsPerGenePerCell

            #Pick isoforms
            choice = PickRank(tmpFreqs, rand())

            #Make two indexes. index is for use with rescaled (and shrinking)
            #temporary variable arrays used in this loop
            index = findall(x -> x==choice, tmpRankings)[1]
            #ObsExIndex is used for indexing into ObservedExpression and needs
            #to not shrink as the temporary variable arrays shrink
            ObsExIndex = staticRankings[index]

            #See whether random number < Pdropout
            if rand() > tmpDropouts[index] #Note if this is not true, ObservedExpression[index] stays as false, which is the default
                ObservedExpression[ObsExIndex] = true
            end

            #Rescaling necessary before next iteration
            if isoform != IsoformsPerGenePerCell
                tmpFreqs = RescaleFreqs(tmpFreqs, choice)
                tmpRankings = RescaleRankings(tmpRankings, choice, index)
                tmpDropouts = deleteat!(tmpDropouts, index)
                staticRankings = deleteat!(staticRankings, index)
            end

        end

        #Simulate quant errors
        for obs in 1:length(ObservedExpression)
            if ObservedExpression[obs] == false
                if rand() < pFP
                     ObservedExpression[obs] = true
                end
            else
                if rand() < pFN
                    ObservedExpression[obs] = false
                end
            end
        end

        #Save the observed num isoforms per cell
        push!(ObservedIsoformsPerCell, sum(ObservedExpression))
    end

    #Return the mean
    return(mean(ObservedIsoformsPerCell))
end

#For GPU version should use curand from CuArrays.CURAND to make random number array before starting
#Or maybe implement a GPU version of rand

################################################################################
#globalArraySimulation() helper functions

#Function that returns 4D array of pDropouts corresponding to isoform choice
function pickRandIsoforms( IsoChoiceProbabilities, NumSimulations,
    genesList, pDropoutDict, numCells, NumIsoformsToSimulate)

    numGenes = length(genesList) #convenient binding
    IsoChoiceArr = Array{Float64}(undef, NumSimulations, numGenes,numCells,NumIsoformsToSimulate)

    for sim in 1:NumSimulations

        for i in 1:length(genesList)

            #get dropout info
            dropoutArr::Array{Float64, 1} = pDropoutDict[genesList[i]]

            #for cell in cells
            for j in 1:numCells

                choices = sample(collect(1:length(IsoChoiceProbabilities[i, :])),
                ProbabilityWeights(IsoChoiceProbabilities[i, :]),
                NumIsoformsToSimulate, replace = false)

                for k in 1:length(choices)
                    IsoChoiceArr[sim,i,j,k] = dropoutArr[choices[k]]
                end
            end
        end
    end
    return IsoChoiceArr
end

function PickUnifObsIsoforms(numIsoforms, genesList, pObsDict,
     pDropoutDict, numCells, NumSimulations)

     numGenes = length(genesList)
     output = Array{Float64}(undef, NumSimulations, numGenes,numCells,numIsoforms)

     for sim in 1:NumSimulations

         for i in 1:numGenes

             #get dropout info and get pObs info
             dropoutArr::Array{Float64, 1} = pDropoutDict[genesList[i]]
             pObsArr::Array{Float64} = pObsDict[genesList[i]]

             #Get pChoices
             pChoiceArr::Array{Float64} = []
             for iso in 1:length(pObsArr)
                 pChoice = pObsArr[iso] / (((1-dropoutArr[iso]) * (1 - pFN)) + pFP)

                 #if pChoice >1, set to 1 as pChoice is a probability
                 #I think this sometimes happens for stochastic reasons or
                 #because we have estimated parameters poorly.
                 if pChoice > 1
                     pChoice = 1.0
                end
                 push!(pChoiceArr, pChoice)
             end

             for j in 1:numCells

                 choices = sample(collect(1:length(pChoiceArr)),
                 ProbabilityWeights(pChoiceArr),
                 numIsoforms, replace = false)

                 for k in 1:length(choices)
                     output[sim,i,j,k] = dropoutArr[choices[k]]
                end
            end
        end
    end
    return output
end

function PickCellVarIsoforms(numIsoforms, genesList, pObsDict,
     pDropoutDict, numCells, NumSimulations)

     numGenes = length(genesList)
     output = Array{Float64}(undef, NumSimulations, numGenes,numCells,numIsoforms)
     sigma = 0.002

     for sim in 1:NumSimulations

         for i in 1:numGenes

             #get dropout info and get pObs info
             dropoutArr::Array{Float64, 1} = pDropoutDict[genesList[i]]
             pObsArr::Array{Float64} = pObsDict[genesList[i]]

             #Get pChoices
             pChoiceArr::Array{Float64} = []
             for iso in 1:length(pObsArr)
                 pChoice = pObsArr[iso] / (((1-dropoutArr[iso]) * (1 - pFN)) + pFP)

                 #if pChoice >1, set to 1 as pChoice is a probability
                 #I think this sometimes happens for stochastic reasons or
                 #because we have estimated parameters poorly.
                 if pChoice > 1
                     pChoice = 1.0
                end
                 push!(pChoiceArr, pChoice)
             end

             for j in 1:numCells

                 cellpChoiceArr::Array{Float64} = []

                 for iso in 1:length(pChoiceArr)

                     mu = pChoiceArr[iso]

                     #This condition ensures alpha and beta are positive.
                     #If it is not fulfilled the beta distrubution is not
                     #appropriate, so we just use the original probability.
                     if mu > 0 && mu < 0.99999

                         #find alpha and beta
                         #@show mu
                         alpha = (((1 - mu)/(sigma^2)) - (1/mu)) * mu^2
                         beta = alpha * ((1/mu) - 1)
                         #@show alpha
                         #@show beta

                         #find p(choice) for this cell by pulling from beta distribution
                         CellChoice = rand(Beta(alpha, beta), 1)
                         push!(cellpChoiceArr, CellChoice[1])
                    else
                        push!(cellpChoiceArr, mu)
                    end
                end

                #@show pObsArr
                #@show pChoiceArr
                #@show cellpChoiceArr
                #@show numIsoforms


                 choices = sample(collect(1:length(cellpChoiceArr)),
                 ProbabilityWeights(cellpChoiceArr),
                 numIsoforms, replace = false)

                 for k in 1:length(choices)
                     output[sim,i,j,k] = dropoutArr[choices[k]]
                end
            end
        end
    end
    return output
end

function PickEqualProbabilities(numIsoforms, genesList, pObsDict,
     pDropoutDict, numCells, NumSimulations, nameDistribution)

     numGenes = length(genesList)
     output = Array{Float64}(undef, NumSimulations, numGenes,numCells,numIsoforms)
     sigma = 0.002

     for sim in 1:NumSimulations

         for i in 1:numGenes

             #get dropout info and get pObs info
             dropoutArr::Array{Float64, 1} = pDropoutDict[genesList[i]]

             #Get pChoices
             pChoiceArr::Array{Float64} = [0.25,0.25,0.25,0.25] #Equal probabilities

             for j in 1:numCells

                 cellpChoiceArr::Array{Float64} = []

                 for iso in 1:length(pChoiceArr)

                     mu = pChoiceArr[iso]

                     #This condition ensures alpha and beta are positive.
                     #If it is not fulfilled the beta distrubution is not
                     #appropriate, so we just use the original probability.
                     if nameDistribution == "Bernoulli"

                         #find p(choice) for this cell by pulling from beta distribution
                         CellChoice = rand(Bernoulli(0.25))
                         push!(cellpChoiceArr, Float64(CellChoice))
                    end

                    if nameDistribution == "Normal"
                        push!(cellpChoiceArr, rand(TruncatedNormal(mu, 0.06, 0.0, 1.0))) #sd chosen to get a large range but >0
                    end
                    if nameDistribution == "Uniform"
                        push!(cellpChoiceArr, mu)
                    end
                end

                #@show pObsArr
                #@show pChoiceArr
                #@show cellpChoiceArr
                #@show numIsoforms

                choices::Array{Int64} = []

                if nameDistribution == "Bernoulli" && sum(cellpChoiceArr) == 1
                    choices = findall(x -> x==1, cellpChoiceArr)
                elseif nameDistribution == "Bernoulli" && sum(cellpChoiceArr) == 0
                    cellpChoiceArr = [0.25,0.25,0.25,0.25]  #collapse back to uniform
                    choices = sample(collect(1:length(cellpChoiceArr)),
                    ProbabilityWeights(cellpChoiceArr),
                    numIsoforms, replace = false)
                else
                    choices = sample(collect(1:length(cellpChoiceArr)),
                    ProbabilityWeights(cellpChoiceArr),
                    numIsoforms, replace = false)
                end

                 for k in 1:length(choices)
                     output[sim,i,j,k] = dropoutArr[choices[k]]
                end
            end
        end
    end
    return output
end

#Function that returns 4D array of pDropouts corresponding to isoform choice
function pickIsoforms(IsoChoiceProbabilities, NumSimulations,
    genesList, rankingDict, pDropoutDict, numCells,
    NumIsoformsToSimulate)

    numGenes = length(genesList) #convenient binding
    IsoChoiceArr = Array{Float64}(undef, NumSimulations, numGenes,numCells,NumIsoformsToSimulate)

    #for gene in genes
    for sim in 1:NumSimulations

        for i in 1:length(genesList)

            #get ranking info
            rankArr::Array{Int64, 1} = rankingDict[genesList[i]]
            dropoutArr::Array{Float64, 1} = pDropoutDict[genesList[i]]

            #for cell in cells
            for j in 1:numCells

                #Sample NumIsoformsToSimulate from IsoChoiceProbabilities
                #@show IsoChoiceProbabilities
                ranks = sample(collect(1:length(IsoChoiceProbabilities)),
                ProbabilityWeights(IsoChoiceProbabilities),
                NumIsoformsToSimulate, replace = false)

                #Find the index of the isoform at the chosen ranks, then write
                #p(dropout) to the appropriate location in the array
                for k in 1:length(ranks)
                    index = findall(x -> x==ranks[k], rankArr)[1]
                    IsoChoiceArr[sim,i,j,k] = dropoutArr[index]
                end
            end
        end
    end

    return IsoChoiceArr
end

#Make a modified version of this for IsoformDependence model
function pickUnexpressedErrors(MaxUnexpressedDetectedIsoforms, NumSimulations,
    genesList, numCells, unexpr_dict, NumIsoformsToSimulate, filteringThreshold)

    numGenes = length(genesList) #convenient binding
    UnexpressedDetectedIsoforms = Array{Int64}(undef, NumSimulations, numGenes,
    numCells)

    #Extract correct number of isoforms per gene
    for sim in 1:NumSimulations

        for i in 1:length(genesList)
            numUnexpr = unexpr_dict[genesList[i]] +
            (filteringThreshold - NumIsoformsToSimulate)

            for j in 1:numCells

                numDetected=0

                for k in 1:numUnexpr
                    numDetected += MaxUnexpressedDetectedIsoforms[sim,i,j,k]
                end

                UnexpressedDetectedIsoforms[sim,i,j] = numDetected
            end
        end
    end

    return UnexpressedDetectedIsoforms
end

#pickUnexpressedErrors for IsoformDependence model
function UnexpressedErrorsIsoformDependence(unexprRandArr, NumSimulations,
    genesList, numCells, unexpr_dict, NumIsoformsToSimulate, filteringThreshold)

    numGenes = length(genesList) #convenient binding
    UnexpressedDetectedIsoforms = Array{Int64}(undef, NumSimulations, numGenes,
    numCells)

    for sim in 1:NumSimulations

        for i in 1:length(genesList)

            #set error rate
            numUnexpr = unexpr_dict[genesList[i]] +
            (filteringThreshold - NumIsoformsToSimulate)
            pFP = ((numUnexpr + NumIsoformsToSimulate) * 0.01)

            for j in 1:numCells

                numDetected = 0

                for k in 1:numUnexpr
                    isDetected = unexprRandArr[sim,i,j,k] < pFP
                    numDetected += isDetected
                end

                UnexpressedDetectedIsoforms[sim,i,j] = numDetected
            end
        end
    end

    return UnexpressedDetectedIsoforms
end

function ExpressedErrorsIsoformDependence(readCapturedIsoforms, randIsoChoiceArr,
    NumSimulations, genesList, numCells, unexpr_dict, NumIsoformsToSimulate,
    filteringThreshold)

    numGenes = length(genesList) #convenient binding
    expressedDetectedIsoforms = Array{Bool}(undef, NumSimulations, numGenes,
    numCells, NumIsoformsToSimulate)
    pFN = 0.04

    for sim in 1:NumSimulations

        for i in 1:numGenes

            #set error rate
            numUnexpr = unexpr_dict[genesList[i]] +
            (filteringThreshold - NumIsoformsToSimulate)
            pFP = ((numUnexpr + NumIsoformsToSimulate) * 0.01)

            for j in 1:numCells
                for k in 1:NumIsoformsToSimulate
                    if readCapturedIsoforms[sim,i,j,k] == true
                        if randIsoChoiceArr[sim,i,j,k] < pFN
                            expressedDetectedIsoforms[sim,i,j,k] = false
                        else
                            expressedDetectedIsoforms[sim,i,j,k] = true
                        end
                    else
                        if randIsoChoiceArr[sim,i,j,k] < pFP
                            expressedDetectedIsoforms[sim,i,j,k] = true
                        else
                            expressedDetectedIsoforms[sim,i,j,k] = false
                        end
                    end
                end
            end
        end
    end
    return expressedDetectedIsoforms
end




#Global simulation function using arrays (hopefully faster)
function globalArraySimulation(NumSimulations, genesList, rankingDict, pObsDict,
    pDropoutDict, numCells, NumIsoformsToSimulate, filteringThreshold,
    unexpr_dict, MaxNumIsos, error_model, choice_model)

    #@show choice_model

    #make a 3D array of dims numGenes by numCells by NumIsoforms and fill with
    #random numbers
    numGenes = length(genesList) #convenient binding
    randIsoChoiceArr = rand(NumSimulations, numGenes,numCells,NumIsoformsToSimulate)
    #@show size(randIsoChoiceArr)

    local WeibullIsoChoiceProbabilities = Array{Float64}(undef, NumIsoformsToSimulate)
    local IsoChoiceProbabilities = Array{Float64}(undef, numGenes, NumIsoformsToSimulate)
    local IsoChoiceArr = Array{Float64}(undef, NumSimulations, numGenes,
    numCells, NumIsoformsToSimulate)

    if choice_model == "Weibull"
        #Get general isoform choice probabilities
        WeibullIsoChoiceProbabilities = GetIsoformChoiceProbabilities(filteringThreshold)

        #Pick isoforms
        IsoChoiceArr = pickIsoforms(WeibullIsoChoiceProbabilities,
        NumSimulations, genesList, rankingDict, pDropoutDict, numCells,
        NumIsoformsToSimulate)

    elseif choice_model =="Random"

        #Get probabilities and pick isoforms
        IsoChoiceProbabilities = RandomIsoformChoiceProbabilities(
        filteringThreshold, length(genesList))

        IsoChoiceArr = pickRandIsoforms(IsoChoiceProbabilities,
        NumSimulations, genesList, pDropoutDict, numCells, NumIsoformsToSimulate)

    elseif choice_model == "UniformObserved"

        IsoChoiceArr = PickUnifObsIsoforms(NumIsoformsToSimulate, genesList,
        pObsDict, pDropoutDict, numCells, NumSimulations)

    elseif choice_model == "CellVariable"

        IsoChoiceArr = PickCellVarIsoforms(NumIsoformsToSimulate, genesList,
        pObsDict, pDropoutDict, numCells, NumSimulations)

    elseif choice_model == "Uniform"

        IsoChoiceArr = PickEqualProbabilities(NumIsoformsToSimulate, genesList,
        pObsDict, pDropoutDict, numCells, NumSimulations, "Uniform")

    elseif choice_model == "Normal"

        IsoChoiceArr = PickEqualProbabilities(NumIsoformsToSimulate, genesList,
        pObsDict, pDropoutDict, numCells, NumSimulations, "Normal")

    elseif choice_model == "Bernoulli"

        IsoChoiceArr = PickEqualProbabilities(NumIsoformsToSimulate, genesList,
        pObsDict, pDropoutDict, numCells, NumSimulations, "Bernoulli")

    else
        throw(ArgumentError("choice_model not recognised."))
    end

    #Do a more than comparison to eliminate dropouts
    #This line returns false when the random number is less than the probability of dropout
    readCapturedIsoforms::Array{Bool} = randIsoChoiceArr .> IsoChoiceArr

    #Fill with new random numbers
    randIsoChoiceArr = rand(NumSimulations, numGenes,numCells,NumIsoformsToSimulate)

    local expressedDetectedIsoforms::Array{Bool}

    #Quant errors for expressed isoforms
    if error_model == "none"
        expressedDetectedIsoforms = map(QuantErrorsForExpressed,
        readCapturedIsoforms, randIsoChoiceArr)
    elseif error_model == "random"
        expressedDetectedIsoforms = map(RandomQuantErrorsForExpressed,
        readCapturedIsoforms, randIsoChoiceArr)
    elseif error_model == "IsoformDependence"
        expressedDetectedIsoforms = ExpressedErrorsIsoformDependence(
        readCapturedIsoforms, randIsoChoiceArr, NumSimulations, genesList,
        numCells, unexpr_dict, NumIsoformsToSimulate, filteringThreshold)
    else
        throw(ArgumentError("Unrecognised error_model"))
    end

    #Find number of expressed detected isoforms per gene per cell
    totalExpressedDetected = reshape(sum(expressedDetectedIsoforms, dims = 4),
    NumSimulations, numGenes, numCells)

    #Now we need to do quant errors for unexpressed isoforms
    #First, make the biggest array of random numbers we could possibly need
    unexprRandArr = rand(NumSimulations, numGenes,numCells,
    MaxNumIsos - NumIsoformsToSimulate) #Do something about hard coded numbers

    local MaxUnexpressedDetectedIsoforms::Array{Bool}

    #Apply quant errors
    if error_model == "none"
        MaxUnexpressedDetectedIsoforms = map(x -> x < pFP,
        unexprRandArr)
        UnexpressedDetectedIsoforms = pickUnexpressedErrors(
        MaxUnexpressedDetectedIsoforms, NumSimulations, genesList, numCells,
        unexpr_dict, NumIsoformsToSimulate, filteringThreshold)
    elseif error_model == "random"
        MaxUnexpressedDetectedIsoforms = map(
        RandomQuantErrorsForUnexpressed, unexprRandArr)
        UnexpressedDetectedIsoforms = pickUnexpressedErrors(
        MaxUnexpressedDetectedIsoforms, NumSimulations, genesList, numCells,
        unexpr_dict, NumIsoformsToSimulate, filteringThreshold)
    elseif error_model == "IsoformDependence"
        UnexpressedDetectedIsoforms = UnexpressedErrorsIsoformDependence(
        unexprRandArr, NumSimulations, genesList, numCells, unexpr_dict,
        NumIsoformsToSimulate, filteringThreshold)
    else
        throw(ArgumentError("Unrecognised error_model"))
    end

    #Find total num detected isoforms and convert to array of means
    DetectedIsoforms = map(+, totalExpressedDetected, UnexpressedDetectedIsoforms)

    MeanIsoformsPerGenePerCell = reshape(mean(DetectedIsoforms, dims = 3),
    NumSimulations * numGenes)

    return MeanIsoformsPerGenePerCell
end
