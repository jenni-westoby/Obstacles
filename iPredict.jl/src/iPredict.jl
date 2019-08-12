module iPredict
include("TableProcessing.jl")
include("FindKm.jl")
include("Simulate.jl")
include("Plots.jl")
include("pValues.jl")
import Cairo
import Fontconfig
import DataFrames
using Distributions
using Statistics
using DelimitedFiles



#Function that performs global level predictions
function globalPredict(countsMatrixPath::String, GeneIsoTablePath::String,
     NumSimulations::Int64, output::String, NumIsoformsToSimulate::Int64,
     filteringThreshold::Int64, test::Bool, alpha, beta, pFNeg::Float64,
     pFPos::Float64, error_model::String, choice_model::String, path::String)

    #Set probability(FP|N) and probability(FN|P)
    if error_model == "none"
        global pFN = pFNeg
        global pFP = pFPos
    end

    #Read in counts matrix
    println("Data wrangling commencing!")
    println("Reading in counts matrix.") #Quick
    countsMatrix = ReadCountsMatrix(countsMatrixPath)

    #Read in gene-isoform table
    println("Reading in gene-isoform relationships")
    GeneIsoformRelationships = ReadGeneIsoformRelationships(GeneIsoTablePath)

    #Make sure gene-isoform table matches counts matrix
    if CheckGenesMatchIsoforms(countsMatrix, GeneIsoformRelationships) == false
        throw(ArgumentError("Your counts matrix doesn't match your list of genes
        and isoforms - there is something wrong with your input files."))
    end

    #Add gene information to isoforms table
    println("Adding Gene information to countsMatrix")
    countsMatrixWithGenes = AddGeneInfo(countsMatrix, GeneIsoformRelationships)

    #Only find Km if no Beta distribution parameters provided
    if isnothing(alpha) || isnothing(beta)

        #Clean data for finding Km
        println("Cleaning data for finding Km")
        MMdf = CleanDataForMM(countsMatrix)

        #Find Km
        println("Calculating Km....")
        MMParams = FindKm(MMdf)
        Km = MMParams[1]
        sigma = MMParams[2]
    end

    #Extract genes with filteringThreshold expressed isoforms
    println("Converting counts to dict")
    GenesWithNIsoforms = KeepGenesWithNIsoforms(countsMatrixWithGenes,
    filteringThreshold)
    workingDict = GenesWithNIsoforms[1]
    unexpr_dict = GenesWithNIsoforms[2]
    MaxNumIsos = GenesWithNIsoforms[3]

    #Get ranking and dropout info for each gene
    rankingDict = Dict()
    pObsDict = Dict()
    pDropoutDict = Dict()
    genesList = []

    #How we calculate dropouts depends on whether or not beta parameters are provided
    if isnothing(alpha) || isnothing(beta)
        for (key, value) in workingDict
            rankingDict[key] = FindRank(value)
            pObsDict[key] = FindPObs(value)
            pDropoutDict[key] = FindPDropout(value, Km)
            push!(genesList, key)
        end
    else
        for (key, value) in workingDict
            rankingDict[key] = FindRank(value)
            pObsDict[key] = FindPObs(value)
            pDropoutDict[key] = rand(Beta(alpha, beta), length(value))
            push!(genesList, key)
        end
    end

    plot_arr = []
    log_plot_arr = []
    mean_arr = []
    median_arr = []
    std_arr = []

    println("Beginning simulations!")

    for i in 1:NumIsoformsToSimulate

        #I pre-allocate arr in the hope it will make my code faster
        arr = Array{Float64}(undef, NumSimulations * length(genesList))
        arr_index = 1

        #Need to add a for loop with num simulations here
        for j in 1:NumSimulations

            #@show arr_index

            #Get an array of mean isoforms per gene per cell
            arr[arr_index:(arr_index + length(genesList) - 1)] =
            globalArraySimulation(1, genesList, rankingDict, pObsDict,
            pDropoutDict, ncol(countsMatrix) - 1, i, filteringThreshold,
            unexpr_dict, MaxNumIsos, error_model, choice_model)

            #Increment arr_index
            arr_index += length(genesList)
        end

        #Save histogram of data
        if i == 1
            p1 = makeNormalisedHistogram(arr, string(i) * " Isoform")
        else
            p1 = makeNormalisedHistogram(arr, string(i) * " Isoforms")
        end
        push!(plot_arr, p1[1])
        push!(log_plot_arr, p1[2])
        push!(mean_arr, mean(arr))
        push!(std_arr, std(arr))
        push!(median_arr, median(arr))

        #Write to file if required
        if path != "none"
            writedlm(path * "_" * string(i) * "_isoforms_array.txt", arr)
        end
    end

    #Get an array and a dict of mean isoforms per gene per cell in real data
    println("Finding means in real data.")
    findRealMeans = MeanNumIsoformsPerGenePerCellReal(workingDict)
    realMeansArr = findRealMeans[1]

    #Write to file if required
    if path != "none"
        writedlm(path * "_" * "real_isoforms_array.txt", realMeansArr)
    end

    #Plot real data
    real_plot = makeNormalisedHistogramReal(realMeansArr, "Real")
    push!(plot_arr, real_plot[1])
    push!(log_plot_arr, real_plot[2])

    #Arrange plots and write to output
    graph = plot(plot_arr..., layout=(length(plot_arr),1),legend=false, size = (300, 600), bottom_margin = 1.0mm)
    savefig(graph, output)
    graph = plot(log_plot_arr..., layout=(length(plot_arr),1),legend=false, size = (300, 600))
    savefig(graph, output * "_log.png")

    return (mean_arr, std_arr, median_arr)


end

function oneGenePrediction(countsMatrixPath::String, GeneIsoTablePath::String,
    gene::String, NumSimulations::Int64, output::String, test::Bool)

    #Read in counts matrix
    println("Data wrangling commencing!")
    println("Reading in counts matrix.")
    countsMatrix = ReadCountsMatrix(countsMatrixPath)

    #Read in gene-isoform table
    println("Reading in gene-isoform relationships")
    GeneIsoformRelationships = ReadGeneIsoformRelationships(GeneIsoTablePath)

    #Make sure gene-isoform table matches counts matrix
    if CheckGenesMatchIsoforms(countsMatrix, GeneIsoformRelationships) == false
        throw(ArgumentError("Your counts matrix doesn't match your list of genes
        and isoforms - there is something wrong with your input files."))
    end

    #Add gene information to isoforms table
    println("Adding Gene information to countsMatrix")
    countsMatrixWithGenes = AddGeneInfo(countsMatrix, GeneIsoformRelationships)

    #Make a dictionary just for the expressed isoforms in gene - no point
    #simulating unexpressed isoforms.
    tmp = filter(r -> any(occursin.([gene], r.Genes)), countsMatrixWithGenes)
    NumIsos = nrow(tmp)  #Get number of isoforms (incl. unexpressed) for gene of interest
    countsDict = ConvertExpressedIsosToDict(tmp)

    #Clean data for finding Km
    println("Cleaning data for finding Km")
    MMdf = CleanDataForMM(countsMatrix)

    #Find Km
    println("Calculating Km....")
    MMParams = FindKm(MMdf)
    Km = MMParams[1]
    sigma = MMParams[2]

    println("Preparing for simulations...")

    #Find ranking info
    cumulativeFreqs = GetIsoformChoiceProbabilities(length(countsDict[gene]))

    meansDict = Dict()

    #For each gene - bit unnecessary since just 1 gene, but convenient way of
    #getting value.
    for (key, value) in countsDict

        #Rank isoforms
        rankArr::Array{Int64, 1} = FindRank(value)

        #Find p(dropouts)
        dropoutArr::Array{Float64, 1} = FindPDropout(value, Km)

        for i in 1:length(countsDict[gene])

            meansArr = []

            for j in 1:NumSimulations
                mean = Simulate(cumulativeFreqs, rankArr, dropoutArr,
                ncol(countsMatrix) - 1, i, NumIsos, test)
                push!(meansArr, mean)
            end

            meansDict[i] = meansArr

        end
    end

    #Make a plot and write it to file
    OneGeneHistPlot(meansDict, tmp, output)

    #Find and return p-values
    mean = MeanIsoPerGenePerCell(tmp)
    return FindpValues(meansDict, mean, 40) #40 is an arbitrary number of required unique values

end

export globalPredict
export oneGenePrediction

end #end of module

#GPU structure of kernel
#find p(dropout)
#for i in 1:10000
#simulation process
#save average
#return an array of averages
