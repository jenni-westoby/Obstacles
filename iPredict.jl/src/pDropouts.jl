include("TableProcessing.jl")
include("FindKm.jl")
include("Simulate.jl")
include("Plots.jl")
include("pValues.jl")
import Cairo
import Fontconfig
import DataFrames
using Distributions
using DelimitedFiles
using HDF5, JLD

#Set probability(FP|N) and probability(FN|P)
pFN = 0.04
pFP = 0.01

#Function that performs global level predictions
function pDropoutInfo(countsMatrixPath::String, GeneIsoTablePath::String,
      output::String, filteringThreshold::Int64)

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

    #Clean data for finding Km
    println("Cleaning data for finding Km")
    MMdf = CleanDataForMM(countsMatrix)

    #Find Km
    println("Calculating Km....")
    MMParams = FindKm(MMdf)
    Km = MMParams[1]
    sigma = MMParams[2]

    println("Converting counts to dict")
    GenesWithNIsoforms = KeepGenesWithNIsoforms(countsMatrixWithGenes,
    filteringThreshold)
    workingDict = GenesWithNIsoforms[1]
    unexpr_dict = GenesWithNIsoforms[2]
    MaxNumIsos = GenesWithNIsoforms[3]

    #Get ranking and dropout info for each gene
    pDropoutArr::Array{Float64} = []
    pDropoutDict = Dict()
    ExprArr::Array{Float64} =[]
    ExprDict = Dict()


    for (key, value) in workingDict
        ExprArr = []
        for (k,v) in value
            S = mean(v)
            push!(ExprArr, S)
        end

        append!(pDropoutArr, FindPDropout(value, Km))
        pDropoutDict[key] = FindPDropout(value, Km)
        ExprDict[key] = ExprArr
    end

    save(output * "_dict.jld", "data", pDropoutDict)
    save(output * "expr_dict.jld", "data", ExprDict)

    pDropoutDistributionPlot = pDropoutDistribution(pDropoutArr, output)
    params = fit(Beta, pDropoutArr)
    return (params.α, params.β, pDropoutDistributionPlot)
end

H9_96 = pDropoutInfo("data/H9_96_scnorm_counts.txt", "data/gencode_human_gene_isoform_relationships.txt", "figures/pdrop_scnorm_H9_96.png", 4)
H9_24 = pDropoutInfo("data/H9_24_scnorm_counts.txt", "data/gencode_human_gene_isoform_relationships.txt", "figures/pdrop_scnorm_H9_24.png", 4)
H1_96 = pDropoutInfo("data/H1_96_scnorm_counts.txt", "data/gencode_human_gene_isoform_relationships.txt", "figures/pdrop_scnorm_H1_96.png", 4)
H1_24 = pDropoutInfo("data/H1_24_scnorm_counts.txt", "data/gencode_human_gene_isoform_relationships.txt", "figures/pdrop_scnorm_H1_24.png", 4)

p1 = plot([H1_24[1], H1_96[1], H9_96[1], H9_24[1]], [H1_24[2], H1_96[2], H9_96[2], H9_24[2]], seriestype=:scatter)
savefig(p1, "figures/scnorm_beta_params.png")

v0 = pdf(Beta(0.45, 1.45), 0:0.01:1)
v1 = pdf(Beta(0.7, 1.2), 0:0.01:1)
v2 = pdf(Beta(0.95, 0.95), 0:0.01:1)
v3 = pdf(Beta(1.2, 0.7), 0:0.01:1)
v4 = pdf(Beta(1.45, 0.45), 0:0.01:1)
p_all = plot(0:0.01:1, [v0,v1,v2,v3,v4], label = ["alpha=0.45, beta=1.45",
"alpha=0.7, beta=1.2",
"alpha=0.95, beta=0.95",
"alpha=1.2, beta=0.7",
"alpha=1.45, beta=0.45"])
savefig(p_all, "figures/proposed_beta_distributions.png")
