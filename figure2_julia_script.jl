include("iPredict.jl/src/iPredict.jl")
using .iPredict
using DataFrames
using Plots, Plots.Measures
using JLD2
using DelimitedFiles

#Initialise variables
GeneIsoformRelationships = "data/gencode_human_gene_isoform_relationships.txt"
NumSimulations = 100
NumIsoformsToSimulate = 4
filteringThreshold = 4
alpha_vals = [0.45, 0.7, 0.95, 1.2, 1.45]
beta_vals = reverse(alpha_vals)
scnorm_counts = ["H1_24_scnorm_counts.txt",
"H1_96_scnorm_counts.txt",
"H9_24_scnorm_counts.txt",
"H9_96_scnorm_counts.txt"]
pFP = 0.01
pFN = 0.04

df = DataFrame(counts=String[], alpha = Float64[], beta = Float64[],
NumIsos = Int64[], mean = Float64[], std = Float64[], median = Float64[])

for countsMatrix in scnorm_counts

    output = "figures/figure2/" * countsMatrix
    data_out = "figures/figure2_data/" * countsMatrix

    globalPredict("data/" * countsMatrix, GeneIsoformRelationships,
    NumSimulations, output, NumIsoformsToSimulate, filteringThreshold,
    false, nothing,nothing, pFN, pFP, "none", "Weibull", data_out)

    for i in 1:length(alpha_vals)

        output = "figures/figure2/" * string(alpha_vals[i]) * "_" *
        string(beta_vals[i]) * "_" * countsMatrix
        data_out = "figures/figure2_data/" * string(alpha_vals[i]) * "_" *
        string(beta_vals[i]) * "_" * countsMatrix

        out = globalPredict("data/" * countsMatrix, GeneIsoformRelationships,
        NumSimulations, output, NumIsoformsToSimulate, filteringThreshold,
        false, alpha_vals[i], beta_vals[i], pFN, pFP, "none", "Weibull",
        data_out)

        newRow = DataFrame(counts = repeat([countsMatrix], 4),
        alpha = repeat([alpha_vals[i]], 4),
        beta = repeat([beta_vals[i]], 4),
        NumIsos = [1,2,3,4], mean = out[1], std = out[2], median = out[3])

        append!(df, newRow)
    end
end

plot_arr = []

jldopen("figures/figure2_data/means_df.jld2", "w") do file
    file["df"] = df
end

for i in 1:4
    tmp = df[df.NumIsos .== i, :]
    both_plot = Plots.plot(tmp.alpha, tmp.mean, group = tmp.counts, title = string(i))
    push!(plot_arr, both_plot)
end

p1 = plot(plot_arr..., legend = :none, xlabel = "alpha", ylabel = "mean")
savefig(p1, "figures/figure2/mean_scnorm_alpha.png")

################################################################################
# Make figure 2
################################################################################

function make_hist(path, titleArg)
    Arr = readdlm(path)

    if (titleArg == "3 Isoforms" &&
        path == "figures/figure2_data/H1_96_scnorm_counts.txt_3_isoforms_array.txt")
        #Untransformed plot
        p1 = histogram(Arr,
        title = titleArg,
        xlims = (0, 4),
        ylims = (0, 10),
        ylabel = "Normalised Frequency",
        normalize=:pdf,
        bins = 100,
        axis = (font(100)), titlefontsize=100, yticks = 0:5:10,
        margin = 25mm)

    elseif titleArg == "Real" || titleArg == "alpha=0.45, beta=1.45"
        p1 = histogram(Arr,
        title = titleArg,
        xlims = (0, 4),
        ylims = (0, 10),
        normalize=:pdf,
        xlabel = "Mean Isoforms per Gene per Cell",
        bins = 100,
        axis = (font(100)), titlefontsize=100, yticks = 0:5:10,
        margin = 25mm)

    else
        #Untransformed plot
        p1 = histogram(Arr,
        title = titleArg,
        xlims = (0, 4),
        ylims = (0, 10),
        normalize=:pdf,
        bins = 100,
        axis = (font(100)), titlefontsize=100, yticks = 0:5:10,
        margin = 25mm)
    end
    return p1
end

H1_24_1 = make_hist("figures/figure2_data/H1_24_scnorm_counts.txt_1_isoforms_array.txt", "1 Isoform")
H1_24_2 = make_hist("figures/figure2_data/H1_24_scnorm_counts.txt_2_isoforms_array.txt", "2 Isoforms")
H1_24_3 = make_hist("figures/figure2_data/H1_24_scnorm_counts.txt_3_isoforms_array.txt", "3 Isoforms")
H1_24_4 = make_hist("figures/figure2_data/H1_24_scnorm_counts.txt_4_isoforms_array.txt", "4 Isoforms")
H1_24_real = make_hist("figures/figure2_data/H1_24_scnorm_counts.txt_real_isoforms_array.txt", "Real")

H1_96_1 = make_hist("figures/figure2_data/H1_96_scnorm_counts.txt_1_isoforms_array.txt", "1 Isoform")
H1_96_2 = make_hist("figures/figure2_data/H1_96_scnorm_counts.txt_2_isoforms_array.txt", "2 Isoforms")
H1_96_3 = make_hist("figures/figure2_data/H1_96_scnorm_counts.txt_3_isoforms_array.txt", "3 Isoforms")
H1_96_4 = make_hist("figures/figure2_data/H1_96_scnorm_counts.txt_4_isoforms_array.txt", "4 Isoforms")
H1_96_real = make_hist("figures/figure2_data/H1_96_scnorm_counts.txt_real_isoforms_array.txt", "Real")

pA = plot(H1_96_1,H1_24_1,
H1_96_2,H1_24_2,
H1_96_3,H1_24_3,
H1_96_4,H1_24_4,
H1_96_real,H1_24_real,
layout = (5,2), legend = false, size = (500*8, 600*8))

savefig(pA, "figures/figure2/figure2A.png")

################################################################################
# Panel B
################################################################################

include("iPredict.jl/src/TableProcessing.jl")
include("iPredict.jl/src/FindKm.jl")
include("iPredict.jl/src/Simulate.jl")
include("iPredict.jl/src/Plots.jl")
include("iPredict.jl/src/pValues.jl")
import Cairo
import Fontconfig
import DataFrames
using Distributions
using DelimitedFiles

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
    for (key, value) in workingDict
        append!(pDropoutArr, FindPDropout(value, Km))
    end

    pDropoutDistributionPlot = pDropoutDistribution(pDropoutArr, output)
    params = fit(Beta, pDropoutArr)
    return (params.α, params.β, pDropoutDistributionPlot[1], pDropoutDistributionPlot[2])
end

H9_96 = pDropoutInfo("data/H9_96_scnorm_counts.txt", "data/gencode_human_gene_isoform_relationships.txt", "figures/pdrop_scnorm_H9_96.png", 4)
H9_24 = pDropoutInfo("data/H9_24_scnorm_counts.txt", "data/gencode_human_gene_isoform_relationships.txt", "figures/pdrop_scnorm_H9_24.png", 4)
H1_96 = pDropoutInfo("data/H1_96_scnorm_counts.txt", "data/gencode_human_gene_isoform_relationships.txt", "figures/pdrop_scnorm_H1_96.png", 4)
H1_24 = pDropoutInfo("data/H1_24_scnorm_counts.txt", "data/gencode_human_gene_isoform_relationships.txt", "figures/pdrop_scnorm_H1_24.png", 4)

pB = plot(H1_96[3], H1_24[3],
       H1_96[4], H1_24[4],
        size = (500 * 8, 240 * 8), layout = (2,2), legend = false, margin = 35mm)

savefig(pB, "figures/figure2/figure2B.png")

################################################################################
# Panel C
################################################################################

v0 = pdf(Beta(0.45, 1.45), 0:0.01:1)
v1 = pdf(Beta(0.7, 1.2), 0:0.01:1)
v2 = pdf(Beta(0.95, 0.95), 0:0.01:1)
v3 = pdf(Beta(1.2, 0.7), 0:0.01:1)
v4 = pdf(Beta(1.45, 0.45), 0:0.01:1)
pC = plot(0:0.01:1, [v0,v1,v2,v3,v4], label = ["alpha=0.45, beta=1.45",
       "alpha=0.7, beta=1.2",
       "alpha=0.95, beta=0.95",
       "alpha=1.2, beta=0.7",
       "alpha=1.45, beta=0.45"], size = (250 * 8, 240 * 8), axis = (font(60)),
       titlefontsize=80, legendfontsize = 50, linestyle=:solid, linealpha=0.5,
        linewidth=2*8, ylabel = "Normalised Frequency", xlabel = "Probability",
        title = "Beta Distributions", margin = 30mm)
savefig(pC, "figures/figure2/figure2C.png")

################################################################################
# Panel D
################################################################################

v0 =  make_hist("figures/figure2_data/0.45_1.45_H1_24_scnorm_counts.txt_1_isoforms_array.txt", "alpha=0.45, beta=1.45")
v1 =  make_hist("figures/figure2_data/0.7_1.2_H1_24_scnorm_counts.txt_1_isoforms_array.txt", "alpha=0.7, beta=1.2")
v2 =  make_hist("figures/figure2_data/0.95_0.95_H1_24_scnorm_counts.txt_1_isoforms_array.txt", "alpha=0.95, beta=0.95")
v3 =  make_hist("figures/figure2_data/1.2_0.7_H1_24_scnorm_counts.txt_1_isoforms_array.txt", "alpha=1.2, beta=0.7")
v4 =  make_hist("figures/figure2_data/1.45_0.45_H1_24_scnorm_counts.txt_1_isoforms_array.txt", "alpha=1.45, beta=0.45")

pD = plot(v4, v3, v2, v1, v0,
       layout = (5,1), legend = false, size = (250*8, 600*8),
       axis = (font(180)), titlefontsize=160, margin = 30mm)

savefig(pD, "figures/figure2/figure2D.png")

# #Doesn't work sad times
# l = @layout [  grid(7,2){0.66w} [b{0.285h}
#                                          grid(5,1) ]]
#
# plot(H1_96_1,H1_24_1,
# H1_96_2,H1_24_2,
# H1_96_3,H1_24_3,
# H1_96_4,H1_24_4,
# H1_96_real,H1_24_real,
# H1_96[3], H1_24[3],
# H1_96[4], H1_24[4],
# pC,
# v4, v3, v2, v1, v0, layout = l)
