include("iPredict.jl/src/iPredict.jl")
using .iPredict
using DataFrames
using Plots, Plots.Measures
using JLD2
using DelimitedFiles
using StatsBase

function make_hist(path, titleArg)
    Arr = readdlm(path)

    if (titleArg == "3 Isoforms" &&
        path == "figures/figure2_data/H9_96_scnorm_counts.txt_3_isoforms_array.txt")
        #Untransformed plot
        p1 = histogram(Arr,
        title = titleArg,
        xlims = (0, 4),
        ylims = (0, 10),
        ylabel = "Probability Density",
        normalize=:pdf,
        bins = 100,
        axis = (font(100)), titlefontsize=100, yticks = 0:5:10,
        margin = 25mm)

    elseif titleArg == "Real" #|| titleArg == "alpha=0.45, beta=1.45"
        p1 = histogram(Arr,
        title = titleArg,
        xlims = (0, 4),
        ylims = (0, 10),
        normalize=:pdf,
        xlabel = "Mean No. Isoforms\nper Gene per Cell",
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
    vline!([mean(Arr)], color = :black, linewidth = 8)
    return p1
end

H9_24_1 = make_hist("figures/figure2_data/H9_24_scnorm_counts.txt_1_isoforms_array.txt", "1 Isoform")
H9_24_2 = make_hist("figures/figure2_data/H9_24_scnorm_counts.txt_2_isoforms_array.txt", "2 Isoforms")
H9_24_3 = make_hist("figures/figure2_data/H9_24_scnorm_counts.txt_3_isoforms_array.txt", "3 Isoforms")
H9_24_4 = make_hist("figures/figure2_data/H9_24_scnorm_counts.txt_4_isoforms_array.txt", "4 Isoforms")
H9_24_real = make_hist("figures/figure2_data/H9_24_scnorm_counts.txt_real_isoforms_array.txt", "Real")

H9_96_1 = make_hist("figures/figure2_data/H9_96_scnorm_counts.txt_1_isoforms_array.txt", "1 Isoform")
H9_96_2 = make_hist("figures/figure2_data/H9_96_scnorm_counts.txt_2_isoforms_array.txt", "2 Isoforms")
H9_96_3 = make_hist("figures/figure2_data/H9_96_scnorm_counts.txt_3_isoforms_array.txt", "3 Isoforms")
H9_96_4 = make_hist("figures/figure2_data/H9_96_scnorm_counts.txt_4_isoforms_array.txt", "4 Isoforms")
H9_96_real = make_hist("figures/figure2_data/H9_96_scnorm_counts.txt_real_isoforms_array.txt", "Real")

pA = plot(H9_96_1,H9_24_1,
H9_96_2,H9_24_2,
H9_96_3,H9_24_3,
H9_96_4,H9_24_4,
H9_96_real,H9_24_real,
layout = (5,2), legend = false, size = (500*8, 600*8),
margin = 45mm)

savefig(pA, "figures/supplementary_figs/figure1A.png")

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

pB = plot(H9_96[3], H9_24[3],
       H9_96[4], H9_24[4],
        size = (500 * 8, 240 * 8), layout = (2,2), legend = false, margin = 35mm)

        @show H9_96[1]
        @show H9_96[2]
        @show H9_24[1]
        @show H9_24[2]

savefig(pB, "figures/supplementary_figs/figure1B.png")

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
       "alpha=1.45, beta=0.45"], size = (250 * 8, 240 * 8), axis = (font(60)), margin = 35mm,
       titlefontsize=80, legendfontsize = 60, linestyle=:solid, linealpha=0.5, linewidth=2*8, ylabel = "Normalised Frequency", xlabel = "Probability", title = "Beta Distributions")
savefig(pC, "figures/supplementary_figs/figure1C.png")

################################################################################
# Panel D
################################################################################

function make_hist_MM(path, titleArg)
    Arr = readdlm(path)

    if (titleArg == "3 Isoforms" &&
        path == "figures/figure2_data/H9_96_scnorm_counts.txt_3_isoforms_array.txt")
        #Untransformed plot
        p1 = histogram(Arr,
        title = titleArg,
        xlims = (0, 1.2),
        ylims = (0, 5),
        ylabel = "Probability Density",
        normalize=:pdf,
        bins = 100,
        axis = (font(100)), titlefontsize=100, yticks = 0:5:10,
        margin = 25mm)

    elseif titleArg == "Real" #|| titleArg == "alpha=0.45, beta=1.45"
        p1 = histogram(Arr,
        title = titleArg,
        xlims = (0, 1.2),
        ylims = (0, 5),
        normalize=:pdf,
        xlabel = "Mean No. Isoforms\nper Gene per Cell",
        bins = 100,
        axis = (font(100)), titlefontsize=100, yticks = 0:5:10,
        margin = 25mm)

    else
        #Untransformed plot
        p1 = histogram(Arr,
        title = titleArg,
        xlims = (0, 1.2),
        ylims = (0, 5),
        normalize=:pdf,
        bins = 100,
        axis = (font(100)), titlefontsize=100, yticks = 0:5:10,
        margin = 25mm)
    end
    vline!([mean(Arr)], color = :black, linewidth = 8)
    return p1
end


v0 =  make_hist_MM("figures/figure2_data/0.45_1.45_H9_24_scnorm_counts.txt_1_isoforms_array.txt", "alpha=0.45, beta=1.45")
v1 =  make_hist_MM("figures/figure2_data/0.7_1.2_H9_24_scnorm_counts.txt_1_isoforms_array.txt", "alpha=0.7, beta=1.2")
v2 =  make_hist_MM("figures/figure2_data/0.95_0.95_H9_24_scnorm_counts.txt_1_isoforms_array.txt", "alpha=0.95, beta=0.95")
v3 =  make_hist_MM("figures/figure2_data/1.2_0.7_H9_24_scnorm_counts.txt_1_isoforms_array.txt", "alpha=1.2, beta=0.7")
v4 =  make_hist_MM("figures/figure2_data/1.45_0.45_H9_24_scnorm_counts.txt_1_isoforms_array.txt", "alpha=1.45, beta=0.45")

pD = plot(v4, v3, v2, v1, v0,
       layout = (5,1), legend = false, size = (250*8, 600*8),
       axis = (font(180)), titlefontsize=160, margin = 25mm)

savefig(pD, "figures/supplementary_figs/figure1D.png")

####################################
# Overlap plots
####################################

function make_overlap_hist(path, titleArg)

    Arr = readdlm(path)

    if titleArg == "4 Isoforms"

        p1 = histogram(Arr,
        title = titleArg,
        color = palette(:default)[2],
        xlims = (0, 1),
        ylims = (0, 10),
        ylabel = "Probability Density",
        xlabel = "% Overlap with Ground Truth",
        normalize=:pdf,
        bins = 100,
        axis = (font(100)), titlefontsize=100, yticks = 0:5:10,
        margin = 25mm)
        vline!([mean(Arr)], color = :black, linewidth = 8)
        return p1

    else

        p1 = histogram(Arr,
        title = titleArg,
        color = palette(:default)[2],
        xlims = (0, 1),
        ylims = (0, 10),
        ylabel = "Probability Density",
        normalize=:pdf,
        bins = 100,
        axis = (font(100)), titlefontsize=100, yticks = 0:5:10,
        margin = 25mm)
        vline!([mean(Arr)], color = :black, linewidth = 8)
        return p1
    end


end

H9_24_1 = make_overlap_hist("figures/figure2_data/H9_24_scnorm_counts.txt_1_isoforms_overlap_array.txt", "1 Isoform")
H9_24_2 = make_overlap_hist("figures/figure2_data/H9_24_scnorm_counts.txt_2_isoforms_overlap_array.txt", "2 Isoforms")
H9_24_3 = make_overlap_hist("figures/figure2_data/H9_24_scnorm_counts.txt_3_isoforms_overlap_array.txt", "3 Isoforms")
H9_24_4 = make_overlap_hist("figures/figure2_data/H9_24_scnorm_counts.txt_4_isoforms_overlap_array.txt", "4 Isoforms")

H9_96_1 = make_overlap_hist("figures/figure2_data/H9_96_scnorm_counts.txt_1_isoforms_overlap_array.txt", "1 Isoform")
H9_96_2 = make_overlap_hist("figures/figure2_data/H9_96_scnorm_counts.txt_2_isoforms_overlap_array.txt", "2 Isoforms")
H9_96_3 = make_overlap_hist("figures/figure2_data/H9_96_scnorm_counts.txt_3_isoforms_overlap_array.txt", "3 Isoforms")
H9_96_4 = make_overlap_hist("figures/figure2_data/H9_96_scnorm_counts.txt_4_isoforms_overlap_array.txt", "4 Isoforms")

pA = plot(H9_96_1,H9_24_1,
H9_96_2,H9_24_2,
H9_96_3,H9_24_3,
H9_96_4,H9_24_4,
layout = (4,2), legend = false, size = (400*8, 600*8))

savefig(pA, "figures/supplementary_figs/figure2_overlap.png")

########################
# Old 1D overlap hists
########################

v0 =  make_overlap_hist("figures/figure2_data/0.45_1.45_H9_24_scnorm_counts.txt_1_isoforms_overlap_array.txt", "alpha=0.45, beta=1.45")
v1 =  make_overlap_hist("figures/figure2_data/0.7_1.2_H9_24_scnorm_counts.txt_1_isoforms_overlap_array.txt", "alpha=0.7, beta=1.2")
v2 =  make_overlap_hist("figures/figure2_data/0.95_0.95_H9_24_scnorm_counts.txt_1_isoforms_overlap_array.txt", "alpha=0.95, beta=0.95")
v3 =  make_overlap_hist("figures/figure2_data/1.2_0.7_H9_24_scnorm_counts.txt_1_isoforms_overlap_array.txt", "alpha=1.2, beta=0.7")
v4 =  make_overlap_hist("figures/figure2_data/1.45_0.45_H9_24_scnorm_counts.txt_1_isoforms_overlap_array.txt", "alpha=1.45, beta=0.45")

pD = plot(v4, v3, v2, v1, v0,
       layout = (5,1), legend = false, size = (250*8, 600*8),
       axis = (font(180)), titlefontsize=160, margin = 30mm)

savefig(pD, "figures/supplementary_figs/figure2D_overlap.png")
