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
scnorm_counts = [#"H1_24_scnorm_counts.txt",
"H1_96_scnorm_counts.txt",
"H9_24_scnorm_counts.txt",
"H9_96_scnorm_counts.txt"]
pFP = 0.01
pFN = 0.04

for countsMatrix in scnorm_counts

    output = "figures/negControl/" * countsMatrix
    data_out = "figures/negControl/" * countsMatrix

    #NegativeControl("data/" * countsMatrix, GeneIsoformRelationships,
    #NumSimulations, output, NumIsoformsToSimulate, filteringThreshold,
    #data_out)

    NegativeControlDrop("data/" * countsMatrix, GeneIsoformRelationships,
    NumSimulations, output * "_no_dropouts", NumIsoformsToSimulate, filteringThreshold, false,
    pFN, pFP, "none", "Weibull", data_out * "_no_dropouts")

    globalPredict("data/" * countsMatrix, GeneIsoformRelationships,
    NumSimulations, output * "_no_quant", NumIsoformsToSimulate, filteringThreshold,
    false, nothing,nothing, 0.0, 0.0, "none", "Weibull", data_out * "_no_quant")
end

################################################################################
# Make figures
################################################################################

function make_hist(path, titleArg)
    Arr = readdlm(path)

    if (titleArg == "3 Isoforms" &&
        path == "figures/negControl/H1_96_scnorm_counts.txt_3_isoforms_array.txt")
        #Untransformed plot
        p1 = histogram(Arr,
        title = titleArg,
        xlims = (0, 4.1),
        ylims = (0, 10),
        ylabel = "Probability Density",
        normalize=:pdf,
        bar_width = 0.1,
        bins = collect(0:0.01:4.1),
        axis = (font(100)), titlefontsize=100, yticks = 0:5:10,
        margin = 30mm)

    elseif titleArg == "Real" || titleArg == "alpha=0.45, beta=1.45"
        p1 = histogram(Arr,
        title = titleArg,
        xlims = (0, 4.1),
        ylims = (0, 10),
        normalize=:pdf,
        xlabel = "Mean No. Isoforms\nper Gene per Cell",
        bins = 100,
        axis = (font(100)), titlefontsize=100, yticks = 0:5:10,
        margin = 30mm)

    else
        #Untransformed plot
        p1 = histogram(Arr,
        title = titleArg,
        xlims = (0, 4.1),
        ylims = (0, 10),
        normalize=:pdf,
        bar_width = 0.1,
        bins = collect(0:0.01:4.1),
        axis = (font(100)), titlefontsize=100, yticks = 0:5:10,
        margin = 30mm)
    end
    return p1
end

H1_24_1 = make_hist("figures/negControl/H1_24_scnorm_counts.txt_1_isoforms_array.txt", "1 Isoform")
H1_24_2 = make_hist("figures/negControl/H1_24_scnorm_counts.txt_2_isoforms_array.txt", "2 Isoforms")
H1_24_3 = make_hist("figures/negControl/H1_24_scnorm_counts.txt_3_isoforms_array.txt", "3 Isoforms")
H1_24_4 = make_hist("figures/negControl/H1_24_scnorm_counts.txt_4_isoforms_array.txt", "4 Isoforms")
H1_24_real = make_hist("figures/negControl/H1_24_scnorm_counts.txt_no_quant_real_isoforms_array.txt", "Real")

H1_96_1 = make_hist("figures/negControl/H1_96_scnorm_counts.txt_1_isoforms_array.txt", "1 Isoform")
H1_96_2 = make_hist("figures/negControl/H1_96_scnorm_counts.txt_2_isoforms_array.txt", "2 Isoforms")
H1_96_3 = make_hist("figures/negControl/H1_96_scnorm_counts.txt_3_isoforms_array.txt", "3 Isoforms")
H1_96_4 = make_hist("figures/negControl/H1_96_scnorm_counts.txt_4_isoforms_array.txt", "4 Isoforms")
H1_96_real = make_hist("figures/negControl/H1_96_scnorm_counts.txt_no_quant_real_isoforms_array.txt", "Real")

pA = plot(H1_96_1,H1_24_1,
H1_96_2,H1_24_2,
H1_96_3,H1_24_3,
H1_96_4,H1_24_4,
H1_96_real,H1_24_real,
layout = (5,2), legend = false, size = (500*8, 600*8))

savefig(pA, "figures/negControl/H1_figs.png")

###################################
# no drop
###################################

function make_hist(path, titleArg)
    Arr = readdlm(path)

    if (titleArg == "3 Isoforms" &&
        path == "figures/negControl/H1_96_scnorm_counts.txt_3_isoforms_array.txt")
        #Untransformed plot
        p1 = histogram(Arr,
        title = titleArg,
        xlims = (0, 4.1),
        ylims = (0, 10),
        ylabel = "Probability Density",
        normalize=:pdf,
        bins =100,
        axis = (font(100)), titlefontsize=100, yticks = 0:5:10,
        margin = 30mm)

    elseif titleArg == "Real" || titleArg == "alpha=0.45, beta=1.45"
        p1 = histogram(Arr,
        title = titleArg,
        xlims = (0, 4.1),
        ylims = (0, 10),
        normalize=:pdf,
        xlabel = "Mean No. Isoforms\nper Gene per Cell",
        bins = 100,
        axis = (font(100)), titlefontsize=100, yticks = 0:5:10,
        margin = 30mm)

    else
        #Untransformed plot
        p1 = histogram(Arr,
        title = titleArg,
        xlims = (0, 4.1),
        ylims = (0, 10),
        normalize=:pdf,
        bins = 100,
        axis = (font(100)), titlefontsize=100, yticks = 0:5:10,
        margin = 30mm)
    end
    return p1
end

H1_24_1 = make_hist("figures/negControl/H1_24_scnorm_counts.txt_no_dropouts_1_isoforms_array.txt", "1 Isoform")
H1_24_2 = make_hist("figures/negControl/H1_24_scnorm_counts.txt_no_dropouts_2_isoforms_array.txt", "2 Isoforms")
H1_24_3 = make_hist("figures/negControl/H1_24_scnorm_counts.txt_no_dropouts_3_isoforms_array.txt", "3 Isoforms")
H1_24_4 = make_hist("figures/negControl/H1_24_scnorm_counts.txt_no_dropouts_4_isoforms_array.txt", "4 Isoforms")
H1_24_real = make_hist("figures/negControl/H1_24_scnorm_counts.txt_no_quant_real_isoforms_array.txt", "Real")

H1_96_1 = make_hist("figures/negControl/H1_96_scnorm_counts.txt_no_dropouts_1_isoforms_array.txt", "1 Isoform")
H1_96_2 = make_hist("figures/negControl/H1_96_scnorm_counts.txt_no_dropouts_2_isoforms_array.txt", "2 Isoforms")
H1_96_3 = make_hist("figures/negControl/H1_96_scnorm_counts.txt_no_dropouts_3_isoforms_array.txt", "3 Isoforms")
H1_96_4 = make_hist("figures/negControl/H1_96_scnorm_counts.txt_no_dropouts_4_isoforms_array.txt", "4 Isoforms")
H1_96_real = make_hist("figures/negControl/H1_96_scnorm_counts.txt_no_quant_real_isoforms_array.txt", "Real")

pA = plot(H1_96_1,H1_24_1,
H1_96_2,H1_24_2,
H1_96_3,H1_24_3,
H1_96_4,H1_24_4,
H1_96_real,H1_24_real,
layout = (5,2), legend = false, size = (500*8, 600*8))

savefig(pA, "figures/negControl/H1_figs_no_dropouts.png")

###################################
# no quant
###################################

H1_24_1 = make_hist("figures/negControl/H1_24_scnorm_counts.txt_no_quant_1_isoforms_array.txt", "1 Isoform")
H1_24_2 = make_hist("figures/negControl/H1_24_scnorm_counts.txt_no_quant_2_isoforms_array.txt", "2 Isoforms")
H1_24_3 = make_hist("figures/negControl/H1_24_scnorm_counts.txt_no_quant_3_isoforms_array.txt", "3 Isoforms")
H1_24_4 = make_hist("figures/negControl/H1_24_scnorm_counts.txt_no_quant_4_isoforms_array.txt", "4 Isoforms")
H1_24_real = make_hist("figures/negControl/H1_24_scnorm_counts.txt_no_quant_real_isoforms_array.txt", "Real")

H1_96_1 = make_hist("figures/negControl/H1_96_scnorm_counts.txt_no_quant_1_isoforms_array.txt", "1 Isoform")
H1_96_2 = make_hist("figures/negControl/H1_96_scnorm_counts.txt_no_quant_2_isoforms_array.txt", "2 Isoforms")
H1_96_3 = make_hist("figures/negControl/H1_96_scnorm_counts.txt_no_quant_3_isoforms_array.txt", "3 Isoforms")
H1_96_4 = make_hist("figures/negControl/H1_96_scnorm_counts.txt_no_quant_4_isoforms_array.txt", "4 Isoforms")
H1_96_real = make_hist("figures/negControl/H1_96_scnorm_counts.txt_no_quant_real_isoforms_array.txt", "Real")

pA = plot(H1_96_1,H1_24_1,
H1_96_2,H1_24_2,
H1_96_3,H1_24_3,
H1_96_4,H1_24_4,
H1_96_real,H1_24_real,
layout = (5,2), legend = false, size = (500*8, 600*8))

savefig(pA, "figures/negControl/H1_figs_no_quant.png")

###################################
# H9
####################################

H9_24_1 = make_hist("figures/negControl/H9_24_scnorm_counts.txt_1_isoforms_array.txt", "1 Isoform")
H9_24_2 = make_hist("figures/negControl/H9_24_scnorm_counts.txt_2_isoforms_array.txt", "2 Isoforms")
H9_24_3 = make_hist("figures/negControl/H9_24_scnorm_counts.txt_3_isoforms_array.txt", "3 Isoforms")
H9_24_4 = make_hist("figures/negControl/H9_24_scnorm_counts.txt_4_isoforms_array.txt", "4 Isoforms")
H9_24_real = make_hist("figures/figure2_data/H9_24_scnorm_counts.txt_real_isoforms_array.txt", "Real")

H9_96_1 = make_hist("figures/negControl/H9_96_scnorm_counts.txt_1_isoforms_array.txt", "1 Isoform")
H9_96_2 = make_hist("figures/negControl/H9_96_scnorm_counts.txt_2_isoforms_array.txt", "2 Isoforms")
H9_96_3 = make_hist("figures/negControl/H9_96_scnorm_counts.txt_3_isoforms_array.txt", "3 Isoforms")
H9_96_4 = make_hist("figures/negControl/H9_96_scnorm_counts.txt_4_isoforms_array.txt", "4 Isoforms")
H9_96_real = make_hist("figures/figure2_data/H9_96_scnorm_counts.txt_real_isoforms_array.txt", "Real")

pB = plot(H9_96_1,H9_24_1,
H9_96_2,H9_24_2,
H9_96_3,H9_24_3,
H9_96_4,H9_24_4,
H9_96_real,H9_24_real,
layout = (5,2), legend = false, size = (500*8, 600*8))

savefig(pB, "figures/negControl/H9_figs.png")
