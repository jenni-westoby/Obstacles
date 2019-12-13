include("iPredict.jl/src/iPredict.jl")
using .iPredict
using DataFrames
using Plots, Plots.Measures
using JLD2
using DelimitedFiles
using StatsBase

#Initialise variables
GeneIsoformRelationships = "data/gencode_human_gene_isoform_relationships.txt"
NumSimulations = 100
NumIsoformsToSimulate = 4
filteringThreshold = 4
alpha_vals = [0.45, 0.7, 0.95, 1.2, 1.45]
beta_vals = reverse(alpha_vals)
scnorm_counts = ["H1_24_scnorm_counts.txt"]#,
"H1_96_scnorm_counts.txt",
"H9_24_scnorm_counts.txt",
"H9_96_scnorm_counts.txt"]
pFP = 0.01
pFN = 0.04
isoChoiceModels = ["Random", "UniformObserved", "CellVariable"]




for countsMatrix in scnorm_counts

    output = "figures/figure2/" * countsMatrix
    data_out = "figures/figure2_data/" * countsMatrix

    for model in isoChoiceModels

        for i in 1:length(alpha_vals)

            output = "figures/figure2/" * model * string(alpha_vals[i]) * "_" *
            string(beta_vals[i]) * "_" * countsMatrix
            data_out = "figures/figure2_data/" * model * string(alpha_vals[i]) * "_" *
            string(beta_vals[i]) * "_" * countsMatrix

            out = globalPredict("data/" * countsMatrix, GeneIsoformRelationships,
            NumSimulations, output, NumIsoformsToSimulate, filteringThreshold,
            false, alpha_vals[i], beta_vals[i], pFN, pFP, "none", model,
            data_out)
        end
    end
end

function make_hist(path, titleArg)
    Arr = readdlm(path)

    if (titleArg == "3 Isoforms" &&
        path == "figures/figure2_data/H1_96_scnorm_counts.txt_3_isoforms_array.txt")
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

v0 =  make_hist("figures/figure2_data/Random0.45_1.45_H1_24_scnorm_counts.txt_1_isoforms_array.txt", "alpha=0.45, beta=1.45")
v1 =  make_hist("figures/figure2_data/Random0.7_1.2_H1_24_scnorm_counts.txt_1_isoforms_array.txt", "alpha=0.7, beta=1.2")
v2 =  make_hist("figures/figure2_data/Random0.95_0.95_H1_24_scnorm_counts.txt_1_isoforms_array.txt", "alpha=0.95, beta=0.95")
v3 =  make_hist("figures/figure2_data/Random1.2_0.7_H1_24_scnorm_counts.txt_1_isoforms_array.txt", "alpha=1.2, beta=0.7")
v4 =  make_hist("figures/figure2_data/Random1.45_0.45_H1_24_scnorm_counts.txt_1_isoforms_array.txt", "alpha=1.45, beta=0.45")

pD = plot(v4, v3, v2, v1, v0,
       layout = (5,1), legend = false, size = (250*8, 600*8),
       axis = (font(180)), titlefontsize=160, margin = 30mm)

savefig(pD, "figures/figure2/Random_figure2D.png")

###############

v0 =  make_hist("figures/figure2_data/UniformObserved0.45_1.45_H1_24_scnorm_counts.txt_1_isoforms_array.txt", "alpha=0.45, beta=1.45")
v1 =  make_hist("figures/figure2_data/UniformObserved0.7_1.2_H1_24_scnorm_counts.txt_1_isoforms_array.txt", "alpha=0.7, beta=1.2")
v2 =  make_hist("figures/figure2_data/UniformObserved0.95_0.95_H1_24_scnorm_counts.txt_1_isoforms_array.txt", "alpha=0.95, beta=0.95")
v3 =  make_hist("figures/figure2_data/UniformObserved1.2_0.7_H1_24_scnorm_counts.txt_1_isoforms_array.txt", "alpha=1.2, beta=0.7")
v4 =  make_hist("figures/figure2_data/UniformObserved1.45_0.45_H1_24_scnorm_counts.txt_1_isoforms_array.txt", "alpha=1.45, beta=0.45")

pD = plot(v4, v3, v2, v1, v0,
       layout = (5,1), legend = false, size = (250*8, 600*8),
       axis = (font(180)), titlefontsize=160, margin = 30mm)

savefig(pD, "figures/figure2/UniformObserved_figure2D.png")

#################

v0 =  make_hist("figures/figure2_data/CellVariable0.45_1.45_H1_24_scnorm_counts.txt_1_isoforms_array.txt", "alpha=0.45, beta=1.45")
v1 =  make_hist("figures/figure2_data/CellVariable0.7_1.2_H1_24_scnorm_counts.txt_1_isoforms_array.txt", "alpha=0.7, beta=1.2")
v2 =  make_hist("figures/figure2_data/CellVariable0.95_0.95_H1_24_scnorm_counts.txt_1_isoforms_array.txt", "alpha=0.95, beta=0.95")
v3 =  make_hist("figures/figure2_data/CellVariable1.2_0.7_H1_24_scnorm_counts.txt_1_isoforms_array.txt", "alpha=1.2, beta=0.7")
v4 =  make_hist("figures/figure2_data/CellVariable1.45_0.45_H1_24_scnorm_counts.txt_1_isoforms_array.txt", "alpha=1.45, beta=0.45")

pD = plot(v4, v3, v2, v1, v0,
       layout = (5,1), legend = false, size = (250*8, 600*8),
       axis = (font(180)), titlefontsize=160, margin = 30mm)

savefig(pD, "figures/figure2/CellVariable_figure2D.png")

####################################################
# OVERLAP
####################################################

function make_hist(path, titleArg)
    Arr = readdlm(path)

    if (titleArg == "3 Isoforms" &&
        path == "figures/figure2_data/H1_96_scnorm_counts.txt_3_isoforms_array.txt")
        #Untransformed plot
        p1 = histogram(Arr,
        title = titleArg,
        color = palette(:default)[2],
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
        color = palette(:default)[2],
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
        color = palette(:default)[2],
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

v0 =  make_hist("figures/figure2_data/Random0.45_1.45_H1_24_scnorm_counts.txt_1_isoforms_overlap_array.txt", "alpha=0.45, beta=1.45")
v1 =  make_hist("figures/figure2_data/Random0.7_1.2_H1_24_scnorm_counts.txt_1_isoforms_overlap_array.txt", "alpha=0.7, beta=1.2")
v2 =  make_hist("figures/figure2_data/Random0.95_0.95_H1_24_scnorm_counts.txt_1_isoforms_overlap_array.txt", "alpha=0.95, beta=0.95")
v3 =  make_hist("figures/figure2_data/Random1.2_0.7_H1_24_scnorm_counts.txt_1_isoforms_overlap_array.txt", "alpha=1.2, beta=0.7")
v4 =  make_hist("figures/figure2_data/Random1.45_0.45_H1_24_scnorm_counts.txt_1_isoforms_overlap_array.txt", "alpha=1.45, beta=0.45")

pD = plot(v4, v3, v2, v1, v0,
       layout = (5,1), legend = false, size = (250*8, 600*8),
       axis = (font(180)), titlefontsize=160, margin = 30mm)

savefig(pD, "figures/figure2/Random_figure2D_overlap.png")

###############

v0 =  make_hist("figures/figure2_data/UniformObserved0.45_1.45_H1_24_scnorm_counts.txt_1_isoforms_overlap_array.txt", "alpha=0.45, beta=1.45")
v1 =  make_hist("figures/figure2_data/UniformObserved0.7_1.2_H1_24_scnorm_counts.txt_1_isoforms_overlap_array.txt", "alpha=0.7, beta=1.2")
v2 =  make_hist("figures/figure2_data/UniformObserved0.95_0.95_H1_24_scnorm_counts.txt_1_isoforms_overlap_array.txt", "alpha=0.95, beta=0.95")
v3 =  make_hist("figures/figure2_data/UniformObserved1.2_0.7_H1_24_scnorm_counts.txt_1_isoforms_overlap_array.txt", "alpha=1.2, beta=0.7")
v4 =  make_hist("figures/figure2_data/UniformObserved1.45_0.45_H1_24_scnorm_counts.txt_1_isoforms_overlap_array.txt", "alpha=1.45, beta=0.45")

pD = plot(v4, v3, v2, v1, v0,
       layout = (5,1), legend = false, size = (250*8, 600*8),
       axis = (font(180)), titlefontsize=160, margin = 30mm)

savefig(pD, "figures/figure2/UniformObserved_figure2D_overlap.png")

#################

v0 =  make_hist("figures/figure2_data/CellVariable0.45_1.45_H1_24_scnorm_counts.txt_1_isoforms_overlap_array.txt", "alpha=0.45, beta=1.45")
v1 =  make_hist("figures/figure2_data/CellVariable0.7_1.2_H1_24_scnorm_counts.txt_1_isoforms_overlap_array.txt", "alpha=0.7, beta=1.2")
v2 =  make_hist("figures/figure2_data/CellVariable0.95_0.95_H1_24_scnorm_counts.txt_1_isoforms_overlap_array.txt", "alpha=0.95, beta=0.95")
v3 =  make_hist("figures/figure2_data/CellVariable1.2_0.7_H1_24_scnorm_counts.txt_1_isoforms_overlap_array.txt", "alpha=1.2, beta=0.7")
v4 =  make_hist("figures/figure2_data/CellVariable1.45_0.45_H1_24_scnorm_counts.txt_1_isoforms_overlap_array.txt", "alpha=1.45, beta=0.45")

pD = plot(v4, v3, v2, v1, v0,
       layout = (5,1), legend = false, size = (250*8, 600*8),
       axis = (font(180)), titlefontsize=160, margin = 30mm)

savefig(pD, "figures/figure2/CellVariable_figure2D_overlap.png")
