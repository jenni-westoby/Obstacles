include("iPredict.jl/src/iPredict.jl")
using .iPredict
using DataFrames
using Plots, Plots.Measures
using DelimitedFiles

#Initialise variables
GeneIsoformRelationships = "data/unfiltered_gencode_M20_gene_iso_relationships.txt"
NumSimulations = 100
NumIsoformsToSimulate = 4
filteringThreshold = 4
alpha = nothing
beta = nothing
scnorm_counts = ["BLUEPRINT_Kallisto_Counts_B.txt",
"E-MTAB-2600-standard_2i.txt"]
pFP = 0.01
pFN = 0.04
isoChoiceModels = "Weibull"


#Analyse BLUEPRINT
out = globalPredict("data/clean_BLUEPRINT_counts.txt",
GeneIsoformRelationships, NumSimulations,
"figures/figure1/BLUEPRINT_Bs.txt", NumIsoformsToSimulate,
filteringThreshold, false, alpha, beta, pFN, pFP, "none", isoChoiceModels,
"figures/figure1_data/BLUEPRINT")

#Analyse Kolod et al ES cells
out = globalPredict("data/E-MTAB-2600-standard_2i.txt",
GeneIsoformRelationships, NumSimulations,
"figures/figure1/E-MTAB-2600-standard_2i.txt", NumIsoformsToSimulate,
filteringThreshold, false, alpha, beta, pFN, pFP, "none", isoChoiceModels,
"figures/figure1_data/E-MTAB-2600")

function make_hist(path, titleArg)
    Arr = readdlm(path)

    if titleArg == "3 Isoforms" && path == "figures/figure1_data/BLUEPRINT_3_isoforms_array.txt"
        #Untransformed plot
        p1 = histogram(Arr,
        title = titleArg,
        xlims = (0, 4),
        ylims = (0, 10),
        ylabel = "Normalised Frequency",
        normalize=:pdf,
        bins = 100,
        axis = (font(80)), titlefontsize=80, yticks = 0:5:10,
        margin = 25mm)

    elseif titleArg == "Real"
        p1 = histogram(Arr,
        title = titleArg,
        xlims = (0, 4),
        ylims = (0, 10),
        normalize=:pdf,
        xlabel = "Mean Isoforms per Gene per Cell",
        bins = 100,
        axis = (font(80)), titlefontsize=80, yticks = 0:5:10,
        margin = 25mm)

    else
        #Untransformed plot
        p1 = histogram(Arr,
        title = titleArg,
        xlims = (0, 4),
        ylims = (0, 10),
        normalize=:pdf,
        bins = 100,
        axis = (font(80)), titlefontsize=80, yticks = 0:5:10,
        margin = 25mm)
    end
    return p1
end

BLUE_1 = make_hist("figures/figure1_data/BLUEPRINT_1_isoforms_array.txt", "1 Isoform")
BLUE_2 = make_hist("figures/figure1_data/BLUEPRINT_2_isoforms_array.txt", "2 Isoforms")
BLUE_3 = make_hist("figures/figure1_data/BLUEPRINT_3_isoforms_array.txt", "3 Isoforms")
BLUE_4 = make_hist("figures/figure1_data/BLUEPRINT_4_isoforms_array.txt", "4 Isoforms")
BLUE_real = make_hist("figures/figure1_data/BLUEPRINT_real_isoforms_array.txt", "Real")

kolod_1 = make_hist("figures/figure1_data/E-MTAB-2600_1_isoforms_array.txt", "1 Isoform")
kolod_2 = make_hist("figures/figure1_data/E-MTAB-2600_2_isoforms_array.txt", "2 Isoforms")
kolod_3 = make_hist("figures/figure1_data/E-MTAB-2600_3_isoforms_array.txt", "3 Isoforms")
kolod_4 = make_hist("figures/figure1_data/E-MTAB-2600_4_isoforms_array.txt", "4 Isoforms")
kolod_real = make_hist("figures/figure1_data/E-MTAB-2600_real_isoforms_array.txt", "Real")

out_plot = plot(BLUE_1,kolod_1, BLUE_2,kolod_2, BLUE_3, kolod_3,BLUE_4,kolod_4, BLUE_real,
  kolod_real, layout = (5, 2), legend = false, size = (500*8, 600*8))

savefig(out_plot, "figures/figure1/figure1.png")
