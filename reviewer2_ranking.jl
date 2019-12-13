include("iPredict.jl/src/iPredict.jl")
using .iPredict
using DataFrames
using Plots
using JLD2
using DelimitedFiles
using StatsPlots

#Initialise variables
GeneIsoformRelationships = "data/gencode_human_gene_isoform_relationships.txt"
NumSimulations = 100
NumIsoformsToSimulate = 4
filteringThreshold = 4
alpha = nothing
beta = nothing
countsMatrix = "H1_24_scnorm_counts.txt"
error = collect(0:0.1:0.5)
pFP = 0.01
pFN = 0.04


for i in 1:length(error) # i = pFN
    for j in 1:length(error) # j = pFP

        #output directories
        output = "figures/reviewer2/" * string(error[i]) * "_" * string(error[j]) * "_" * countsMatrix * "fig"
        data_out = "figures/reviewer2/" * string(error[i]) * "_" * string(error[j]) * "_" * countsMatrix

        #simulate
        rankingPredict("data/" * countsMatrix, GeneIsoformRelationships,
        NumSimulations, output, NumIsoformsToSimulate, filteringThreshold,
       false, nothing,nothing, pFN, pFP, "none", "Weibull", data_out)


    end

end
    # globalPredict("data/" * countsMatrix, GeneIsoformRelationships,
    # NumSimulations, "figures/figure3/random_" * countsMatrix , NumIsoformsToSimulate, filteringThreshold,
    # false, alpha, beta, 0.04, 0.01, "random", "Weibull", "none")
    # globalPredict("data/" * countsMatrix, GeneIsoformRelationships,
    # NumSimulations, "figures/figure3/IsoformDependence_" * countsMatrix , NumIsoformsToSimulate, filteringThreshold,
    # false, alpha, beta, 0.04, 0.01, "IsoformDependence", "Weibull", "none")

#FN
Arr0 = readdlm("figures/reviewer2/0.0_0.0_H1_24_scnorm_counts.txt_spearmans_array.txt")
Arr1 = readdlm("figures/reviewer2/0.1_0.0_H1_24_scnorm_counts.txt_spearmans_array.txt")
Arr2 = readdlm("figures/reviewer2/0.2_0.0_H1_24_scnorm_counts.txt_spearmans_array.txt")
Arr3 = readdlm("figures/reviewer2/0.3_0.0_H1_24_scnorm_counts.txt_spearmans_array.txt")
Arr4 = readdlm("figures/reviewer2/0.4_0.0_H1_24_scnorm_counts.txt_spearmans_array.txt")
Arr5 = readdlm("figures/reviewer2/0.5_0.0_H1_24_scnorm_counts.txt_spearmans_array.txt")

test = ["0.0", "0.1", "0.2", "0.3", "0.4", "0.5"]
p1 = boxplot(test, [Arr0,Arr1,Arr2,Arr3,Arr4,Arr5], legend = false, xlabel = "pFN", ylabel = "Spearman's Rho")

savefig(p1, "figures/reviewer2/pFN_boxplots.png")

#FP
Arr0 = readdlm("figures/reviewer2/0.0_0.0_H1_24_scnorm_counts.txt_spearmans_array.txt")
Arr1 = readdlm("figures/reviewer2/0.0_0.1_H1_24_scnorm_counts.txt_spearmans_array.txt")
Arr2 = readdlm("figures/reviewer2/0.0_0.2_H1_24_scnorm_counts.txt_spearmans_array.txt")
Arr3 = readdlm("figures/reviewer2/0.0_0.3_H1_24_scnorm_counts.txt_spearmans_array.txt")
Arr4 = readdlm("figures/reviewer2/0.0_0.4_H1_24_scnorm_counts.txt_spearmans_array.txt")
Arr5 = readdlm("figures/reviewer2/0.0_0.5_H1_24_scnorm_counts.txt_spearmans_array.txt")

test = ["0.0", "0.1", "0.2", "0.3", "0.4", "0.5"]
p2 = boxplot(test, [Arr0,Arr1,Arr2,Arr3,Arr4,Arr5], legend = false, xlabel = "pFP", ylabel = "Spearman's Rho")

savefig(p2, "figures/reviewer2/pFP_boxplots.png")

##############################################################

rankingOriginal("data/" * countsMatrix, GeneIsoformRelationships,
NumSimulations, "figures/reviewer2/original", NumIsoformsToSimulate, filteringThreshold,
false, nothing,nothing, pFN, pFP, "none", "Weibull", "figures/reviewer2/original")
