include("iPredict.jl/src/iPredict.jl")
using .iPredict
using DataFrames
using Plots, Plots.Measures
using HDF5
using DelimitedFiles
using Distributions
using StatsPlots

#Initialise variables
GeneIsoformRelationships = "data/gencode_human_gene_isoform_relationships.txt"
NumSimulations = 100
NumIsoformsToSimulate = 1
filteringThreshold = 4
alpha = nothing
beta = nothing
scnorm_counts = ["H1_24_scnorm_counts.txt",
"H1_96_scnorm_counts.txt",
"H9_24_scnorm_counts.txt",
"H9_96_scnorm_counts.txt"]
pFP = 0.01
pFN = 0.04
isoChoiceModels = ["Bernoulli", "Normal","Uniform"]

for countsMatrix in scnorm_counts

    output = "figures/figure6/" * countsMatrix

    for model in 1:length(isoChoiceModels)

        output = "figures/figure6/" * isoChoiceModels[model] * "_" * countsMatrix

        out = globalPredict("data/" * countsMatrix, GeneIsoformRelationships,
        NumSimulations, output, NumIsoformsToSimulate, filteringThreshold,
        false, alpha, beta, pFN, pFP, "none", isoChoiceModels[model],
        "figures/figure6_data/" * isoChoiceModels[model] * "_" * countsMatrix)
    end
end

###########################################################
# Make figures
function make_hist(path, titleArg)
    Arr = readdlm(path)

    if titleArg == "3 Isoforms" && path == "figures/figure5_data/Weibull_H1_24_scnorm_counts.txt_3_isoforms_array.txt"
        #Untransformed plot
        p1 = histogram(Arr,
        title = titleArg,
        xlims = (0, 4),
        ylims = (0, 5),
        ylabel = "Probability Density",
        normalize=:pdf,
        bins = 100)

    elseif titleArg == "Real"
        p1 = histogram(Arr,
        title = titleArg,
        xlims = (0, 4),
        ylims = (0, 5),
        normalize=:pdf,
        xlabel = "Mean No. Isoforms\nper Gene per Cell",
        bins = 100)

    else
        #Untransformed plot
        p1 = histogram(Arr,
        title = titleArg,
        xlims = (0, 4),
        ylims = (0, 5),
        normalize=:pdf,
        bins = 100,
        xlabel = "Mean No. Isoforms\nper Gene per Cell")
    end
    return p1
end

function make_drop_plots(path, output)

    test = load(path, "data")

    arr1 = []
    arr2 = []
    arr3 = []
    arr4 = []

    for (key, value) in test

        new_val = sort(value)
        push!(arr1, new_val[1])
        push!(arr2, new_val[2])
        push!(arr3, new_val[3])
        push!(arr4, new_val[4])
    end

    p1 = histogram(arr1, normalize = :pdf, title = "Rank 1", xlim = (0,1), ylim = (0,10) , ylabel="Probability Density" )
    p2 = histogram(arr2, normalize = :pdf, title = "Rank 2", xlim = (0,1), ylim = (0,10) )
    p3 = histogram(arr3, normalize = :pdf, title = "Rank 3", xlim = (0,1), ylim = (0,10) , ylabel="Probability Density", xlabel = "p(Dropout)")
    p4 = histogram(arr4, normalize = :pdf, title = "Rank 4", xlim = (0,1), ylim = (0,10), xlabel = "p(Dropout)" )

    pSave = plot(p1,p2,p3,p4, layout = (2,2), legend = :false,
    size = (500 * 8,500 *8), axis = (font(120)), titlefontsize=120, margin = 40mm)
    savefig(pSave, output)
    return pSave
end

function make_expr_plots(path, output)

    test = load(path, "data")

    arr1 = []
    arr2 = []
    arr3 = []
    arr4 = []

    for (key, value) in test

        new_val = sort(value, rev = true)
        push!(arr1, new_val[1])
        push!(arr2, new_val[2])
        push!(arr3, new_val[3])
        push!(arr4, new_val[4])
    end

    p1 = histogram(log.(arr1 .+ 1), normalize = :pdf, title = "Rank 1", xlim = (0,10), ylim = (0,2), ylabel="Probability Density" )
    p2 = histogram(log.(arr2 .+ 1), normalize = :pdf, title = "Rank 2", xlim = (0,10), ylim = (0,2) )
    p3 = histogram(log.(arr3 .+ 1), normalize = :pdf, title = "Rank 3", xlim = (0,10), ylim = (0,2), xlabel="log(counts + 1)", ylabel="Probability Density" )
    p4 = histogram(log.(arr4 .+ 1), normalize = :pdf, title = "Rank 4", xlim = (0,10), ylim = (0,2), xlabel="log(counts + 1)")

    pSave = plot(p1,p2,p3,p4, layout = (2,2), legend = :false,
    size = (500 * 8,500 *8), axis = (font(120)), titlefontsize=120, margin = 40mm)
    savefig(pSave, output)
    return pSave
end


normal_hist = make_hist("figures/figure6_data/Normal_H1_96_scnorm_counts.txt_1_isoforms_array.txt", "")
normal_hist = plot(normal_hist, ylabel = "Probability Density")
bernoulli_hist = make_hist("figures/figure6_data/Bernoulli_H1_96_scnorm_counts.txt_1_isoforms_array.txt", "")
uniform_hist = make_hist("figures/figure6_data/Uniform_H1_96_scnorm_counts.txt_1_isoforms_array.txt", "")

bernoulli_dist = bar(Bernoulli(0.25), bar_width = 0.1, ylim =(0,1), legend = false, title = "Bernoulli", xticks = 0:0.5:1, xlabel = "p(choice)")
normal_dist = plot(TruncatedNormal(0.25, 0.06, 0.0, 1.0), legend = false, title = "Normal", linestyle=:solid, linewidth=2*8, xticks = 0:0.5:1, ylabel = "Probability Density", xlabel = "p(choice)")
uniform_dist = histogram([ 0.25], ylim = (0,1), legend = false, title = "p=0.25", xlim = (0,1), bins = 0:0.1:1, xticks = 0:0.5:1, xlabel = "p(choice)")

p_A = plot(normal_dist, bernoulli_dist, uniform_dist,
normal_hist, bernoulli_hist, uniform_hist,
layout = (2,3), legend = false,
size = (1000 * 8,500 *8), axis = (font(120)), titlefontsize=120, margin = 80mm)

savefig(p_A, "figures/figure6/different_dists.png")

expr_H1_24 = make_expr_plots("figures/pdrop_scnorm_H1_24.pngexpr_dict.jld", "figures/figure6/H1_24_expr.png")
expr_H1_96 = make_expr_plots("figures/pdrop_scnorm_H1_96.pngexpr_dict.jld", "figures/figure6/H1_96_expr.png")
expr_H9_24 = make_expr_plots("figures/pdrop_scnorm_H9_24.pngexpr_dict.jld", "figures/figure6/H9_24_expr.png")
expr_H9_96 = make_expr_plots("figures/pdrop_scnorm_H9_96.pngexpr_dict.jld", "figures/figure6/H9_96_expr.png")

drop_H1_24 = make_drop_plots("figures/pdrop_scnorm_H1_24.png_dict.jld", "figures/figure6/H1_24_drop.png")
drop_H1_96 = make_drop_plots("figures/pdrop_scnorm_H1_96.png_dict.jld", "figures/figure6/H1_96_drop.png")
drop_H9_24 = make_drop_plots("figures/pdrop_scnorm_H9_24.png_dict.jld", "figures/figure6/H9_24_drop.png")
drop_H9_96 = make_drop_plots("figures/pdrop_scnorm_H9_96.png_dict.jld", "figures/figure6/H9_96_drop.png")


#make a supplementary fig
H1_24_normal_hist = make_hist("figures/figure6_data/Normal_H1_24_scnorm_counts.txt_1_isoforms_array.txt", "")
H1_24_normal_hist = plot(H1_24_normal_hist, ylabel = "Probability Density")
H1_24_bernoulli_hist = make_hist("figures/figure6_data/Bernoulli_H1_24_scnorm_counts.txt_1_isoforms_array.txt", "")
H1_24_uniform_hist = make_hist("figures/figure6_data/Uniform_H1_24_scnorm_counts.txt_1_isoforms_array.txt", "")

H9_24_normal_hist = make_hist("figures/figure6_data/Normal_H9_24_scnorm_counts.txt_1_isoforms_array.txt", "")
H9_24_normal_hist = plot(H9_24_normal_hist, ylabel = "Probability Density")
H9_24_bernoulli_hist = make_hist("figures/figure6_data/Bernoulli_H9_24_scnorm_counts.txt_1_isoforms_array.txt", "")
H9_24_uniform_hist = make_hist("figures/figure6_data/Uniform_H9_24_scnorm_counts.txt_1_isoforms_array.txt", "")

H9_96_normal_hist = make_hist("figures/figure6_data/Normal_H9_96_scnorm_counts.txt_1_isoforms_array.txt", "")
H9_96_normal_hist = plot(H9_96_normal_hist, ylabel = "Probability Density")
H9_96_bernoulli_hist = make_hist("figures/figure6_data/Bernoulli_H9_96_scnorm_counts.txt_1_isoforms_array.txt", "")
H9_96_uniform_hist = make_hist("figures/figure6_data/Uniform_H9_96_scnorm_counts.txt_1_isoforms_array.txt", "")

H1_96_normal_hist = make_hist("figures/figure6_data/Normal_H1_96_scnorm_counts.txt_1_isoforms_array.txt", "")
H1_96_normal_hist = plot(H1_96_normal_hist, ylabel = "Probability Density")
H1_96_bernoulli_hist = make_hist("figures/figure6_data/Bernoulli_H1_96_scnorm_counts.txt_1_isoforms_array.txt", "")
H1_96_uniform_hist = make_hist("figures/figure6_data/Uniform_H1_96_scnorm_counts.txt_1_isoforms_array.txt", "")


p_supp = plot(normal_dist, bernoulli_dist, uniform_dist,
H1_96_normal_hist, H1_96_bernoulli_hist, H1_96_uniform_hist,
H1_24_normal_hist, H1_24_bernoulli_hist, H1_24_uniform_hist,
H9_96_normal_hist, H9_96_bernoulli_hist, H9_96_uniform_hist,
H9_24_normal_hist, H9_24_bernoulli_hist, H9_24_uniform_hist,
layout = (5,3), legend = false,
size = (1000 * 8,1200 *8), axis = (font(120)), titlefontsize=120, margin = 80mm)

savefig(p_supp, "figures/supplementary_figs/different_dists.png")

###############################################################
using JLD2
using DelimitedFiles
using HypothesisTests
using DataFrames

df = DataFrame(IsoformsSimulated=Int64[], Comparison = String[], pValue = Float64[])

normal_H1_96 = readdlm("figures/figure6_data/Normal_H1_96_scnorm_counts.txt_1_isoforms_array.txt")
bernoulli_H1_96 = readdlm("figures/figure6_data/Bernoulli_H1_96_scnorm_counts.txt_1_isoforms_array.txt")
uniform_hist_H1_96 = readdlm("figures/figure6_data/Uniform_H1_96_scnorm_counts.txt_1_isoforms_array.txt")
append!(df, DataFrame(IsoformsSimulated = [1], Comparison = ["H1_96"], pValue=[pvalue(KSampleADTest(vec(normal_H1_96), vec(bernoulli_H1_96), vec(uniform_hist_H1_96)))]))

normal_H1_24 = readdlm("figures/figure6_data/Normal_H1_24_scnorm_counts.txt_1_isoforms_array.txt")
bernoulli_H1_24 = readdlm("figures/figure6_data/Bernoulli_H1_24_scnorm_counts.txt_1_isoforms_array.txt")
uniform_hist_H1_24 = readdlm("figures/figure6_data/Uniform_H1_24_scnorm_counts.txt_1_isoforms_array.txt")
append!(df, DataFrame(IsoformsSimulated = [1], Comparison = ["H1_24"], pValue=[pvalue(KSampleADTest(vec(normal_H1_24), vec(bernoulli_H1_24), vec(uniform_hist_H1_24)))]))

normal_H9_96 = readdlm("figures/figure6_data/Normal_H9_96_scnorm_counts.txt_1_isoforms_array.txt")
bernoulli_H9_96 = readdlm("figures/figure6_data/Bernoulli_H9_96_scnorm_counts.txt_1_isoforms_array.txt")
uniform_hist_H9_96 = readdlm("figures/figure6_data/Uniform_H9_96_scnorm_counts.txt_1_isoforms_array.txt")
append!(df, DataFrame(IsoformsSimulated = [1], Comparison = ["H9_96"], pValue=[pvalue(KSampleADTest(vec(normal_H9_96), vec(bernoulli_H9_96), vec(uniform_hist_H9_96)))]))

normal_H9_24 = readdlm("figures/figure6_data/Normal_H9_24_scnorm_counts.txt_1_isoforms_array.txt")
bernoulli_H9_24 = readdlm("figures/figure6_data/Bernoulli_H9_24_scnorm_counts.txt_1_isoforms_array.txt")
uniform_hist_H9_24 = readdlm("figures/figure6_data/Uniform_H9_24_scnorm_counts.txt_1_isoforms_array.txt")
append!(df, DataFrame(IsoformsSimulated = [1], Comparison = ["H9_24"], pValue=[pvalue(KSampleADTest(vec(normal_H9_24), vec(bernoulli_H9_24), vec(uniform_hist_H9_24)))]))

jldopen("figures/figure6_data/Supfig_AD_test_results.jld2", "w") do file
    file["df"] = df
end

#############################################
# Overlap
#############################################

function make_hist(path, titleArg)
    Arr = readdlm(path)

    if titleArg == "3 Isoforms" && path == "figures/figure5_data/Weibull_H1_24_scnorm_counts.txt_3_isoforms_array.txt"
        #Untransformed plot
        p1 = histogram(Arr,
        title = titleArg,
        xlims = (0, 1),
        ylims = (0, 10),
        ylabel = "Probability Density",
        normalize=:pdf,
        bins = 50,
        color = palette(:default)[2],)

    elseif titleArg == "Real"
        p1 = histogram(Arr,
        title = titleArg,
        xlims = (0, 1),
        ylims = (0, 10),
        normalize=:pdf,
        xlabel = "Overlap With Ground Truth",
        bins = 50,
        color = palette(:default)[2],)

    else
        #Untransformed plot
        p1 = histogram(Arr,
        title = titleArg,
        xlims = (0, 1),
        ylims = (0, 10),
        normalize=:pdf,
        bins = 50,
        color = palette(:default)[2],
        xlabel = "Overlap With Ground Truth")
    end
    return p1
end

#make a supplementary fig
H1_24_normal_hist = make_hist("figures/figure6_data/Normal_H1_24_scnorm_counts.txt_1_isoforms_overlap_array.txt", "")
H1_24_normal_hist = plot(H1_24_normal_hist, ylabel = "Probability Density")
H1_24_bernoulli_hist = make_hist("figures/figure6_data/Bernoulli_H1_24_scnorm_counts.txt_1_isoforms_overlap_array.txt", "")
H1_24_uniform_hist = make_hist("figures/figure6_data/Uniform_H1_24_scnorm_counts.txt_1_isoforms_overlap_array.txt", "")

H9_24_normal_hist = make_hist("figures/figure6_data/Normal_H9_24_scnorm_counts.txt_1_isoforms_overlap_array.txt", "")
H9_24_normal_hist = plot(H9_24_normal_hist, ylabel = "Probability Density")
H9_24_bernoulli_hist = make_hist("figures/figure6_data/Bernoulli_H9_24_scnorm_counts.txt_1_isoforms_overlap_array.txt", "")
H9_24_uniform_hist = make_hist("figures/figure6_data/Uniform_H9_24_scnorm_counts.txt_1_isoforms_overlap_array.txt", "")

H9_96_normal_hist = make_hist("figures/figure6_data/Normal_H9_96_scnorm_counts.txt_1_isoforms_overlap_array.txt", "")
H9_96_normal_hist = plot(H9_96_normal_hist, ylabel = "Probability Density")
H9_96_bernoulli_hist = make_hist("figures/figure6_data/Bernoulli_H9_96_scnorm_counts.txt_1_isoforms_overlap_array.txt", "")
H9_96_uniform_hist = make_hist("figures/figure6_data/Uniform_H9_96_scnorm_counts.txt_1_isoforms_overlap_array.txt", "")

H1_96_normal_hist = make_hist("figures/figure6_data/Normal_H1_96_scnorm_counts.txt_1_isoforms_overlap_array.txt", "")
H1_96_normal_hist = plot(H1_96_normal_hist, ylabel = "Probability Density")
H1_96_bernoulli_hist = make_hist("figures/figure6_data/Bernoulli_H1_96_scnorm_counts.txt_1_isoforms_overlap_array.txt", "")
H1_96_uniform_hist = make_hist("figures/figure6_data/Uniform_H1_96_scnorm_counts.txt_1_isoforms_overlap_array.txt", "")


p_supp = plot(normal_dist, bernoulli_dist, uniform_dist,
H1_96_normal_hist, H1_96_bernoulli_hist, H1_96_uniform_hist,
H1_24_normal_hist, H1_24_bernoulli_hist, H1_24_uniform_hist,
H9_96_normal_hist, H9_96_bernoulli_hist, H9_96_uniform_hist,
H9_24_normal_hist, H9_24_bernoulli_hist, H9_24_uniform_hist,
layout = (5,3), legend = false,
size = (1000 * 8,1200 *8), axis = (font(120)), titlefontsize=120, margin = 80mm)

savefig(p_supp, "figures/supplementary_figs/different_dists_overlap.png")
