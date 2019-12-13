include("iPredict.jl/src/iPredict.jl")
using .iPredict
using DataFrames
using Plots
using DelimitedFiles
using Plots.Measures


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
        bins = 100,
        axis = (font(150)), titlefontsize=150, yticks = 0:5:5,
        margin = 35mm,
        left_margin = 50mm,
        right_margin = 50mm)

    elseif titleArg == "Real"
        p1 = histogram(Arr,
        title = titleArg,
        xlims = (0, 4),
        ylims = (0, 5),
        normalize=:pdf,
        xlabel = "Mean No. Isoforms\nper Gene per Cell",
        bins = 100,
        axis = (font(150)), titlefontsize=150, yticks = 0:5:5,
        margin = 35mm,
        left_margin = 50mm,
        right_margin = 50mm)

    else
        #Untransformed plot
        p1 = histogram(Arr,
        title = titleArg,
        xlims = (0, 4),
        ylims = (0, 5),
        normalize=:pdf,
        bins = 100,
        axis = (font(150)), titlefontsize=150, yticks = 0:5:5,
        margin = 35mm,left_margin = 50mm,
        right_margin = 50mm)
    end
    #vline!([mean(Arr)], color = :black, linewidth = 8)
    return p1
end
#for countsMatrix in scnorm_counts
# "Weibull", "Random", "UniformObserved", "CellVariable"
Weibull_1 = make_hist("figures/figure5_data/Weibull_H9_96_scnorm_counts.txt_1_isoforms_array.txt", "1 Isoform")
Weibull_2 = make_hist("figures/figure5_data/Weibull_H9_96_scnorm_counts.txt_2_isoforms_array.txt", "2 Isoforms")
Weibull_3 = make_hist("figures/figure5_data/Weibull_H9_96_scnorm_counts.txt_3_isoforms_array.txt", "3 Isoforms")
Weibull_4 = make_hist("figures/figure5_data/Weibull_H9_96_scnorm_counts.txt_4_isoforms_array.txt", "4 Isoforms")
Weibull_real = make_hist("figures/figure5_data/Weibull_H9_96_scnorm_counts.txt_real_isoforms_array.txt", "Real")

Random_1 = make_hist("figures/figure5_data/Random_H9_96_scnorm_counts.txt_1_isoforms_array.txt", "1 Isoform")
Random_2 = make_hist("figures/figure5_data/Random_H9_96_scnorm_counts.txt_2_isoforms_array.txt", "2 Isoforms")
Random_3 = make_hist("figures/figure5_data/Random_H9_96_scnorm_counts.txt_3_isoforms_array.txt", "3 Isoforms")
Random_4 = make_hist("figures/figure5_data/Random_H9_96_scnorm_counts.txt_4_isoforms_array.txt", "4 Isoforms")
Random_real = make_hist("figures/figure5_data/Random_H9_96_scnorm_counts.txt_real_isoforms_array.txt", "Real")

UniformObserved_1 = make_hist("figures/figure5_data/UniformObserved_H9_96_scnorm_counts.txt_1_isoforms_array.txt", "1 Isoform")
UniformObserved_2 = make_hist("figures/figure5_data/UniformObserved_H9_96_scnorm_counts.txt_2_isoforms_array.txt", "2 Isoforms")
UniformObserved_3 = make_hist("figures/figure5_data/UniformObserved_H9_96_scnorm_counts.txt_3_isoforms_array.txt", "3 Isoforms")
UniformObserved_4 = make_hist("figures/figure5_data/UniformObserved_H9_96_scnorm_counts.txt_4_isoforms_array.txt", "4 Isoforms")
UniformObserved_real = make_hist("figures/figure5_data/UniformObserved_H9_96_scnorm_counts.txt_real_isoforms_array.txt", "Real")

CellVariable_1 = make_hist("figures/figure5_data/CellVariable_H9_96_scnorm_counts.txt_1_isoforms_array.txt", "1 Isoform")
CellVariable_2 = make_hist("figures/figure5_data/CellVariable_H9_96_scnorm_counts.txt_2_isoforms_array.txt", "2 Isoforms")
CellVariable_3 = make_hist("figures/figure5_data/CellVariable_H9_96_scnorm_counts.txt_3_isoforms_array.txt", "3 Isoforms")
CellVariable_4 = make_hist("figures/figure5_data/CellVariable_H9_96_scnorm_counts.txt_4_isoforms_array.txt", "4 Isoforms")
CellVariable_real = make_hist("figures/figure5_data/CellVariable_H9_96_scnorm_counts.txt_real_isoforms_array.txt", "Real")

p3 = plot( Weibull_1, Random_1, UniformObserved_1, CellVariable_1,
       Weibull_2, Random_2,UniformObserved_2, CellVariable_2,
       Weibull_3, Random_3,UniformObserved_3, CellVariable_3,
       Weibull_4, Random_4,UniformObserved_4, CellVariable_4,
       Weibull_real, Random_real,UniformObserved_real, CellVariable_real,
       layout = (5,4), legend = false, size = (750 *8, 750*8), axis = (font(90)), titlefontsize=90,
       margin = 45mm)


savefig(p3, "figures/supplementary_figs/figure4_sup.png")

###############################################################
using JLD2
using DelimitedFiles
using HypothesisTests
using DataFrames

Weibull_1 = readdlm("figures/figure5_data/Weibull_H9_96_scnorm_counts.txt_1_isoforms_array.txt")
Weibull_2 = readdlm("figures/figure5_data/Weibull_H9_96_scnorm_counts.txt_2_isoforms_array.txt")
Weibull_3 = readdlm("figures/figure5_data/Weibull_H9_96_scnorm_counts.txt_3_isoforms_array.txt")
Weibull_4 = readdlm("figures/figure5_data/Weibull_H9_96_scnorm_counts.txt_4_isoforms_array.txt")
Weibull_real = readdlm("figures/figure5_data/Weibull_H9_96_scnorm_counts.txt_real_isoforms_array.txt")

Random_1 = readdlm("figures/figure5_data/Random_H9_96_scnorm_counts.txt_1_isoforms_array.txt")
Random_2 = readdlm("figures/figure5_data/Random_H9_96_scnorm_counts.txt_2_isoforms_array.txt")
Random_3 = readdlm("figures/figure5_data/Random_H9_96_scnorm_counts.txt_3_isoforms_array.txt")
Random_4 = readdlm("figures/figure5_data/Random_H9_96_scnorm_counts.txt_4_isoforms_array.txt")
Random_real = readdlm("figures/figure5_data/Random_H9_96_scnorm_counts.txt_real_isoforms_array.txt")

UniformObserved_1 = readdlm("figures/figure5_data/UniformObserved_H9_96_scnorm_counts.txt_1_isoforms_array.txt")
UniformObserved_2 = readdlm("figures/figure5_data/UniformObserved_H9_96_scnorm_counts.txt_2_isoforms_array.txt")
UniformObserved_3 = readdlm("figures/figure5_data/UniformObserved_H9_96_scnorm_counts.txt_3_isoforms_array.txt")
UniformObserved_4 = readdlm("figures/figure5_data/UniformObserved_H9_96_scnorm_counts.txt_4_isoforms_array.txt")
UniformObserved_real = readdlm("figures/figure5_data/UniformObserved_H9_96_scnorm_counts.txt_real_isoforms_array.txt")

CellVariable_1 = readdlm("figures/figure5_data/CellVariable_H9_96_scnorm_counts.txt_1_isoforms_array.txt")
CellVariable_2 = readdlm("figures/figure5_data/CellVariable_H9_96_scnorm_counts.txt_2_isoforms_array.txt")
CellVariable_3 = readdlm("figures/figure5_data/CellVariable_H9_96_scnorm_counts.txt_3_isoforms_array.txt")
CellVariable_4 = readdlm("figures/figure5_data/CellVariable_H9_96_scnorm_counts.txt_4_isoforms_array.txt")
CellVariable_real = readdlm("figures/figure5_data/CellVariable_H9_96_scnorm_counts.txt_real_isoforms_array.txt")

df = DataFrame(IsoformsSimulated=Int64[], Comparison = String[], pValue = Float64[])

append!(df, DataFrame(IsoformsSimulated = [1], Comparison = ["All"], pValue=[pvalue(KSampleADTest(vec(Weibull_1), vec(Random_1), vec(UniformObserved_1), vec(CellVariable_1)))]))
append!(df, DataFrame(IsoformsSimulated = [2], Comparison = "All", pValue=[pvalue(KSampleADTest(vec(Weibull_2), vec(Random_2), vec(UniformObserved_2), vec(CellVariable_2)))]))
append!(df, DataFrame(IsoformsSimulated = [3], Comparison = "All", pValue=[pvalue(KSampleADTest(vec(Weibull_3), vec(Random_3), vec(UniformObserved_3), vec(CellVariable_3)))]))
append!(df, DataFrame(IsoformsSimulated = [4], Comparison = "All", pValue=[pvalue(KSampleADTest(vec(Weibull_4), vec(Random_4), vec(UniformObserved_4), vec(CellVariable_4)))]))

append!(df, DataFrame(IsoformsSimulated = [1], Comparison = "CellVarUnif", pValue=[pvalue(KSampleADTest(vec(UniformObserved_1), vec(CellVariable_1)))]))
append!(df, DataFrame(IsoformsSimulated = [2], Comparison = "CellVarUnif", pValue=[pvalue(KSampleADTest( vec(UniformObserved_2), vec(CellVariable_2)))]))
append!(df, DataFrame(IsoformsSimulated = [3], Comparison = "CellVarUnif", pValue=[pvalue(KSampleADTest(vec(UniformObserved_3), vec(CellVariable_3)))]))
append!(df, DataFrame(IsoformsSimulated = [4], Comparison = "CellVarUnif", pValue=[pvalue(KSampleADTest(vec(UniformObserved_4), vec(CellVariable_4)))]))

jldopen("figures/figure4_data/Supfig4_AD_test_results.jld2", "w") do file
    file["df"] = df
end
#####################
#overlap

function make_hist(path, titleArg)

     Arr = readdlm(path)
     p1 = histogram(Arr, title = titleArg, normalize=:pdf, bins = 100,
     ylims = (0,4), color = palette(:default)[2],
     xlims = (0,1),
     axis = (font(150)), titlefontsize=150, yticks = 0:2:4,
     xticks = 0:0.5:1,
     margin = 35mm,left_margin = 50mm,
     right_margin = 50mm)
     return p1

end

Weibull_1 = make_hist("figures/figure5_data/Weibull_H9_96_scnorm_counts.txt_1_isoforms_overlap_array.txt", "1 Isoform")
Weibull_2 = make_hist("figures/figure5_data/Weibull_H9_96_scnorm_counts.txt_2_isoforms_overlap_array.txt", "2 Isoforms")
Weibull_3 = make_hist("figures/figure5_data/Weibull_H9_96_scnorm_counts.txt_3_isoforms_overlap_array.txt", "3 Isoforms")
Weibull_4 = make_hist("figures/figure5_data/Weibull_H9_96_scnorm_counts.txt_4_isoforms_overlap_array.txt", "4 Isoforms")

Random_1 = make_hist("figures/figure5_data/Random_H9_96_scnorm_counts.txt_1_isoforms_overlap_array.txt", "1 Isoform")
Random_2 = make_hist("figures/figure5_data/Random_H9_96_scnorm_counts.txt_2_isoforms_overlap_array.txt", "2 Isoforms")
Random_3 = make_hist("figures/figure5_data/Random_H9_96_scnorm_counts.txt_3_isoforms_overlap_array.txt", "3 Isoforms")
Random_4 = make_hist("figures/figure5_data/Random_H9_96_scnorm_counts.txt_4_isoforms_overlap_array.txt", "4 Isoforms")

UniformObserved_1 = make_hist("figures/figure5_data/UniformObserved_H9_96_scnorm_counts.txt_1_isoforms_overlap_array.txt", "1 Isoform")
UniformObserved_2 = make_hist("figures/figure5_data/UniformObserved_H9_96_scnorm_counts.txt_2_isoforms_overlap_array.txt", "2 Isoforms")
UniformObserved_3 = make_hist("figures/figure5_data/UniformObserved_H9_96_scnorm_counts.txt_3_isoforms_overlap_array.txt", "3 Isoforms")
UniformObserved_4 = make_hist("figures/figure5_data/UniformObserved_H9_96_scnorm_counts.txt_4_isoforms_overlap_array.txt", "4 Isoforms")

CellVariable_1 = make_hist("figures/figure5_data/CellVariable_H9_96_scnorm_counts.txt_1_isoforms_overlap_array.txt", "1 Isoform")
CellVariable_2 = make_hist("figures/figure5_data/CellVariable_H9_96_scnorm_counts.txt_2_isoforms_overlap_array.txt", "2 Isoforms")
CellVariable_3 = make_hist("figures/figure5_data/CellVariable_H9_96_scnorm_counts.txt_3_isoforms_overlap_array.txt", "3 Isoforms")
CellVariable_4 = make_hist("figures/figure5_data/CellVariable_H9_96_scnorm_counts.txt_4_isoforms_overlap_array.txt", "4 Isoforms")

p3 = plot( Weibull_1, Random_1, UniformObserved_1, CellVariable_1,
       Weibull_2, Random_2,UniformObserved_2, CellVariable_2,
       Weibull_3, Random_3,UniformObserved_3, CellVariable_3,
       Weibull_4, Random_4,UniformObserved_4, CellVariable_4,
       layout = (4,4), legend = false, size = (750 *8, 750*8), axis = (font(90)), titlefontsize=90,
       margin = 45mm)


savefig(p3, "figures/supplementary_figs/figure4_sup_overlap.png")
