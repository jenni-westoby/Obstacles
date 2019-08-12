include("iPredict.jl/src/iPredict.jl")
using .iPredict
using DataFrames
using Plots
using DelimitedFiles
using Plots.Measures


function make_hist(path, titleArg)
    Arr = readdlm(path)

    if titleArg == "3 Isoforms" && path == "figures/figure5_data/Weibull_H1_96_scnorm_counts.txt_3_isoforms_array.txt"
        #Untransformed plot
        p1 = histogram(Arr,
        title = titleArg,
        xlims = (0, 4),
        ylims = (0, 5),
        ylabel = "Normalised Frequency",
        normalize=:pdf,
        bins = 100,
        axis = (font(150)), titlefontsize=150, yticks = 0:5:5,
        left_margin = 75mm,
        right_margin = 75mm)

    elseif titleArg == "Real"
        p1 = histogram(Arr,
        title = titleArg,
        xlims = (0, 4),
        ylims = (0, 5),
        normalize=:pdf,
        xlabel = "Mean Isoforms\nper Gene per Cell",
        bins = 100,
        axis = (font(150)), titlefontsize=150, yticks = 0:5:5,
        left_margin = 75mm,
        right_margin = 75mm)

    else
        #Untransformed plot
        p1 = histogram(Arr,
        title = titleArg,
        xlims = (0, 4),
        ylims = (0, 5),
        normalize=:pdf,
        bins = 100,
        axis = (font(150)), titlefontsize=150, yticks = 0:5:5,
        left_margin = 75mm,
        right_margin = 75mm)
    end
    return p1
end

#for countsMatrix in scnorm_counts
# "Weibull", "Random", "UniformObserved", "CellVariable"
Weibull_1 = make_hist("figures/figure5_data/Weibull_H1_96_scnorm_counts.txt_1_isoforms_array.txt", "1 Isoform")
Weibull_2 = make_hist("figures/figure5_data/Weibull_H1_96_scnorm_counts.txt_2_isoforms_array.txt", "2 Isoforms")
Weibull_3 = make_hist("figures/figure5_data/Weibull_H1_96_scnorm_counts.txt_3_isoforms_array.txt", "3 Isoforms")
Weibull_4 = make_hist("figures/figure5_data/Weibull_H1_96_scnorm_counts.txt_4_isoforms_array.txt", "4 Isoforms")
Weibull_real = make_hist("figures/figure5_data/Weibull_H1_96_scnorm_counts.txt_real_isoforms_array.txt", "Real")

Random_1 = make_hist("figures/figure5_data/Random_H1_96_scnorm_counts.txt_1_isoforms_array.txt", "1 Isoform")
Random_2 = make_hist("figures/figure5_data/Random_H1_96_scnorm_counts.txt_2_isoforms_array.txt", "2 Isoforms")
Random_3 = make_hist("figures/figure5_data/Random_H1_96_scnorm_counts.txt_3_isoforms_array.txt", "3 Isoforms")
Random_4 = make_hist("figures/figure5_data/Random_H1_96_scnorm_counts.txt_4_isoforms_array.txt", "4 Isoforms")
Random_real = make_hist("figures/figure5_data/Random_H1_96_scnorm_counts.txt_real_isoforms_array.txt", "Real")

UniformObserved_1 = make_hist("figures/figure5_data/UniformObserved_H1_96_scnorm_counts.txt_1_isoforms_array.txt", "1 Isoform")
UniformObserved_2 = make_hist("figures/figure5_data/UniformObserved_H1_96_scnorm_counts.txt_2_isoforms_array.txt", "2 Isoforms")
UniformObserved_3 = make_hist("figures/figure5_data/UniformObserved_H1_96_scnorm_counts.txt_3_isoforms_array.txt", "3 Isoforms")
UniformObserved_4 = make_hist("figures/figure5_data/UniformObserved_H1_96_scnorm_counts.txt_4_isoforms_array.txt", "4 Isoforms")
UniformObserved_real = make_hist("figures/figure5_data/UniformObserved_H1_96_scnorm_counts.txt_real_isoforms_array.txt", "Real")

CellVariable_1 = make_hist("figures/figure5_data/CellVariable_H1_96_scnorm_counts.txt_1_isoforms_array.txt", "1 Isoform")
CellVariable_2 = make_hist("figures/figure5_data/CellVariable_H1_96_scnorm_counts.txt_2_isoforms_array.txt", "2 Isoforms")
CellVariable_3 = make_hist("figures/figure5_data/CellVariable_H1_96_scnorm_counts.txt_3_isoforms_array.txt", "3 Isoforms")
CellVariable_4 = make_hist("figures/figure5_data/CellVariable_H1_96_scnorm_counts.txt_4_isoforms_array.txt", "4 Isoforms")
CellVariable_real = make_hist("figures/figure5_data/CellVariable_H1_96_scnorm_counts.txt_real_isoforms_array.txt", "Real")

p3 = plot( Weibull_1, Random_1, UniformObserved_1, CellVariable_1,
       Weibull_2, Random_2,UniformObserved_2, CellVariable_2,
       Weibull_3, Random_3,UniformObserved_3, CellVariable_3,
       Weibull_4, Random_4,UniformObserved_4, CellVariable_4,
       Weibull_real, Random_real,UniformObserved_real, CellVariable_real,
       layout = (5,4), legend = false, size = (750 *8, 750*8), axis = (font(90)), titlefontsize=90,
       margin = 45mm)


savefig(p3, "figures/supplementary_figs/figure2.png")
