include("iPredict.jl/src/MixtureModelling.jl")
using DelimitedFiles
using Distributions
using StatsPlots
using DataFrames
using Plots, Plots.Measures

scnorm_counts = ["H1_24_scnorm_counts.txt",
"H1_96_scnorm_counts.txt",
"H9_24_scnorm_counts.txt",
"H9_96_scnorm_counts.txt"]
isoChoiceModels = ["Weibull", "Random", "UniformObserved", "CellVariable"]

function make_k_plot(p1,p2,p3,p4,pReal)

    #initialise mu, sigma and relationships
    mu_arr = Float64[p1[2], p2[2], p3[2], p4[2]]
    sigma_arr = Float64[p1[3], p2[3], p3[3], p4[3]]
    k_arr = Float64[0.25, 0.25, 0.25, 0.25]
    relationships = Array{Float64}(undef, length(pReal[4]), 4)

    df = DataFrame(iteration = Int64[0, 0, 0, 0], NumIsos = Int64[1,2,3,4],
    k = Float64[0.25,0.25,0.25,0.25])

    #Do expectation maximisation for k.
    for i in 1:100
        expectation!(relationships, length(pReal[4]), 4, mu_arr, sigma_arr, k_arr, pReal[4])
        maximisation_k!(relationships, k_arr)

        #save results
        newRow = DataFrame(iteration = repeat([i], 4), NumIsos = [1,2,3,4],
        k = k_arr)
        append!(df, newRow)
    end

    l = @layout [a{0.001h}; b c{0.13w}]

    p1 = Plots.plot(df.iteration, df.k, group = df.NumIsos,
    xlabel = "Iteration", ylabel = "Mixing Fraction", legend = :none,
    axis = (font(20)))

    p2 = Plots.plot(df.iteration, df.k, group = df.NumIsos, grid=false,
    xlims=(-1,-0.5), showaxis=false, legendfontsize = 20)


    p0=plot(title ="", grid=false, showaxis=false)

    #make a plot of k vs iterations to see whether we reached convergence
    out_plot = plot(p0,p1,p2,layout=l, size = (250 * 8, 200 * 8),
    linestyle=:solid, linealpha=0.5, linewidth=2*8, axis = (font(60)),
    legendfontsize=60, margin = 20mm)

    return out_plot
end



for countsMatrix in scnorm_counts

    for model in 1:length(isoChoiceModels)

        #Get the data
        data_store = "figures/figure5_data/" * isoChoiceModels[model] * "_" * countsMatrix
        out_plot_path = "figures/figure5/" * isoChoiceModels[model] * "_" * countsMatrix

        #Extract a plot, mu, sigma and the data
        pReal = GetLogNormal(data_store * "_real_isoforms_array.txt", "Real")
        p1 = GetLogNormal(data_store * "_1_isoforms_array.txt", "1 Isoform")
        p2 = GetLogNormal(data_store * "_2_isoforms_array.txt", "2 Isoforms")
        p3 = GetLogNormal(data_store * "_3_isoforms_array.txt", "3 Isoforms")
        p4 = GetLogNormal(data_store * "_4_isoforms_array.txt", "4 Isoforms")

        #make a plot of k vs iterations to see whether we reached convergence
        out_plot = make_k_plot(p1,p2,p3,p4,pReal)
        savefig(out_plot, out_plot_path)

        #See how good the log normal fit is
        ln_plot = plot([p1[1], p2[1], p3[1], p4[1], pReal[1]]...,
        legend = :none, layout=(5,1),  size = (250*8, 600*8),
        axis = (font(180)), titlefontsize=180, margin = 30mm,
 linestyle=:solid, linealpha=0.5, linewidth=2*10,)
        savefig(ln_plot, out_plot_path * "logNormalFits.png")

    end
end


H1_96_1 =  GetLogNormal(
"figures/figure5_data/Weibull_H1_96_scnorm_counts.txt_1_isoforms_array.txt",
"1 Isoform")
H1_96_2 =  GetLogNormal(
"figures/figure5_data/Weibull_H1_96_scnorm_counts.txt_2_isoforms_array.txt",
"2 Isoforms")
H1_96_3 =  GetLogNormal(
"figures/figure5_data/Weibull_H1_96_scnorm_counts.txt_3_isoforms_array.txt",
"3 Isoforms")
H1_96_4 =  GetLogNormal(
"figures/figure5_data/Weibull_H1_96_scnorm_counts.txt_4_isoforms_array.txt",
"4 Isoforms")
H1_96_real =  GetLogNormal(
"figures/figure5_data/Weibull_H1_96_scnorm_counts.txt_real_isoforms_array.txt",
"Real")

H1_96_kplot = make_k_plot( H1_96_1, H1_96_2, H1_96_3, H1_96_4, H1_96_real)


H1_24_1 =  GetLogNormal(
"figures/figure5_data/Weibull_H1_24_scnorm_counts.txt_1_isoforms_array.txt",
"1 Isoform")
H1_24_2 =  GetLogNormal(
"figures/figure5_data/Weibull_H1_24_scnorm_counts.txt_2_isoforms_array.txt",
"2 Isoforms")
H1_24_3 =  GetLogNormal(
"figures/figure5_data/Weibull_H1_24_scnorm_counts.txt_3_isoforms_array.txt",
"3 Isoforms")
H1_24_4 =  GetLogNormal(
"figures/figure5_data/Weibull_H1_24_scnorm_counts.txt_4_isoforms_array.txt",
"4 Isoforms")
H1_24_real =  GetLogNormal(
"figures/figure5_data/Weibull_H1_24_scnorm_counts.txt_real_isoforms_array.txt",
"Real")

H1_24_kplot = make_k_plot( H1_24_1, H1_24_2, H1_24_3, H1_24_4, H1_24_real)

pAB = plot( H1_96_1[1], H1_24_1[1],
              H1_96_2[1], H1_24_2[1],
              H1_96_3[1], H1_24_3[1],
              H1_96_4[1], H1_24_4[1],
              H1_96_real[1], H1_24_real[1],
              #H1_96_kplot, H1_24_kplot,
              layout=(5,2), size = (500*8, 600*8),
              axis = (font(100)), titlefontsize=100, margin = 30mm,
       linestyle=:solid, linealpha=0.5, linewidth=2*10,)

savefig(pAB, "figures/figure5/figure5AB.png")
savefig(H1_96_kplot, "figures/figure5/figure5C.png")
savefig(H1_24_kplot, "figures/figure5/figure5D.png")
