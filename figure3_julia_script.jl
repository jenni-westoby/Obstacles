include("iPredict.jl/src/iPredict.jl")
using .iPredict
using DataFrames
using Plots
using JLD2
using DelimitedFiles

#Initialise variables
GeneIsoformRelationships = "data/gencode_human_gene_isoform_relationships.txt"
NumSimulations = 100
NumIsoformsToSimulate = 4
filteringThreshold = 4
alpha = nothing
beta = nothing
scnorm_counts = ["H1_24_scnorm_counts.txt",
"H1_96_scnorm_counts.txt",
"H9_24_scnorm_counts.txt",
"H9_96_scnorm_counts.txt"]
error = collect(0:0.1:0.5)

df = DataFrame(counts=String[], pFP = Float64[], pFN = Float64[],
NumIsos = Int64[], mean = Float64[], std = Float64[], median = Float64[])

for countsMatrix in scnorm_counts

    for i in 1:length(error) # i = pFN
        for j in 1:length(error) # j = pFP

            #output directories
            output_both = "figures/figure3/" * string(error[i]) * "_" * string(error[j]) * "_" * countsMatrix
            output_both_data = "figures/figure3_data/" * string(error[i]) * "_" * string(error[j]) * "_" * countsMatrix

            #simulate
            out = globalPredict("data/" * countsMatrix, GeneIsoformRelationships,
            NumSimulations, output_both, NumIsoformsToSimulate, filteringThreshold,
            false, alpha, beta, error[i], error[j], "none", "Weibull", output_both_data)

            newRow = DataFrame(counts = repeat([countsMatrix], 4),
            pFP = repeat([error[j]], 4), pFN = repeat([error[i]], 4),
            NumIsos = [1,2,3,4], mean = out[1], std = out[2], median = out[3])
            append!(df, newRow)

        end

    end
    # globalPredict("data/" * countsMatrix, GeneIsoformRelationships,
    # NumSimulations, "figures/figure3/random_" * countsMatrix , NumIsoformsToSimulate, filteringThreshold,
    # false, alpha, beta, 0.04, 0.01, "random", "Weibull", "none")
    # globalPredict("data/" * countsMatrix, GeneIsoformRelationships,
    # NumSimulations, "figures/figure3/IsoformDependence_" * countsMatrix , NumIsoformsToSimulate, filteringThreshold,
    # false, alpha, beta, 0.04, 0.01, "IsoformDependence", "Weibull", "none")
end


jldopen("figures/figure3_data/means_df.jld2", "w") do file
    file["df"] = df
end




# p1 = Plots.plot(df.error, df.mean, markershapes =[:circle, :star5, :x, :square], group = df.NumIsos)
# savefig(p1, "figures/figure3/mean_scnorm_err.png")
# p1 = Plots.plot(df.error, df.median, markershapes =[:circle, :star5, :x, :square], group = df.NumIsos)
# savefig(p1, "figures/figure3/median_scnorm_err.png")
# p1 = Plots.plot(df.error, df.std, markershapes =[:circle, :star5, :x, :square], group = df.NumIsos)
# savefig(p1, "figures/figure3/std_scnorm_err.png")


using DataFrames
using Plots, Plots.Measures
using JLD2
using DelimitedFiles

error = collect(0:0.1:0.5)

function make_hist(path)
    Arr = readdlm(path)

    #Untransformed plot
    p1 = histogram(Arr,
    xlims = (0, 4),
    ylims = (0, 5),
    normalize=:pdf,
    bins = 100,
    axis = (font(80)), titlefontsize=80, yticks = 0:5:10, xticks = 0:1:4,
    margin = 25mm)

    return p1
end


plot_arr = []

for i in 1:length(error) # i = pFN
    for j in 1:length(error) # j = pFP

        push!(plot_arr, make_hist("figures/figure3_data/" *
        string(error[i]) * "_" * string(error[j]) * "_" *
        "H1_24_scnorm_counts.txt_1_isoforms_array.txt"))
    end
end

tmp = plot(plot_arr..., layout =(6,6),
size = (720 * 8, 720 * 8), legend = false, margin = 25mm)

savefig(tmp, "figures/figure3/figure3A.png")

c = jldopen("figures/figure3_data/means_df.jld2", "r")
df = c["df"]


both = df[df.pFP .== df.pFN, :]
pFP = df[df.pFN .== 0.00, :]
pFN = df[df.pFP .== 0.00, :]

both_arr =[]
pFN_arr =[]
pFP_arr =[]

for i in 1:4
    tmp = both[both.NumIsos .== i, :]
    both_plot = Plots.plot(tmp.pFP, tmp.mean, group = tmp.counts, title = string(i))
    push!(both_arr, both_plot)

    tmp = pFN[pFN.NumIsos .== i, :]
    pFN_plot = Plots.plot(tmp.pFN, tmp.mean, group = tmp.counts, title = string(i))
    push!(pFN_arr, pFN_plot)

    tmp = pFP[pFP.NumIsos .== i, :]
    pFP_plot = Plots.plot(tmp.pFP, tmp.mean, group = tmp.counts, title = string(i))
    push!(pFP_arr, pFP_plot)
end

p1 = plot(both_arr[1], legend = :none, xlabel = "pFP and pFN", ylabel = "mean",
       size = (180 * 8, 180*8),titlefontsize = 80, title = "1 Isoform",
       linestyle=:solid, linealpha=0.5, linewidth=2*8, axis = font(80), margin = 10mm)
savefig(p1, "figures/figure3/both_means.png")
p1 = plot(pFN_arr[1], legend = :none, xlabel = "pFN", ylabel = "mean",
       size = (180 * 8, 180*8),titlefontsize = 80, title = "1 Isoform",
       linestyle=:solid, linealpha=0.5, linewidth=2*8, axis = font(80), margin = 10mm)
savefig(p1, "figures/figure3/pFN_means.png")
p1 = plot(pFP_arr[1], legend = :none, xlabel = "pFP", ylabel = "mean",
       size = (180 * 8, 180*8),titlefontsize = 80, title = "1 Isoform",
       linestyle=:solid, linealpha=0.5, linewidth=2*8, axis = font(80), margin = 10mm)
savefig(p1, "figures/figure3/pFP_means.png")

#get a legend
#get a legend
       p1 = plot(both_arr[1], legend = :inside,
              labels = ["H1 4 million reads", "H1 1 million reads", "H9 4 million reads",
              "H9 1 million reads"], size = (180 * 8, 180*8),titlefontsize = 50, legendfontsize=80, grid=false,
                  xlims=(-1,-0.5), showaxis=false,linestyle=:solid, linealpha=0.5, linewidth=2*8, xlabel = "", ylabel = "", title="")
savefig(p1, "figures/figure3/means_legend.png")
