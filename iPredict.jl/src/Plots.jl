#using Gadfly
using DataFrames
using Statistics
import Cairo, Fontconfig
using StatsBase
using LinearAlgebra
using Plots
import GR
using Plots.PlotMeasures

#Function to convert dictionary to dataframe
function DictToDf(InputDict)
    return DataFrame(; [Symbol(k)=>v for (k,v) in InputDict]...)
end

#Find mean number of isoforms per gene per cell in real data
function MeanIsoPerGenePerCell(singleGeneCounts::DataFrame)

    arr = []

    for i in 3:ncol(singleGeneCounts)
        isoPerGenePerCell = sum(x->x>0, singleGeneCounts[:, i])
        push!(arr, isoPerGenePerCell)
    end

    return mean(arr)
end

#Histogram plot with real mean value marked on it for a single gene
function OneGeneHistPlot(InputDict, singleGeneCounts, output)

    #Prep df for histogram part of plot
    df = DictToDf(InputDict)
    df = stack(df, 1:ncol(df))

    #Find mean for line on plot
    mean_val = MeanIsoPerGenePerCell(singleGeneCounts)

    test = Gadfly.plot(df, x = "value", color = "variable", xintercept = [mean_val],
    Geom.histogram(position=:dodge), Geom.vline(color = "black"), Guide.xlabel("Mean isoforms per gene per cell"),
    Guide.ylabel("Frequency"), Guide.colorkey(title = "Number of Isoforms"))

    img = PNG(output)
    draw(img, test)

end

#Function that finds the means number of isoforms per gene
#This has to be re-written due to a bullshit memory bug
#Function that finds the means number of isoforms per gene
function MeanNumIsoformsPerGenePerCellReal(countsDict)

    meansArr = Float64[]
    meansDict = Dict()

    for (key, value) in countsDict

        tmp_df = DictToDf(value)
        cellArr =[]

        for row in eachrow(tmp_df)
            cellMean = sum(row .> 0)
            push!(cellArr, cellMean)
        end

        push!(meansArr, mean(cellArr))
        meansDict[key] = mean(cellArr)

    end

    return [meansArr, meansDict]
end

function makeNormalisedHistogram(Arr, titleArg)

    #Do some magic to make the histogram as x and y coordinates
    hist = fit(Histogram, Arr, nbins=100, closed=:right);
    hist = normalize(hist, mode=:pdf)

    if titleArg == "3 Isoforms"
        #Untransformed plot
        p1 = histogram(Arr,
        title = titleArg,
        xlims = (0, 4),
        ylims = (0, 15),
        ylabel = "Normalised Frequency",
        normalize=:pdf,
        bins = 100,
        axis = (font(18)), titlefontsize=18, yticks = 0:5:15, left_margin=1.0mm)

    else
        #Untransformed plot
        p1 = histogram(Arr,
        title = titleArg,
        xlims = (0, 4),
        ylims = (0, 15),
        normalize=:pdf,
        bins = 100,
        axis = (font(18)), titlefontsize=18, yticks = 0:5:15, left_margin=1.0mm)
    end

    #Log transformed plot
    p2 = plot(hist.edges[1][2:end],
    log10.(hist.weights.+1),
    title = titleArg,
    xlims = (0, 4),
    ylims = (0, 2),
    ylabel = "Log Normalised Frequency")

    return (p1, p2)
end

function makeNormalisedHistogramReal(Arr, titleArg)

    #Do some magic to make the histogram as x and y coordinates
    hist = fit(Histogram, Arr, nbins=100, closed=:right);
    hist = normalize(hist, mode=:pdf)

    println("Making real histogram")

    p1 = histogram(Arr,
    title = titleArg,
    xlims = (0, 4),
    xlabel = "Mean Isoforms per Gene per Cell",
    ylims = (0, 15),
    normalize=:pdf,
    bins = 100,
    axis = (font(18)), titlefontsize=18, yticks = 0:5:15, left_margin=1.0mm,
    bottom_margin = 1.0mm)

    #Log transformed plot
    p2 = plot(hist.edges[1][2:end],
    log10.(hist.weights.+1),
    title = titleArg,
    xlims = (0, 4),
    xlabel = "Mean Isoforms per Gene per Cell",
    ylims = (0, 2),
    ylabel = "Log Normalised Frequency",
    normalize=:pdf)

    return (p1, p2)
end


function makeHistogramforNIsoforms(df, condition, titleArg)

    #Extract data from dataframe and make plot
    get_data = df[df.type .== condition, :]
    p1 = makeNormalisedHistogram(get_data.means, titleArg)
    return p1
end


function GlobalGeneDensityPlot(globalMeans, realArr, output)

    prinln("Making plots")

    #Make plots
    pReal = makeNormalisedHistogramReal(realArr, "Real")
    p1 = makeHistogramforNIsoforms(globalMeans, 1, "One Isoform")
    p2 = makeHistogramforNIsoforms(globalMeans, 2, "Two Isoforms")
    p3 = makeHistogramforNIsoforms(globalMeans, 3, "Three Isoforms")
    p4 = makeHistogramforNIsoforms(globalMeans, 4, "Four Isoforms")

    #Make panel of normalised histograms
    graph = plot(p1[1],p2[1],p3[1],p4[1],pReal[1], layout=(5,1),legend=false, size = (300, 600) )
    savefig(graph, output)

    #Make panel of normalised log histograms
    log_graph = plot(p1[2],p2[2],p3[2],p4[2],pReal[2],layout=(5,1),legend=false, size = (300, 600))
    log_out = "log_" * output
    savefig(log_graph, log_out)

end

function pDropoutDistribution(Arr, output)

    #Default dropout distribution
    p1 = histogram(Arr, title = "Dropout Distribution", axis = (font(70)),
    titlefontsize=70, margin = 10.0mm)

    #Beta distribution
    params = fit(Beta, Arr)
    @show params

    #plot Beta distribution
    p2 = plot(0:0.01:1, pdf(Beta(params.α, params.β), 0:0.01:1),
    title = "Beta Distribution", axis = (font(70)), titlefontsize=70,
    margin = 10.0mm, linestyle=:solid, linealpha=0.5, linewidth=2*8,
    xlabel = "p(Dropout)")

    #merge plots
    p_all = plot(p1,p2, layout = (2,1))

    #Save
    savefig(p_all, output)
    return (p1, p2)
end

function meanDistribution(Arr)
    p1 = histogram(Arr, title = "Mean Expression Distribution")
    savefig(p1, "expression_distribution.png" )
end
