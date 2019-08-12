using DelimitedFiles
using Distributions
using StatsPlots
using Plots.PlotMeasures

#function that reads in data and return log normal parameters + plot
function GetLogNormal(path, titleArg)
    tmp = readdlm(path)
    params = fit(LogNormal, tmp)
    mu = params.μ
    sigma = params.σ
    if titleArg == "Real"
        p1 = histogram(tmp, normalize = :pdf, xlims = (0, 4), ylims = (0,10),
        title = titleArg, xlabel = "Mean Isoforms per Gene per Cell",
        axis = (font(18)), titlefontsize=18, yticks = 0:5:15, left_margin=1.0mm,
        legend = false)
    elseif titleArg == "3 Isoforms" && path != "figures/figure5_data/Weibull_H1_24_scnorm_counts.txt_3_isoforms_array.txt"
        p1 = histogram(tmp, normalize = :pdf, xlims = (0, 4), ylims = (0,10),
        title = titleArg, ylabel = "Normalised Frequency", axis = (font(18)),
        titlefontsize=18, yticks = 0:5:15, left_margin=1.0mm, legend = false)
    else
        p1 = histogram(tmp, normalize = :pdf, xlims = (0, 4), ylims = (0,10),
        title = titleArg, axis = (font(18)), titlefontsize=18, yticks = 0:5:15,
        left_margin=1.0mm, legend = false)
    end
    plot!(LogNormal(mu, sigma), linewidth = 2)

    return [p1, mu, sigma, tmp]
end

################################################################################
# Implementation of expectation maximisation - we are only looking for k
function find_relationships(mu_arr, sigma_arr, k_arr, k, x)
    top = k_arr[k] * pdf(LogNormal(mu_arr[k], sigma_arr[k]), x)

    bottom = 0.0

    for i in 1:length(mu_arr)
        bottom += k_arr[i] * pdf(LogNormal(mu_arr[i], sigma_arr[i]), x)
    end

    #@show top/bottom
    return top/bottom
end

function expectation!(relationships, dim1, dim2, mu_arr, sigma_arr, k_arr, real_arr)
    for point in 1:dim1
        for k in 1:dim2
            #@show point
            #@show k
            relationships[point, k] = find_relationships(mu_arr, sigma_arr,
            k_arr, k, real_arr[point])
        end
    end
end


function maximisation_k!(relationships, k_arr)

    for i in 1:length(k_arr)
        k_arr[i] = sum(relationships[:, i]) / length(relationships[:,i])
    end
end

################################################################################
#These functions are not for actual use but are convenient for testing
function maximisation_mu!(relationships, mu_arr, x)

    for i in 1:length(mu_arr)
        mu_arr[i] = (relationships[:, i] .* x) ./ (relationships[:, i])
    end
end

function maximisation_sigma!(relationships, sigma_arr, mu_arr, x)

    for i in 1:length(sigma_arr)
        sigma_arr[i] = (relationships[:, i] .* (x .- mu_arr[i]).^2) ./
        relationships[:, i]
    end
end
