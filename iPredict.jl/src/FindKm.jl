using DataFrames
using Statistics
using Rmath
using MultipleTesting
using Optim


################################################################################
# Functions that normalise data like M3Drop expects
################################################################################

#Find no. detected or undetected genes/cell. This is a faff in Julia :(
function FindNumDetectedGenesPerCell(countsMatrix::DataFrame, Detected::Bool)

    #This allows us to use the same function for detected + undetected genes
    if Detected == true
        GreaterThanZero = colwise(x -> x .> 0 , countsMatrix[:,2:ncol(countsMatrix)])
    else
        GreaterThanZero = colwise(x -> x .== 0 , countsMatrix[:,2:ncol(countsMatrix)])
    end

    tmp_arr = Float64[]

    for arr in GreaterThanZero
        push!(tmp_arr, sum(arr))
    end

    outputDf = DataFrame(Cells = names(countsMatrix)[2:ncol(countsMatrix)],
    DetectedGenes = tmp_arr)

    return outputDf

end

#Find mean
function FindMean(DetectedGenesDf::DataFrame)
    return mean(DetectedGenesDf.DetectedGenes)
end

#Find standard deviation
function FindStd(DetectedGenesDf::DataFrame)
    return std(DetectedGenesDf.DetectedGenes)
end


#Remove low quality cells based on no. undetected genes
function RemoveLowQuality(DetectedGenesDf, mean, sd)

    pnorms = (DetectedGenesDf.DetectedGenes .- mean) ./ sd
    pnorms = pnorm.(pnorms, false)
    padjust = adjust(pnorms, BenjaminiHochberg())

    DetectedGenesDf.pvals = padjust
    filtered_df = DetectedGenesDf[DetectedGenesDf.pvals .> 0.05, :]

    return(filtered_df[:, 1:2])
end

#Only keep filtered cells
function KeepFiltered(DetectedGenesDf::DataFrame, countsMatrix::DataFrame)

    arr = [Symbol("Isoforms")]
    append!(arr, DetectedGenesDf.Cells)

    return countsMatrix[arr]

end

#Helper function to make an empty counts matrix
function MakeEmptyDf(countsMatrix)

    #Make an empty df to store results in
    namelist = names(countsMatrix)
    df = DataFrame()

    for (i, name) in enumerate(namelist)
        if i == 1
            df[name] = String[]
        else
            df[name] = Float64[]
        end
    end

    return df
end

#Remove undetected genes (separate function)
function RemoveUndetectedGenes(countsMatrix::DataFrame)

    #Make an empty df to store results in
    df = MakeEmptyDf(countsMatrix)

    #Iterate over rows, keep rows where isoforms is expressed in >3 cells
    for row in eachrow(countsMatrix)

        if sum(row[2:ncol(countsMatrix)] .> 0 ) > 3
            push!(df, row)
        end
    end

    return df
end


#Convert to cpm
function ConvertToCpm(countsMatrix::DataFrame)

    #Find cpms
    cpms = colwise(x -> (x / sum(x)) * 1000000, countsMatrix[2:ncol(countsMatrix)])

    #Make back into dataframe...
    namelist = names(countsMatrix[2:ncol(countsMatrix)])
    df = DataFrame(Isoforms = countsMatrix.Isoforms)

    for (i, name) in enumerate(namelist)
        df[name] =  cpms[i]
    end

    return df

end

#Remove lowly expressed genes (function)
function RemoveLowlyExpressedGenes(countsMatrix::DataFrame)

    #Make an empty df to store results in
    df = MakeEmptyDf(countsMatrix)

    #Iterate over rows, keep rows where isoforms isn't lowly expressed
    for row in eachrow(countsMatrix)

        if sum(row[2:ncol(countsMatrix)]) > 0.00001
            push!(df, row)
        end
    end

    return df
end

#Make a function to join everything up
function CleanDataForMM(countsMatrix::DataFrame)

    #Find no. undetected genes/cell
    NumZero = FindNumDetectedGenesPerCell(countsMatrix, false)

    #Find mean and sd no. undetected genes per cell
    Mu = FindMean(NumZero)
    Sigma = FindStd(NumZero)

    #Remove low quality cells based on no. undetected genes per cell
    highQualityCells = RemoveLowQuality(NumZero, Mu, Sigma)

    #Only keep high quality cells in countsMatrix
    outputDf = KeepFiltered(highQualityCells, countsMatrix)

    #Remove undetected genes
    outputDf = RemoveLowlyExpressedGenes(outputDf)

    #Convert counts to cpm
    outputDf = ConvertToCpm(outputDf)

    #Remove lowly expressed genes
    outputDf = RemoveLowlyExpressedGenes(outputDf)

    return outputDf

end

################################################################################
# Functions that calculate Km
################################################################################

#Calculate p (M-M parameter)
function Findp(countsMatrix::DataFrame)

    df = DataFrame()
    df[:Isoforms] = countsMatrix[:Isoforms]
    arr = Float64[]

    for row in eachrow(countsMatrix)
        p = 1 - (sum(row[2:ncol(countsMatrix)] .> 0) / (ncol(countsMatrix) - 1))
        push!(arr, p)
    end

    df[:p] = arr
    return df

end

#Calculate s (M-M parameter)
function Finds(df::DataFrame, countsMatrix::DataFrame)

    arr = Float64[]

    for row in eachrow(countsMatrix)
        s = mean(row[2:ncol(countsMatrix)])
        push!(arr, s)
    end

    df[:s] = arr
    return df
end

#Function that makes function to calculate log likelihood
function LogLikelihood(p, s)
    f(params) = -sum(dnorm.(p .- (1.0 .- (s ./ (params[1] .+ s))), 0.0, params[2], true))
    f
end

#MLE
function MaxLikelihoodEst(p, s)

    #Check p and s make sense
    if length(p) != length(s)
        throw(ArgumentError("Length of p and s don't match, something went wrong\n"))
    end

    #NB. MLE starts get inconsistent below 5,000 genes
    if length(p) < 5000
        println("Warning: Calculating Km for less than 5000 isoforms... MLE
        estimate of Km may not be accurate.")
    end

    #This line basically internalises p and s - useful for optimize function below
    LogLikelihoodFn = LogLikelihood(p,s)

    #Perform MLE and return parameters
    res = optimize(LogLikelihoodFn, [3.0, 0.25], method = NelderMead())
    parameters = Optim.minimizer(res)
    return parameters
end

#Function to join up all the other functions
function FindKm(countsMatrix::DataFrame)

    df = Findp(countsMatrix)
    df = Finds(df, countsMatrix)
    parameters = MaxLikelihoodEst(df.p, df.s)
    return parameters
end
