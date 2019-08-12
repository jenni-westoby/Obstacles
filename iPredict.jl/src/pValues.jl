
#Get p-values for each number of isoforms for gene. Will this work for global?
function FindpValues(InputDict, mean_val, threshold)

    OutputDict = Dict()

    for (key, value) in InputDict

        #Remove p-values for insufficiently continuous distributions
        if length(unique(value)) < threshold
            OutputDict[key] = missing
        else
            pValue = sum(x -> x > mean_val, value) / length(value)
            OutputDict[key] = pValue
        end
    end

    return OutputDict
end

#Get p-value given an array and mean value
function GlobalpValues(meanArr, mean_val, threshold)

    if length(unique(meanArr)) < threshold
        pValue = missing
    else
        pValue = sum(x -> x > mean_val, meanArr) / length(meanArr)
    end

    return pValue
end
