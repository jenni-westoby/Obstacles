using DataFrames
using CSV

#Function that reads in counts matrix
function ReadCountsMatrix(filename)
    matrix = CSV.read(filename, header = true, delim = '\t')
    return matrix
end

#Function that reads in gene-isoform table
function ReadGeneIsoformRelationships(filename)
    relationships = CSV.read(filename, header = false, delim = ' ')
    return relationships
end


#Function that checks gene-isoform table matches counts matrix
function CheckGenesMatchIsoforms(countsMatrix, GeneIsoformRelationship)
    arr1 = sort(countsMatrix[:Column1])
    arr2 = sort(GeneIsoformRelationship[:Column2])
    if arr1 == arr2
        return true
    end
    return false
end

#Function that adds gene information to isoform table
function AddGeneInfo(countsMatrix, GeneIsoformRelationship)

    #Call the columns something sensible
    rename!(countsMatrix, :Column1 => :Isoforms)
    rename!(GeneIsoformRelationship, :Column1 => :Genes)
    rename!(GeneIsoformRelationship, :Column2 => :Isoforms)

    return join(GeneIsoformRelationship, countsMatrix, on = :Isoforms)
end


#Function that converts dataframe to dict, removes unexpressed isoforms and genes,
#and returns a dict of how many isoforms are expressed per gene.
function KeepGenesWithNIsoforms(countsMatrix, N)

    #Check N is set to a sensible value
    if N<1
        throw(ArgumentError("You can't simulate the expression of less than one
        isoform!\n"))
    end

    if N>4
        println("You have asked to simulate more than 4 isoforms. Not many genes
        have more than 4 isoforms, so this might violate the assumptions of some
         of iPredict's statistical models.")
    end

    #sort countsMatrix. The rationale behind this is that in theory we can make
    #this closer to O(n) instead of O(n^2)
    sort!(countsMatrix, (:Genes))

    #Set first gene
    gene = ""

    #Make empty dicts
    local dict_to_add = Dict{String, Array{Float64, 1}}()
    local output_dict = Dict{String, Any}()
    local unexpr_dict = Dict{String, Int64}()
    local NumIsos = 0
    local MaxNumIsos = 0
    local UnexprIsos = 0

    #Iterate once over countsMatrix
    for i in 1:nrow(countsMatrix)

        row = countsMatrix[i,:]

        #Check gene is still the same
        if row.Genes == gene

            #if N > numisos or N < numalreadyexpressed, continue
            if N > NumIsos
                continue
            end
            if N < length(dict_to_add)
                continue
            end

            #See whether row passes cutoff
            cutoff = sum(row[3:ncol(countsMatrix)].>=5)

            if cutoff >= 2
                #Add isoform to dictionary
                dict_to_add[row.Isoforms] = Vector(row[3:ncol(countsMatrix)])
            end

        else #new gene

            if length(dict_to_add) == N #See whether the previous gene passes
                output_dict[gene] = dict_to_add
                unexpr_dict[gene] = NumIsos - N

                if NumIsos > MaxNumIsos
                    MaxNumIsos = NumIsos
                end
            end

            #Reset gene, NumExpressed and dict_to_add
            gene = countsMatrix.Genes[i]
            dict_to_add = Dict{String, Array{Float64, 1}}()


            #Find NumIsos
            NumIsos = 0
            j = deepcopy(i)
            while countsMatrix.Genes[j] == gene
                j += 1
                NumIsos += 1
                #Make sure we don't go beyond the end of the matrix
                if j > nrow(countsMatrix)
                    break
                end
            end

            if NumIsos < N
                continue
            end

            #See whether row passes cutoff
            cutoff = sum(row[3:ncol(countsMatrix)].>=5)

            if cutoff >= 2
                #Add isoform to dictionary
                dict_to_add[row.Isoforms] = Vector(row[3:ncol(countsMatrix)])
            end
        end

        #Catch last gene in table
        if i == nrow(countsMatrix)
            if length(dict_to_add) == N
                output_dict[gene] = dict_to_add
            end
        end
    end

    #Return dictionary
    return (output_dict, unexpr_dict, MaxNumIsos)

end

#Function that converts expressed isoforms for a single gene to a dict. Used by
#oneGenePrediction()
function ConvertExpressedIsosToDict(countsMatrix)

    #Make empty dicts
    local dict_to_add = Dict{String, Array{Float64, 1}}()
    local output_dict = Dict{String, Any}()

    for row in eachrow(countsMatrix)

        #See whether row passes cutoff
        cutoff = sum(row[3:ncol(countsMatrix)].>=5)

        if cutoff >= 2
            #Add isoform to dictionary
            dict_to_add[row.Isoforms] = Vector(row[3:ncol(countsMatrix)])
        end
    end

    #Add gene info and return dictionary
    output_dict[countsMatrix.Genes[1]] = dict_to_add
    return output_dict
end
