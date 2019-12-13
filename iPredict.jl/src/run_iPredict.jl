using ArgParse
include("iPredict.jl")
using .iPredict


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--mode"
            help = "Set to 'global' to do global simulations, set to 'gene' to simulate a single gene"
            required = true
            arg_type = String
        "--gene"
            help = "The gene to simulate in gene-level simulations"
            required = false
            arg_type = String
        "--countsMatrix"
            help = "Path to isoform-level counts matrix"
            required = true
            arg_type = String
        "--GeneIsoformRelationships"
            help = "Path to file containing gene-isoform relationships"
            required = true
            arg_type = String
        "--NumSimulations"
            help = "Number of simulations"
            required = true
            arg_type = Int64
        "--output"
            help="Path to save output at"
            required = true
            arg_type = String
        "--NumIsoformsToSimulate"
            help = "Max number of isoforms to simulate"
            required = false
            arg_type = Int64
        "--filteringThreshold"
            help = "All genes must express at least this number of isoforms"
            required = false
            arg_type = Int64
        "--alpha"
            help = "Beta distribution parameter alpha"
            required = false
            arg_type = Float64
        "--beta"
            help = "Beta distribution parameter beta"
            required = false
            arg_type = Float64
        "--pFP"
            required = false
            default = 0.01
            arg_type = Float64
        "--pFN"
            required = false
            default = 0.04
            arg_type = Float64


    end

    return parse_args(s)
end

function main()

    #Get command line args. Very basic atm, should add some argument validity
    #checking.
    parsed_args = parse_commandline()
    mode = parsed_args["mode"]
    countsMatrix = parsed_args["countsMatrix"]
    GeneIsoformRelationships = parsed_args["GeneIsoformRelationships"]
    NumSimulations = parsed_args["NumSimulations"]
    output = parsed_args["output"]
    NumIsoformsToSimulate = parsed_args["NumIsoformsToSimulate"]
    filteringThreshold = parsed_args["filteringThreshold"]
    alpha::Union{Nothing,Float64} = parsed_args["alpha"]
    beta::Union{Nothing,Float64} = parsed_args["beta"]
    pFP::Float64 = parsed_args["pFP"]
    pFN::Float64 = parsed_args["pFN"]

    if mode == "global"
        globalPredict(countsMatrix, GeneIsoformRelationships,
        NumSimulations, output, NumIsoformsToSimulate, filteringThreshold,
        false, alpha, beta, pFN, pFP)

    elseif mode == "gene"
        gene = parsed_args["gene"]
        oneGenePrediction(countsMatrix, GeneIsoformRelationships, gene,
        NumSimulations, output, false)

    else
        throw(ArgumentError("Unsupported mode - valid modes are gene or global.
        See help manual for details."))
    end
end

main()
