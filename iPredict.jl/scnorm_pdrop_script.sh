#!/bin/bash

input="scnorm_file_list.txt"

while IFS= read -r line
do
    output=`echo $line | awk -F\/ '{print $NF}' | awk -F_ '{print $1"_"$2".png"}'`
    $HOME/julia-1.1.1/bin/julia src/run_iPredict.jl --mode global --countsMatrix $line --GeneIsoformRelationships ../iPredict_data/gencode_human_gene_isoform_relationships.txt --NumSimulations 10 --output $output
    mv $output ../Obstacles_figs/fig2/$output
    mv "log_"$output ../Obstacles_figs/fig2/"log_"$output
    mv expression_distribution.png ../Obstacles_figs/fig2/$output"_mean_expression.png"
    mv pDropout_distribution.png ../Obstacles_figs/fig2/$output"_pdropout.png"
done < "$input"
