#This script will remake the figures from the publication
#Note that the order of the figures changed several times during drafting
#so don't pay too much attention to the names of the output pngs.

tar -xvzf data.tar.gz

mkdir figures/figure1
mkdir figures/figure1_data

mkdir figures/figure2
mkdir figures/figure2_data

mkdir figures/figure3
mkdir figures/figure3_data

mkdir figures/figure4
mkdir figures/figure4_data

mkdir figures/figure5
mkdir figures/figure5_data

mkdir figures/figure6
mkdir figures/figure6_data

julia figure1_julia_script.jl
julia figure2_julia_script.jl
julia figure3_julia_script.jl
julia figure4_julia_script.jl
julia figure5_julia_script.jl
julia figure6_julia_script.jl
julia iPredict.jl/src/pDropouts.jl
