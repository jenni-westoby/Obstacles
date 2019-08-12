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

julia figure1_julia_script.jl
julia figure2_julia_script.jl
julia figure3_julia_script.jl
julia figure4_julia_script.jl
julia figure5_julia_script.jl
julia iPredict.jl/src/pDropouts.jl
