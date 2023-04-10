#!/bin/bash

for x in $(seq 3 4)
do
    sed -i -e "s/trG_.*_Fe_.*_Co_.*.dat/trG_up_Fe_${x}_Co_$((10 - x)).dat/" -e "s/spectrum_.*_Fe_.*_Co_.*.dat/spectrum_up_Fe_$((x*10))_Co_$((10*(10 - x))).dat/" nevanlinna.jl

    ~/julia/julia nevanlinna.jl

 sed -i -e "s/trG_.*_Fe_.*_Co_.*.dat/trG_down_Fe_${x}_Co_$((10 - x)).dat/" -e "s/spectrum_.*_Fe_.*_Co_.*.dat/spectrum_down_Fe_$((x*10))_Co_$((10*(10 - x))).dat/" nevanlinna.jl

    ~/julia/julia nevanlinna.jl

done
