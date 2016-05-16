#!/bin/sh

julia -p 1 --depwarn=no  ~/.julia/v0.4/IDPASE/src/prep_runs.jl -a test_data/sim.reads.grch37.unique.psl -g test_data/sim.gpd -v test_data/sim.vcf -q test_data/sim.reads.fq -d /home/benjamin/temp -c 1 -f 1 -o /home/benjamin/test_out/sim


