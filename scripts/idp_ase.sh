#!/bin/sh


#prepare Gene level data
julia --depwarn=no  ~/.julia/v0.4/IDPASE/src/prep_runs.jl -a test_data/sim.reads.grch37.unique.psl -g test_data/sim.gpd -v test_data/sim.vcf -q test_data/sim.reads.fq -d /home/benjamin/temp -c 1 -f 1 -o /home/benjamin/test_in/sim

#get commands to run each gene individually
julia phase_by_loci.jl -a /home/benjamin/test_in/ -o /home/benjamin/test_out/ -n SGS -m 1 -d /home/benjamin/temp/ -p sim > to_run_curr.sh

#run each gene 
bash to_run_curr.sh

#concatenate all gene level results
find /home/benjamin/test_out/ -name "REAL*" | xargs cat > gene_results.txt

#prepare isoform level data
julia --depwarn=no  ~/.julia/v0.4/IDPASE/src/prep_runs.jl -a test_data/sim.reads.grch37.unique.psl -g test_data/sim.gpd -v test_data/sim.vcf -q test_data/sim.reads.fq -d /home/benjamin/temp -c 1 -f 1 -o /home/benjamin/test_in/sim -l 100 -i -s -e -r gene_results.txt


#get commands to run each isoform individually
julia phase_isoforms_by_loci.jl -i /home/benjamin/test_in/ -o /home/benjamin/test_out/ -b /home/benjamin/.julia/v0.4/IDPASE/ -a -p sim > to_run_isofs.sh

#run each isoform
bash to_run_isofs.sh

#concatenate all isoform level results
find /home/benjamin/test_out/ -name "EXTRA*" | xargs cat > isoform_results.txt

