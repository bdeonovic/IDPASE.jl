# IDPASE

Software is to be run on a linux system. To run IDPASE 

1. Install [julia](http://julialang.org/downloads/)
2. From julia's command line run 

  ```julia
  Pkg.clone("git://github.com/bdeonovic/IDPASE.jl.git")
  ```

3. Install [Bedtools](http://bedtools.readthedocs.io/en/latest/content/installation.html) and make sure it is in your PATH.
4. Prepare necessary input files
  * [VCF file](https://en.wikipedia.org/wiki/Variant_Call_Format)
  * [GPD file](https://genome.ucsc.edu/FAQ/FAQformat.html#format9) in the Extended format
  * PSL alignment file of your hybrid-Seq data
  * FASTQ file of your hybrid-Seq data
5. Prepare Gene level data

  ```bash
  julia ~/.julia/v0.4/IDPASE/src/prep_runs.jl -a test_data/SGS.psl test_data/TGS.psl -g test_data/TDRKH.gpd -v test_data/sim.vcf -q test_data/SGS.fq test_data/TGS.fq -d /home/benjamin/temp/ -c 1 -f 1 1 -o /home/benjamin/gene_files/ -p sim
  ```
  where flag ``-a`` is space separated list of PSL files, ``-g`` is GPD file, ``-v`` is VCF file, ``-q`` is space separated list of FASTQ files, ``-d`` is a directory for intermediate output, ``-c`` is a space separated list of chromosomes of interest, ``-f`` is a space separated list of FASTQ formats corresponding to the FASTQ files of ``-q``, where 1 indicates PHRED+33, and 2 indicates PHRED+64, and ``-o`` indicates output directory and output prefix (so in example /out/ is directory and output files will be prefixed by ``sim``). 
6. Get commands to run each gene individually
  ```bash
  julia  ~/.julia/v0.4/IDPASE/scripts/phase_by_loci.jl -a /home/benjamin/gene_files/ -o /home/benjamin/gene_out/ -n SGS TGS -m 1 0 0 1 1 1 -d /home/benjamin/.julia/v0.4/IDPASE/scripts/ -p sim > to_run_curr.sh
  ```
  where ``-a`` is the ``-o`` flag from command in step 5, ``-o`` is an output directory, ``-n`` are unique names corresponding to PSL files, ``-m`` is a vector indicating which combinations of the seq data to use with IDPASE. In the above example 3 runs of IDPASE will be issued where ``1 0`` indicates SGS only, ``0 1`` indicates TGS only, and ``1 1`` indicates hybrid-Seq. ``p`` is the prefix specified in step 5. The output is a list of commands to run for each gene individually. The flag ``-d`` is the directory where the IDPASE scripts are stored. 
7. Run each gene.
  ```bash
  bash to_run_curr.sh
  ```
8. Concatenate all gene level results
  ```bash
   find /home/benjamin/gene_out/ -name "REAL*" | xargs cat > /home/benjamin/gene_out/gene_results.txt
  ```
9. Prepare isoform level data
  
  ```bash
  julia ~/.julia/v0.4/IDPASE/src/prep_runs.jl -a test_data/SGS.psl test_data/TGS.psl -g test_data/TDRKH.gpd -v test_data/sim.vcf -q test_data/SGS.fq test_data/TGS.fq -d /home/benjamin/temp/ -c 1 -f 1 1 -o /home/benjamin/isoform_files/ -p sim -l 100 -i -s -e -r /home/benjamin/gene_out/gene_results.txt
  ```
  where ``-l`` is read length for short reads, ``-s`` to skip file pre-processing (if using same GPD/VCF files are gene level), ``-e`` to use estimated haplotypes from gene leve, otherwise will use information from VCF, asumming it is phased, ``-r`` is the gene level results file. 
10. Get commands to run each isoform individually
  
  ```bash
  julia ~/.julia/v0.4/IDPASE/scripts/phase_isoforms_by_loci.jl -i /home/benjamin/isoform_files/ -o /home/benjamin/isoform_out/ -b /home/benjamin/.julia/v0.4/IDPASE/scripts -a -p sim > to_run_isofs.sh
  ```

11. Run each isoform

  ```bash
  bash to_run_isofs.sh  
  ```

12. concatenate all isoform level results

  ```bash
  find /isoform_out/ -name "EXTRA*" | xargs cat > isoform_results.txt
  ```
