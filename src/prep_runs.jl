using ArgParse
using IDPASE

function parse_commandline()
  s = ArgParseSettings()
  s.description = "ASE IDP (Allele Specific Expression, Isoform detection and prediction)"
  s.commands_are_required = false
  s.version = "1.0"
  s.add_version = true

  @add_arg_table s begin
    "--psl", "-a"
      arg_type = ASCIIString
      nargs = '+'
      help = "sequence alignment files in PSL format"
      required = true
    "--fastq", "-q"
      arg_type = ASCIIString
      nargs = '+'
      help = "FASTQ files. Order matches order in PSL_FILES"
      required = true
    "--gpd", "-g"
      arg_type = ASCIIString
      help = "GPD file"
      required = true
    "--vcf", "-v"
      arg_type = ASCIIString
      help = "VCF file"
      required = true
    "--temp", "-d"
      arg_type = ASCIIString
      default = "temp"
      help = "Temporary directory"
    "--chr", "-c"
      arg_type = ASCIIString
      nargs = '+'
      help = "Chromosome Names"
      required = true
    "--format","-f"
      arg_type = Int
      nargs = '+'
      help = "Phred format"
    "--out","-o"
      arg_type = ASCIIString
      help = "output dir"
    "--skip","-s"
      action = :store_true
      help = "Skip initial processing"
    "--simulate","-u"
      action = :store_true
      help = "Semi Simulation"
    "--isoform", "-i"
      action = :store_true
      help = "Isoform level estimation"
    "--results", "-r"
      arg_type = ASCIIString
      help = "MCMC results"
    "--phase_psl", "-p"
      arg_type = ASCIIString
      help = "PSL for isoform specific ASE"
    "--type", "-t"
      arg_type = ASCIIString
      default = "SGSTGS"
    "--estimate", "-e"
      action = :store_true
    "--fpkm", "-k"
      arg_type = ASCIIString
    "--read_length", "-l"
      arg_type = Int
      default = 100
    "--only_sim", "-n"
      action = :store_true
    "--print", "-m"
      action = :store_true
  end

  return parse_args(s)
end

function main()
  parsed_args = parse_commandline()
  println("Parsed args:")
  for (arg,val) in parsed_args
    println("  $arg  =>  $val")
  end

  if length(parsed_args["psl"]) != length(parsed_args["fastq"])
    throw(ArgumentError("Number of PSL files needs to match number of FASTQ files"))
  end

  if !parsed_args["skip"]
    println(STDOUT,"Processing GPD...")
    @time gpd_to_bed(parsed_args["gpd"], chr = parsed_args["chr"], output_dir=parsed_args["temp"])
    isoform_dict = process_gpd(parsed_args["gpd"])

    println(STDOUT,"Processing VCF...")
    @time vcf_to_bed(parsed_args["vcf"], chr = parsed_args["chr"], output_dir=parsed_args["temp"])

    num_src_types = length(parsed_args["psl"])
    num_reads = zeros(num_src_types)
    println(STDOUT,"Processing PSL files...")
    for i in 1:num_src_types
      @time num_reads[i] = psl_to_bed2(parsed_args["psl"][i],chr = parsed_args["chr"], output_dir = parsed_args["temp"],output_file_name="psl$i.bed")
    end
    if parsed_args["phase_psl"] != nothing
      @time num_reads = psl_to_bed2(parsed_args["phase_psl"], chr = parsed_args["chr"], output_dir = parsed_args["temp"],output_file_name="psl_phase.bed")
    end
  end

  args = collect(zip(parsed_args["chr"], 
                    [parsed_args for i in 1:length(parsed_args["chr"])], 
                    [parsed_args["simulate"] for i in 1:length(parsed_args["chr"])],
                    [parsed_args["isoform"] for i in 1:length(parsed_args["chr"])] ))
  @time pmap(make_X_Q, args)
end

main()
