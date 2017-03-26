using ArgParse
using IDPASE
using Mamba

function parse_commandline()
  s = ArgParseSettings()
  s.description = "ASE IDP (Allele Specific Expression, Isoform detection and prediction)"
  s.commands_are_required = false
  s.version = "0.0.1"
  s.add_version = true

  @add_arg_table s begin
    "--true", "-t"
      arg_type = String
    "--reads", "-a"
    "--iters", "-i"
      arg_type = Int
      default = 10000
      help = "Number of MCMC iterations"
    "--burnin", "-b"
      arg_type = Int
      default = 1000
      help = "Burnin"
    "--chains", "-c"
      arg_type = Int
      default = 4
      help = "Number of chains"
    "--temp", "-d"
      arg_type = String
      default = "temp"
      help = "Temporary directory"
    "--out", "-o"
      arg_type = String
      default = "out.txt"
      help = "Output directory"
    "--chr", "-r"
      arg_type = String
      default = "chr1"
      help = "Chromosome"
    "--method","-e"
      arg_type = Int
      default = 1
      help = "Default MCMC method"
    "--line", "-l"
      arg_type = Int
      default = 1
      help = "loci line"
    "--simulate", "-u"
      action = :store_true
      help = "Semi Simulation"
    "--names", "-n"
      arg_type = String
      nargs = '+'
      help = "Source names"
      default = ["SGS", "SNY", "SUB"]
    "--matrix", "-m"
      arg_type = Int
      nargs = '+'
      help = "Source matrix"
      default = [1, 0, 0, 1, 1, 0, 0, 1, 0]
    "--schedule", "-s"
      arg_type = Float64
      nargs = '+'
      default = [1.0, 0.5, 0.25, 0.125]
  end

  return parse_args(s)
end

function main()
  parsed_args = parse_commandline()
  println("Parsed args:")
  for (arg,val) in parsed_args
    println("  $arg  =>  $val")
  end

  data = read_data(parsed_args["true"], parsed_args["reads"], parsed_args["line"])
  source_names = parsed_args["names"]
  source_matrix = reshape(convert(Array{Bool,1}, parsed_args["matrix"]), length(source_names), round(Int, length(parsed_args["matrix"])/length(source_names)))

  num_src, num_combs = size(source_matrix)

  
  println(STDERR,"Running MCMC...")
  for j in 1:num_combs, k in parsed_args["schedule"]
    use_src = source_matrix[:,j]
    results = run_MCMC((data, parsed_args["iters"], parsed_args["burnin"], parsed_args["chains"], parsed_args["method"], use_src, parsed_args["names"], k))
    new_gene_name = replace(data.gene_name, r"/", s"|")
    prefix = string(parsed_args["chr"],"::",parsed_args["line"],"::",new_gene_name)
    if parsed_args["simulate"]
      prefix = string("SIM::",prefix)
    else
      prefix = string("REAL::",prefix)
    end
    comb = join(source_names[use_src],"")
    OUT = open(string(parsed_args["out"],prefix,"::","$(comb)::$(k).txt"),"w")
    #writedlm(OUT,results')
    writedlm(OUT, reshape(results, 1, length(results)))
    close(OUT)
  end
end

main()
