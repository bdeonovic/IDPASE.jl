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
    "--prefix", "-p"
      arg_type = String
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
    "--in", "-a"
      arg_type = String
    "--out", "-o"
      arg_type = String
      default = "out.txt"
      help = "Output directory"
    "--names", "-n"
      arg_type = String
      nargs = '+'
      help = "Source names"
      default = ["SGS", "SNY", "SUB"]
    "--matrix", "-m"
      arg_type = Int
      nargs = '+'
      help = "Source matrix"
      default = [1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1]
    "--schedule", "-s"
      arg_type = Float64
      nargs = '+'
      default = [1.0]# [1.0, 0.5, 0.25, 0.125]
  end

  return parse_args(s)
end

function main()
  parsed_args = parse_commandline()
  println(STDERR, "Parsed args:")
  for (arg,val) in parsed_args
    println(STDERR, "  $arg  =>  $val")
  end
  matrix_comb = reshape(parsed_args["matrix"], length(parsed_args["names"]), round(Int64,length(parsed_args["matrix"])/length(parsed_args["names"])))
  name_comb = Array(String, size(matrix_comb,2))

  for i in 1:length(name_comb)
    name_comb[i] = join(parsed_args["names"][find(matrix_comb[:,i])])
  end


  real_completed = Dict{String, Array{Bool,1}}()
  real_data = Dict{String, Array{Bool,1}}()

  sim_completed = Dict{String, Array{Bool,1}}()
  sim_data = Dict{String, Array{Bool,1}}()
  for file in readdir(parsed_args["in"])
    m = match(Regex("$(parsed_args["prefix"])_(true|sim)_([\\w\\[\\]\\/\\.\\-]*).txt"),file)
    if m != nothing
      num_loci = parse(Int64, split(readstring(`wc -l $(parsed_args["in"])/$file`))[1])
      if m.captures[1] == "true"
        real_completed[m.captures[2]] = falses(num_loci)
      elseif m.captures[1] == "sim"
        sim_completed[m.captures[2]] = falses(num_loci)
      end
    end
  end
  for file in readdir(parsed_args["out"])
    m = match(Regex("(REAL|SIM)::([\\w\\[\\]\\/\\.\\-]*)::(\\d*)::([\\w\\[\\]\\/\\.\\-]*)::($(join(name_comb,"|")))::(\\d*).txt"),file)
    if m != nothing
      sim = m.captures[1]
      chr = m.captures[2]
      line_num = parse(Int64, m.captures[3])
      gene_name = m.captures[4]
      gene_name = replace(gene_name, r"/", s"|")
      data_type = m.captures[5]
      subsample = parse(Float64, m.captures[6])
      if sim == "REAL"
        if !haskey(real_data, gene_name)
          real_data[gene_name] = falses(length(name_comb) * length(parsed_args["schedule"]))
        end
        if data_type in name_comb
          idx1 = findfirst((x)-> x == data_type,name_comb)
          idx2 = findfirst((x)-> x == subsample,parsed_args["schedule"])
          real_data[gene_name][(idx2-1)*length(name_comb) + idx1] = true
        else
          warn("unknown data type")
        end
        if sum(real_data[gene_name]) == length(name_comb) * length(parsed_args["schedule"])
          if haskey(real_completed, chr)
            real_completed[chr][line_num] = true
          end
        end
      elseif sim == "SIM"
        if !haskey(sim_data, gene_name)
          sim_data[gene_name] = falses(length(name_comb) * length(parsed_args["schedule"]))
        end
        if data_type in name_comb
          idx1 = findfirst((x)-> x == data_type,name_comb)
          idx2 = findfirst((x)-> x == subsample,parsed_args["schedule"])
          sim_data[gene_name][(idx2-1)*length(name_comb) + idx1] = true
        else
          warn("unknown data type")
        end
        if sum(sim_data[gene_name]) == length(name_comb) * length(parsed_args["schedule"])
          if haskey(sim_completed, chr)
            sim_completed[chr][line_num] = true
          end
        end
      end
    end
  end
  out = STDOUT
  for (chr, lines) in real_completed
    for i in 1:length(lines)
      if !lines[i]
        println(out, "julia -p $(parsed_args["chains"]) $(parsed_args["temp"])/phase_by_loci_sub.jl -t $(parsed_args["in"])/$(parsed_args["prefix"])_true_$chr.txt -a $(parsed_args["in"])/$(parsed_args["prefix"])_reads_$chr.txt -o $(parsed_args["out"]) -l $i -r $chr -i $(parsed_args["iters"]) -b $(parsed_args["burnin"]) -c $(parsed_args["chains"]) -d $(parsed_args["temp"]) -n $(join(parsed_args["names"]," ")) -m $(join(vec(parsed_args["matrix"])," ")) -s $(join(parsed_args["schedule"]," "))")
      end
    end
  end
  for (chr, lines) in sim_completed
    for i in 1:length(lines)
      if !lines[i]
        #println(out,"bash $(parsed_args["dir"])/launch_mcmc.sh $(parsed_args["in"])/$(parsed_args["prefix"])_sim_$chr.txt $(parsed_args["in"])/$(parsed_args["prefix"])_reads_$chr.txt $(parsed_args["out"]) $i $chr u")
        println(out, "julia -p $(parsed_args["chains"]) $(parsed_args["temp"])/phase_by_loci_sub.jl -t $(parsed_args["in"])/$(parsed_args["prefix"])_sim_$chr.txt -a $(parsed_args["in"])/$(parsed_args["prefix"])_reads_$chr.txt -o $(parsed_args["out"]) -l $i -r $chr -i $(parsed_args["iters"]) -b $(parsed_args["burnin"]) -c $(parsed_args["chains"]) -d $(parsed_args["temp"]) -n $(join(parsed_args["names"]," ")) -m $(join(vec(parsed_args["matrix"])," ")) -s $(join(parsed_args["schedule"]," "))")
      end
    end
  end
  close(out)
    

  

end

main()
