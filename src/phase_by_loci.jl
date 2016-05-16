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
    "--in", "-i"
      arg_type = ASCIIString
    "--out", "-o"
      arg_type = ASCIIString
    "--types", "-t"
      arg_type = ASCIIString
      nargs = '+'
      required = true
    "--base", "-b"
      arg_type = ASCIIString
    "--qstat", "-q"
       action = :store_true
  end

  return parse_args(s)
end

function main()
  parsed_args = parse_commandline()
  println(STDERR, "Parsed args:")
  for (arg,val) in parsed_args
    println(STDERR, "  $arg  =>  $val")
  end

  real_completed = Dict{ASCIIString, Array{Bool,1}}()
  real_data = Dict{ASCIIString, Array{Bool,1}}()

  sim_completed = Dict{ASCIIString, Array{Bool,1}}()
  sim_data = Dict{ASCIIString, Array{Bool,1}}()
  for file in readdir(parsed_args["in"])
    m = match(r"out_(true|sim)_(\w*).txt",file)
    if m != nothing
      num_loci = parse(Int64, split(readall(`wc -l $(parsed_args["in"])/$file`))[1])
      if m.captures[1] == "true"
        real_completed[m.captures[2]] = falses(num_loci)
      elseif m.captures[1] == "sim"
        sim_completed[m.captures[2]] = falses(num_loci)
      end
    end
  end
  for file in readdir(parsed_args["out"])
    m = match(Regex("(REAL|SIM)::(\\w*)::(\\d*)::([\\w\\[\\]\\/\\.\\-]*)::($(join(parsed_args["types"],"|"))).txt"),file)
    if m != nothing
      sim = m.captures[1]
      chr = m.captures[2]
      line_num = parse(Int64, m.captures[3])
      gene_name = m.captures[4]
      data_type = m.captures[5]
      if sim == "REAL"
        if !haskey(real_data, gene_name)
          real_data[gene_name] = falses(length(parsed_args["types"]))
        end
        if data_type in parsed_args["types"]
          idx = findfirst((x)-> x == data_type,parsed_args["types"])
          real_data[gene_name][idx] = true
        else
          warn("unknown data type")
        end
        if sum(real_data[gene_name]) == length(parsed_args["types"])
          if haskey(real_completed, chr)
            real_completed[chr][line_num] = true
          end
        end
      elseif sim == "SIM"
        if !haskey(sim_data, gene_name)
          sim_data[gene_name] = falses(length(parsed_args["types"]))
        end
        if data_type in parsed_args["types"]
          idx = findfirst((x)-> x == data_type,parsed_args["types"])
          sim_data[gene_name][idx] = true
        else
          warn("unknown data type")
        end
        if sum(sim_data[gene_name]) == length(parsed_args["types"])
          if haskey(sim_completed, chr)
            sim_completed[chr][line_num] = true
          end
        end
      end
    end
  end
  if parsed_args["qstat"]
    running = readlines(`qstat -u bdeonovic`)
    for i in 3:length(running)
      fields = split(running[i])
      m = match(r"(REAL|SIM)(\w*)-(\d*)", fields[3])
      if m != nothing
        sim = m.captures[1]
        chr = m.captures[2]
        line = parse(Int64, m.captures[3])
        if sim == "REAL"
          real_completed[chr][line] = true
        elseif sim == "SIM"
          sim_completed[chr][line] = true
        else
          warn("problem")
        end
      end
    end
  end

  #out = open(string(parsed_args["base"],"/to_run.sh"),"w")
  out = STDOUT
  for (chr, lines) in real_completed
    for i in 1:length(lines)
      if !lines[i]
        println(out, "qsub -N REAL$(string(chr,"-",i)) /Shared/Au/bdeonovic/work/haplotype/gold_standard7/launch_mcmc.sh $(parsed_args["in"])/out_true_$chr.txt $(parsed_args["in"])/out_reads_$chr.txt $(parsed_args["out"]) $i $chr x")
      end
    end
  end
  for (chr, lines) in sim_completed
    for i in 1:length(lines)
      if !lines[i]
        println(out,"qsub -N SIM$(string(chr,"-",i)) /Shared/Au/bdeonovic/work/haplotype/gold_standard7/launch_mcmc.sh $(parsed_args["in"])/out_sim_$chr.txt $(parsed_args["in"])/out_reads_$chr.txt $(parsed_args["out"]) $i $chr u")
      end
    end
  end
  close(out)

  

end

main()
