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
    "--base", "-b"
      arg_type = ASCIIString
    "--sim", "-s"
      action = :store_true
    "--ase", "-a"
      action = :store_true
    "--prefix", "-p"
      arg_type = ASCIIString
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
  #real_data = Dict{ASCIIString, Bool}()
  prefix = ""
  if parsed_args["ase"]
    prefix = "EXTRA"
  else
    prefix = "NONASE"
  end

  for file in readdir(parsed_args["in"])
    m = nothing
    if parsed_args["ase"]
      m = match(Regex("$(parsed_args["prefix"])_extra_SPECIAL_(\\w*).txt"),file)
    else
      m = match(Regex("$(parsed_args["prefix"])_non_ase_(\\w*).txt"),file)
    end
    if m != nothing
      num_loci = parse(Int64, split(readall(`wc -l $(parsed_args["in"])/$file`))[1])
      real_completed[m.captures[1]] = falses(num_loci)
    end
  end
  for file in readdir(parsed_args["out"])
    m = match(Regex("($prefix)::(\\w*)::(\\d*).txt"),file)
    if m != nothing
      sim = m.captures[1]
      chr = m.captures[2]
      line_num = parse(Int64, m.captures[3])
      if haskey(real_completed, chr)
        real_completed[chr][line_num] = true
      end
    end
  end

  out = STDOUT
  for (chr, lines) in real_completed
    for i in 1:length(lines)
      if !lines[i]
        if parsed_args["ase"]
          if parsed_args["sim"]
            println(out, "qsub -N EXTRA$(string(chr,"-",i)) $(parsed_args["base"])/launch_isoform_phase.sh $(parsed_args["in"])/$(parsed_args["prefix"])_extra_SPECIAL_$chr.txt $(parsed_args["out"])/EXTRA::$(chr)::$i.txt $i 1")
          else
            println(out, "qsub -N EXTRA$(string(chr,"-",i)) $(parsed_args["base"])/launch_isoform_phase.sh $(parsed_args["in"])/$(parsed_args["prefix"])_extra_SPECIAL_$chr.txt $(parsed_args["out"])/EXTRA::$(chr)::$i.txt $i 0")
          end
        else
          if parsed_args["sim"]
            #println(out, "qsub -N $prefix$(string(chr,"-",i)) $(parsed_args["base"])/launch_isoform_phase.sh $(parsed_args["in"])/out_non_ase_$chr.txt $(parsed_args["out"])/$prefix::$(chr)::$i.txt $i 1")
          else
            #println(out, "qsub -N $prefix$(string(chr,"-",i)) $(parsed_args["base"])/launch_isoform_phase.sh $(parsed_args["in"])/out_non_ase_$chr.txt $(parsed_args["out"])/$prefix::$(chr)::$i.txt $i 0")
          end
        end
      end
    end
  end
  close(out)
end

main()
