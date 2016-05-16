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
    "--line", "-l"
      arg_type = Int
      default = 1
      help = "loci line"
    "--sim", "-s"
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

  phase_isoform(parsed_args["in"], parsed_args["line"], parsed_args["out"], parsed_args["sim"])

end

main()
