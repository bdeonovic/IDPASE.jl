if length(ARGS) == 1
  IN = open(ARGS[1],"r")
  t = 1
  for line in eachline(IN)
    fields = split(strip(line),"\t")
    out_fields = Array(ASCIIString,16)
    for i in 1:length(fields)
      out_fields[i] = fields[i]
    end
    num_blocks= parse(Int64,fields[9])
    out_fields[12] = "0"
    out_fields[13] = out_fields[1]
    out_fields[14] = "unk"
    out_fields[15] = "unk"
    out_fields[16] = join([fill("-1",num_blocks);""],",")
    out_fields[1] = "$t"
    println(join(out_fields,"\t"))
    t += 1
  end
  close(IN)
end
