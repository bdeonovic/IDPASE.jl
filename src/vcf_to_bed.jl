function vcf_to_bed(vcf_file_name::String, output_file_name::ASCIIString="vcf.bed")
  vcf = open(vcf_file_name,"r")
  for line in eachline(vcf)
    if ismatch(r"^#",line)
      continue
    end
    fields = split(strip(line),"\t")
    println("$(fields[1])\t$(parse(Int64,fields[2])-1)\t$(parse(Int64,fields[2]))\t"join(fields,"\t"))
  end
  close(vcf)
  return nothing
end

if length(ARGS)==1
  vcf_to_bed(ARGS[1])
end
