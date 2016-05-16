using DataStructures
if length(ARGS) == 2
  isoform_level_input_file = ARGS[1]
  regions_input_file = ARGS[2]

  isoform_level_input = Dict{ASCIIString, ASCIIString}()
  regions_input = Dict{ASCIIString, ASCIIString}()

  for line in readlines(open(isoform_level_input_file, "r"))
    fields = split(strip(line),"\t")
    gene_name = fields[2]
    if !(gene_name in keys(isoform_level_input))
      isoform_level_input[gene_name] = line #[convert(ASCIIString, s) for s in split(fields[10],",")]
    end
  end
  for line in readlines(open(regions_input_file, "r"))
    fields = split(strip(line),"\t")
    gene_name = fields[2]
    if fields[5] == ""
      continue
    end
    if !(gene_name in keys(regions_input))
      regions_input[gene_name] = line
    end
  end
  for (gene_name, isoform_line) in isoform_level_input
    isoform_fields = split(strip(isoform_line),"\t")
    isoform_names = [convert(ASCIIString, s) for s in split(isoform_fields[10],",")]
    K = length(isoform_names)
    num_isoforms = parse(Int64, isoform_fields[3])
    num_regions = parse(Int64, isoform_fields[4])
    if num_regions ==0
      continue
    end
    regions_in_isoforms = reshape([parse(Int64, s) for s in split(isoform_fields[5],",")], num_regions, num_isoforms)
    isoform_map = [parse(Int64, s) for s in split(isoform_fields[7],",")]
    isoform_lengths = [parse(Int64,s) for s in split(isoform_fields[8],",")]
    region_lengths = [parse(Int64,s) for s in split(isoform_fields[9],",")]


    num_snps = 0
    snps_in_isoforms = zeros(Int64, num_snps, K)
    regions = Array{Int64,1}[]
    y = Int64[]
    #snp_idx = Int64[]
    if gene_name in keys(regions_input)
      fields = split(strip(regions_input[gene_name]), "\t")
      snps_in_isoforms_vec = [parse(Int64, s) for s in split(fields[6],",")]
      num_snps = round(Int64, length(snps_in_isoforms_vec) / K)
      snps_in_isoforms = reshape(snps_in_isoforms_vec, num_snps,K)
      regions = [ collect(eval(parse(string(strip(s,[')']),")"))))::Array{Int64,1} for s in split(fields[4],"),")]
      y = [parse(Int64, s) for s in split(fields[5],",")]
      #snp_idx = sort(unique(abs(vcat(regions...)[find(x-> x<0, vcat(regions...))])))
      snp_idx = [parse(Int64, s) for s in split(fields[8],",")]

    else
      warn("$gene_name not in regions_input")
      continue
    end
    new_regions = OrderedDict{Set{Int64},Int64}()
    new_y = Int64[]
    new_lengths = Int64[]
    new_regions_in_isoforms = zeros(Int64, length(new_y), K)
    for i in 1:length(regions)
      new_set = setdiff(abs(regions[i]),snp_idx)
      if Set(new_set) in keys(new_regions)
        idx = new_regions[Set(new_set)]
        new_y[idx] += y[i]
        new_regions_in_isoforms[idx,:] += regions_in_isoforms[i,1:K]
        new_lengths[idx] = max(new_lengths[idx], region_lengths[i])
      else
        new_regions[Set(new_set)] = length(new_y) + 1
        push!(new_y,y[i])
        new_regions_in_isoforms = vcat(new_regions_in_isoforms, regions_in_isoforms[i,1:K])
        push!(new_lengths,region_lengths[i])
      end
    end
    new_regions_in_isoforms = 1*(new_regions_in_isoforms .> 0)

    println(string(isoform_fields[1],"\t",gene_name,"\t", K, "\t", length(new_y), "\t", join(vec(new_regions_in_isoforms),","),"\t", join(new_y,","),"\t",join(collect(1:K),","),"\t",join(isoform_lengths[1:K],","),"\t", join(new_lengths,","),"\t", join(isoform_names,",")))
  end

end
