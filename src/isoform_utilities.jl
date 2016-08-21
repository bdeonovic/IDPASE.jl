function get_fpkm(fpkm_file_name::AbstractString)
  fpkm_file = open(fpkm_file_name,"r")
  fpkm_dict = Dict{ASCIIString, Array{Float64,1}}()

  for line in eachline(fpkm_file)
    fields = split(strip(line),"\t")
    isoform_name = fields[1]
    isoform_fpkm = parse(Float64, fields[2])
    gene_fpkm = parse(Float64, fields[3])
    if !(isoform_name in keys(fpkm_dict))
      fpkm_dict[isoform_name] = [isoform_fpkm, gene_fpkm]
    end
  end
  close(fpkm_file)
  return fpkm_dict
end

function init_loci(gpd_file_name::AbstractString, chr::AbstractString, num_src_types::Integer, fpkm_dict::Dict{ASCIIString, Array{Float64, 1}})
  gpd_file = open(gpd_file_name,"r")
  isoform_dict = Dict{ASCIIString,Array{GPDEntry{ASCIIString,Int64},1}}()
  loci_dict = Dict{ASCIIString,  IsoformPhaseEntry{ASCIIString, Int64}}()
  transcripts = Dict{ASCIIString,Bool}()

  for line in eachline(gpd_file)
    fields = split(strip(line),"\t")
    block_ends   = [parse(Int64,s) for s in split(strip(fields[14],','),",")]
    block_starts = [parse(Int64,s) for s in split(strip(fields[13],','),",")]
    frames = [parse(Int64,s) for s in split(strip(fields[19],','),",")]
    isoform_length = sum(block_ends - block_starts)

    if length(keys(fpkm_dict)) > 0
      if fields[5] in keys(fpkm_dict)
        fpkm = fpkm_dict[fields[5]]
        if fpkm[1] < 0.01 * fpkm[2]
          continue
        end
      else
        continue
      end
    end

    if !haskey(transcripts, fields[5])
      isoform = GPDEntry{ASCIIString,Int64}(parse(Int64,fields[4]),fields[5],fields[6],fields[7],[parse(Int64,s) for s in fields[8:12]]...,block_starts,block_ends,parse(Int64,fields[15]),fields[16],fields[17],fields[18],frames)
      if haskey(isoform_dict,fields[16])
        if !(fields[2] in [isoform_dict[fields[16]][j].name for j in 1:length(isoform_dict[fields[16]])])
          push!(isoform_dict[fields[16]],isoform)
        end
      else
        isoform_dict[fields[16]]   = [isoform]
        loci = IsoformPhaseEntry{ASCIIString, Int64}(chr,Int64[],ASCIIString[],ASCIIString[],Dict{ASCIIString,Bool}(),Dict{ASCIIString,ASCIIString}(),Int64[],Int64[], NaN, true,PSLEntry{ASCIIString,Int64}[],ASCIIString[], ASCIIString[], fill(2,0,0),fill(1.0,0,0),zeros(Int64,num_src_types), fill(0,0,0), fill(0,0,0), fill(0,0,0), fill(0,0,0), Int64[], Int64[], DataStructures.OrderedDict{Set{Int64},Int64}(), Int64[], Int64[], Int64[] )
        loci_dict[fields[16]]  = loci
      end
      transcripts[fields[5]] = true
    end
  end
  close(gpd_file)
  return isoform_dict, loci_dict, transcripts
end

function get_snps!{V<:AbstractString,W<:Integer}(loci_dict::Dict{V,IsoformPhaseEntry{V,W}}, gpd_file_name::V, vcf_file_name::V, temp_dir::V, chr::V, fpkm_dict::Dict{ASCIIString, Array{Float64,1}})
  gpd_file = open(gpd_file_name,"r")
  gpd_with_snps = open("$(temp_dir)/$(chr)_gpd_with_snps.bed","w")

  for line in eachline(open(`bedtools intersect -wo -a $gpd_file_name -b $vcf_file_name`)[1])
    fields = split(strip(line),"\t")
    block_ends   = [parse(Int64,s) for s in split(strip(fields[14],','),",")]
    block_starts = [parse(Int64,s) for s in split(strip(fields[13],','),",")]
    frames = [parse(Int64,s) for s in split(strip(fields[19],','),",")]
    isoform_length = sum(block_ends - block_starts)

    if length(keys(fpkm_dict)) > 0 
      if fields[5] in keys(fpkm_dict)
        fpkm = fpkm_dict[fields[5]]
        if fpkm[1] < 0.01 * fpkm[2]
          continue
        end
      else
        continue
      end
    end

    gpd = GPDEntry{ASCIIString,Int64}(parse(Int64,fields[4]),fields[5],fields[6],fields[7],[parse(Int64,s) for s in fields[8:12]]...,block_starts,block_ends,parse(Int64,fields[15]),fields[16],fields[17],fields[18],frames)
    print_gpd(gpd_with_snps,gpd)

    if haskey(loci_dict,fields[16])
      loci = loci_dict[fields[16]]
      snp_name = fields[25] #50
      if snp_name == "."
        snp_name = string(fields[23],":",fields[24])
      end
      if !haskey(loci.snp_names,snp_name)
        loci.snp_names[snp_name] = true
        push!(loci.ref,fields[26])
        push!(loci.alt,fields[27])
        push!(loci.snps,parse(Int64,fields[24]))
        genotype_fields = split(fields[31],":")
        GT_idx = find((x)-> x=="GT",genotype_fields)
        genotype_values = split(fields[32],":")
        if length(GT_idx) != 1
          continue
        end
        if ismatch(r"\|",genotype_values[GT_idx[1]])
          push!(loci.true_h,parse(Int64,split(genotype_values[GT_idx[1]],"|")[1]))
        elseif ismatch(r"\/",genotype_values[GT_idx[1]])
          loci.phased = false
          push!(loci.true_h,parse(Int64,split(genotype_values[GT_idx[1]],"/")[1]))
        else
          error("bad VCF")
          continue
        end
        #sort snps
        idx = sortperm(loci.snps)
        loci.ref = loci.ref[idx]
        loci.alt = loci.alt[idx]
        loci.snps = loci.snps[idx]
        loci.true_h = loci.true_h[idx]
      end
      #loci_dict[fields[16]] = loci
    end
  end
  close(gpd_file)
  close(gpd_with_snps)
end

function assemble_loci!{V<:AbstractString,W<:Integer}(loci_dict::Dict{V,IsoformPhaseEntry{V,W}}, isoform_dict::Dict{V, Array{GPDEntry{V, W},1}};
  read_length::Integer = 100, min_overlap_length::Integer=10, verbose::Bool = false, mean_insert_length::Integer=100)
  for (loci_name, loci) in loci_dict
    if verbose
      println(loci_name)
    end
    isoforms = haskey(loci_dict,loci_name) ? isoform_dict[loci_name] : continue
    snps = loci.snps
    if length(snps) == 0 || length(loci.est_h) == 0
      delete!(loci_dict, loci_name)
      delete!(isoform_dict, loci_name)
      continue
    end

    regions, regions_in_isoforms, effective_region_lengths, effective_isoform_lengths, exons, snp_idx, exons_in_isoforms, isoforms_with_alt, snps_in_isoforms = get_regions(isoforms, snps, loci.est_h, read_length, min_overlap_length, mean_insert_length)

    num_exons, num_reg_isoforms = size(exons_in_isoforms)
    for i in 1:num_reg_isoforms
      push!(loci.isoform_names, isoforms[i].name)
    end

    num_regions = length(regions)
    num_isoforms = length(effective_isoform_lengths)
    isoform_map = [1:num_reg_isoforms; find(isoforms_with_alt)]

    loci.exons = exons
    loci.exons_in_isoforms = exons_in_isoforms
    loci.snps_in_isoforms = snps_in_isoforms
    loci.regions_in_isoforms = regions_in_isoforms
    loci.effective_region_lengths = effective_region_lengths
    loci.effective_isoform_lengths = effective_isoform_lengths
    loci.regions = regions
    loci.snp_idx = snp_idx
    loci.isoform_map = isoform_map
  end
end
function consecutive_groups(x::Vector{Int})
  result = Array{Int, 1}[]
  for set_size in 0:length(x)
    for i in 1:(length(x) - set_size)
      push!(result, x[i:(i + set_size)])
    end
  end
  return result
end
function get_regions{U<:AbstractString,V<:Integer}(isoforms::Vector{GPDEntry{U,V}}, snps::Vector{V}, haplotype::Vector{V}, read_length::Int, min_overlap_length::Int, mean_insert_length::Int, verbose::Bool=false)
  num_real_isoforms = length(isoforms)
  starts = Int64[]
  ends = Int64[]
  for i in 1:num_real_isoforms
    append!(starts, isoforms[i].exon_starts)
    append!(ends, isoforms[i].exon_ends)
  end

  pnts = sort(unique([starts; ends]))
  positions = union([starts[i]:ends[i] for i in 1:length(starts)]...)
  exons = zeros(Int64,0,2)
  for i in 1:(length(pnts)-1)
    inter = sort(intersect(positions, pnts[i]:pnts[i+1]))
    if length(inter) >0 && inter != [pnts[i]; pnts[i+1]]
      exons = vcat(exons, [pnts[i]; pnts[i+1]]')
    end
  end
  num_exons = size(exons,1)
  num_snps = length(snps)
  exons = exons[sortperm(exons[:,1]),:]

  snps_in_exons = zeros(Int64, num_exons, num_snps)
  for i in 1:num_snps
    for j in 1:num_exons
      if snps[i] >= exons[j,1] && snps[i] <= exons[j,2]
        snps_in_exons[j,i] = 1
      end
    end
  end

  exons_in_isoforms = zeros(Int64, num_exons, num_real_isoforms)
  snps_in_isoforms = zeros(Int64, num_snps, num_real_isoforms)
  isoforms_with_alt = zeros(Int64, num_real_isoforms)
  for j in 1:num_real_isoforms
    for i in 1:num_exons
      if  any((exons[i,1] .>= isoforms[j].exon_starts) & (exons[i,2] .<= isoforms[j].exon_ends))
        exons_in_isoforms[i,j] = 1
      end
    end
    for i in 1:num_snps
      if any( (snps[i] .>= isoforms[j].exon_starts) & (snps[i] .<= isoforms[j].exon_ends))
        snps_in_isoforms[i,j] = 1
        isoforms_with_alt[j] = 1
      end
    end
  end
  snp_idx = num_exons + collect(1:length(snps))
  exon_lengths = exons[:,2] - exons[:,1]

  regions = DataStructures.OrderedDict{Set{Int},Int}()
  num_isoforms = sum(2 * isoforms_with_alt + (1 - isoforms_with_alt))
  regions_in_isoforms = zeros(Int64, 0, num_isoforms)
  effective_region_lengths = Int64[]
  non_junction_regions = Int64[]
  z = 1
  for j in 1:num_real_isoforms
    for curr_set in consecutive_groups(find(exons_in_isoforms[:,j]))
      if !haskey(regions, Set(curr_set))
        l = 0
        l1=l2=l3=0
        if verbose
          println(STDERR, curr_set)
        end
        if length(curr_set) == 1
          id = find(mapslices(sum,snps_in_exons[curr_set,:],1))
          if length(id) > 0
            lox = [exons[curr_set[1],1]; snps[id]; exons[curr_set[1],2]]
            for sn in 2:length(lox)
              l += max(lox[sn] - lox[sn-1] - read_length, 0)
            end
          else
            l = exon_lengths[curr_set[1]] - read_length + 1
          end
        else
          id = find(mapslices(sum,snps_in_exons[curr_set,:],1))
          if length(id) > 0
            if sum(snps_in_exons[curr_set[2:(end-1)],:]) > 0
              l = 0
            elseif sum(snps_in_exons[curr_set[1],:]) > 0
              if sum(snps_in_exons[curr_set[end],:]) > 0
                l1 = exons[curr_set[1],2] - snps[find(mapslices(sum,snps_in_exons[curr_set[1],:],1))[end]]
                l2 = snps[find(mapslices(sum,snps_in_exons[curr_set[end],:],1))[1]] - exons[curr_set[end],1]
                l3 = read_length - sum(exon_lengths[curr_set[2:(end-1)]]) - min_overlap_length
                if l3 + min_overlap_length > l1 + l2
                  l = 0
                else
                  l = min(l1,l2,l3)
                end
              else
                l1 = exons[curr_set[1],2] - snps[find(mapslices(sum,snps_in_exons[curr_set[1],:],1))[end]]
                l2 = exon_lengths[curr_set[end]]
                l3 = read_length - sum(exon_lengths[curr_set[2:(end-1)]]) - min_overlap_length
                if l3 + min_overlap_length > l1 + l2
                  l = 0
                else
                  l = min(l1,l2,l3)
                end
              end
            elseif sum(snps_in_exons[curr_set[end],:]) > 0
              l1 = exon_lengths[curr_set[1]]
              l2 = snps[find(mapslices(sum,snps_in_exons[curr_set[end],:],1))[1]] - exons[curr_set[end],1]
              l3 = read_length - sum(exon_lengths[curr_set[2:(end-1)]]) - min_overlap_length
              if l3 + min_overlap_length > l1 + l2
                l = 0
              else
                l = min(l1,l2,l3)
              end
            end
          else
            l1 = exon_lengths[curr_set[1]]
            l2 = exon_lengths[curr_set[end]]
            l3 = read_length - sum(exon_lengths[curr_set[2:(end-1)]]) - min_overlap_length
            if l3 + min_overlap_length > l1 + l2
              l = 0
            else
              l = min(l1,l2,l3)
            end
          end
        end
        if verbose
          println(STDERR,"$l, $l1, $l2, $l3")
        end
        if l > 0 
          regions[Set(curr_set)] = z
          z += 1
          push!(effective_region_lengths, l)

          if length(curr_set) == 1
            push!(non_junction_regions, 1)
          else
            push!(non_junction_regions, 0)
          end

          line = zeros(Int64, num_isoforms)
          line[j] = 1
          if isoforms_with_alt[j] == 1
            line[num_real_isoforms + sum(isoforms_with_alt[1:j])] = 1 
          end
          regions_in_isoforms = vcat(regions_in_isoforms, line')
        end
      else
        row_i = regions[Set(curr_set)]
        regions_in_isoforms[row_i, j] = 1
        if isoforms_with_alt[j] == 1
          regions_in_isoforms[row_i, num_real_isoforms + sum(isoforms_with_alt[1:j])] = 1 
        end
      end
      for h in [1,0]
        id = find(mapslices(sum,snps_in_exons[curr_set,:],1))
        if length(id) > 0
          cand_snps = -snp_idx[id] .* (2 * (h*haplotype[id] + (1-h)*(1-haplotype[id])) - 1)
          for s in consecutive_groups(cand_snps)
            sugg_set = [curr_set; s]
            snp_exons = find(mapslices(sum,snps_in_exons[:, (abs(s) - minimum(snp_idx) +1)],2))
            if !haskey(regions, Set(sugg_set))
              l = 0
              l_s = 0
              l_e = 0
              r_s = 0
              r_e = 0
              interior = 0
              if length(curr_set) == 1
                l_s = exons[curr_set[1],1]
                r_e = exons[curr_set[1],2]

                l_e = snps[abs(s)[1] - minimum(snp_idx)+1]
                r_s = snps[abs(s)[end] - minimum(snp_idx)+1]

                interior = r_s - l_e

                for sn in 1:num_snps
                  if (snps_in_exons[curr_set[1],sn] == 1) && snps[sn] < snps[abs(s)[1]-minimum(snp_idx)+1]
                    l_s = snps[sn]
                  end
                  if (snps_in_exons[curr_set[1],sn] == 1) && snps[sn] > snps[abs(s)[end]-minimum(snp_idx)+1]
                    r_e = snps[sn]
                  end
                end
                l1 = l_e - l_s
                l2 = r_e - r_s
                l3 = read_length - interior - min_overlap_length
                if l3 + min_overlap_length > l1 + l2
                  l = 0
                else
                  l = min(l1,l2,l3)
                end
              else
                interior = sum(exon_lengths[curr_set[2:(end-1)]])

                #get l_e and r_s
                id1 = find(mapslices(sum,snps_in_exons[curr_set[1],abs(s) - minimum(snp_idx)+1],1))
                if length(id1) > 0
                  l_e = snps[(abs(s) - minimum(snp_idx)+1)[id1[1]]]
                  interior += exons[curr_set[1],2] - l_e
                else
                  l_e = exons[curr_set[1],2]
                end

                id2 = find(mapslices(sum,snps_in_exons[curr_set[end],abs(s) - minimum(snp_idx)+1],1))
                if length(id2) > 0
                  r_s = snps[(abs(s) - minimum(snp_idx)+1)[id2[end]]]
                  interior += r_s - exons[curr_set[end],1]
                else
                  r_s = exons[curr_set[end],1]
                end

                #get l_s and r_e
                id1 = Int64[]
                snps_not_covered = setdiff(snp_idx,abs(s))
                if length(snps_not_covered) > 0
                  id1 = find(mapslices(sum,snps_in_exons[curr_set[1],snps_not_covered - minimum(snp_idx) + 1],1))
                end
                if length(id1) > 0
                  l_s = snps[(snps_not_covered - minimum(snp_idx)+1)[id1[end]]]
                else
                  l_s = exons[curr_set[1],1]
                end

                id2 = Int64[]
                if length(snps_not_covered) > 0
                  id2 = find(mapslices(sum,snps_in_exons[curr_set[end],snps_not_covered - minimum(snp_idx) + 1],1))
                end
                if length(id2) > 0
                  r_e = snps[(snps_not_covered - minimum(snp_idx)+1)[id2[1]]]
                else
                  r_e = exons[curr_set[end],2]
                end
                l1 = l_e-l_s
                l2 = r_e-r_s
                l3 = read_length - interior - min_overlap_length
                if l3 + min_overlap_length > l1 + l2
                  l = 0
                else
                  l = min(l1,l2,l3)
                end

              end
              if l > 0 
                regions[Set(sugg_set)] = z
                z += 1
                push!(effective_region_lengths, l)

                push!(non_junction_regions, 0)

                line = zeros(Int64, num_isoforms)
                if h == 1
                  line[j] = 1
                else
                  line[num_real_isoforms + sum(isoforms_with_alt[1:j])] = 1
                end
                regions_in_isoforms = vcat(regions_in_isoforms, line')
              end
            else
              row_i = regions[Set(sugg_set)]
              if h == 1
                regions_in_isoforms[row_i, j] = 1
              else
                regions_in_isoforms[row_i, num_real_isoforms + sum(isoforms_with_alt[1:j])] = 1
              end
            end
          end
        end
      end
    end
  end
  effective_isoform_lengths = transpose(regions_in_isoforms) * (effective_region_lengths .* non_junction_regions)
  return regions, regions_in_isoforms, effective_region_lengths, effective_isoform_lengths, exons, snp_idx, exons_in_isoforms, isoforms_with_alt, snps_in_isoforms
end

function get_gene_level_results{V<:AbstractString,W<:Integer}(loci_dict::Dict{V,IsoformPhaseEntry{V,W}}, mcmc_result_file::AbstractString, result_type::AbstractString, use_true::Bool=false)
  mcmc_results_vec = readlines(open(mcmc_result_file, "r"))
  num_mcmc_results = length(mcmc_results_vec)
  mcmc_results = Array(Any, length(mcmc_results_vec), 28)
  for i in 1:num_mcmc_results
    fields = split(mcmc_results_vec[i])
    mcmc_results[i,1] = fields[1]
    mcmc_results[i,2] = fields[2]
    mcmc_results[i,3] = parse(Int64, fields[3])
    mcmc_results[i,4] = parse(Int64, fields[4])
    mcmc_results[i,5] = fields[5]
    mcmc_results[i,6] = fields[6]
    mcmc_results[i,7] = fields[7]
    mcmc_results[i,8] = fields[8]
    mcmc_results[i,9] = fields[9]
    mcmc_results[i,10] = fields[10] == "NA" ? NaN : parse(Float64, fields[10]) #cond mean
    mcmc_results[i,11] = fields[11] == "NA" ? NaN : parse(Float64, fields[11]) #mode
    mcmc_results[i,12] = fields[12] == "NA" ? NaN : parse(Float64, fields[12]) #median
    mcmc_results[i,13] = fields[13] == "NA" ? NaN : parse(Float64, fields[13]) #mean
    mcmc_results[i,14] = fields[14] == "NA" ? NaN : parse(Float64, fields[14])
    mcmc_results[i,15] = fields[15] == "NA" ? NaN : parse(Float64, fields[15])
    mcmc_results[i,16] = fields[16] == "NA" ? NaN : parse(Float64, fields[16])
    mcmc_results[i,17] = fields[17] == "NA" ? NaN : parse(Float64, fields[17])
    mcmc_results[i,18] = fields[18] == "NA" ? NaN : parse(Float64, fields[18])
    mcmc_results[i,19] = fields[19] == "NA" ? NaN : parse(Float64, fields[19])
    mcmc_results[i,20] = fields[20] == "NA" ? NaN : parse(Float64, fields[20])
    mcmc_results[i,21] = fields[21] == "NA" ? NaN : parse(Float64, fields[21])
    mcmc_results[i,22] = fields[22] == "NA" ? NaN : parse(Float64, fields[22])
    mcmc_results[i,23] = fields[23] == "NA" ? NaN : parse(Int64, fields[23])
    mcmc_results[i,24] = fields[24] == "NA" ? NaN : parse(Int64, fields[24])
    mcmc_results[i,25] = fields[25]
    mcmc_results[i,26] = fields[26]
    mcmc_results[i,27] = parse(Float64, fields[27])
    mcmc_results[i,28] = fields[28]

    num_sr = [parse(Int64, s) for s in split(mcmc_results[i,2],",")][1]

    if mcmc_results[i,5] in keys(loci_dict)
      if (mcmc_results[i,26] != result_type) || (mcmc_results[i,27] != 1.0) #"SGSSNY_1.0"
        continue
      end
      if (mcmc_results[i,8] == "NA") || (length(loci_dict[mcmc_results[i,5]].snps) != mcmc_results[i,4]) || num_sr > 10000
        println("deleting $(mcmc_results[i,5])")
        delete!(loci_dict, mcmc_results[i,5])
        continue
      end
      snp_coords = [parse(Int64, s) for s in split(mcmc_results[i,25],",")] #gene level
      new_snps = loci_dict[mcmc_results[i,5]].snps
      idx1 = findin(new_snps,snp_coords)

      idx2 = findin(snp_coords,new_snps)
      if length(idx1) == 0 || length(idx2) == 0
        println("deleting $(mcmc_results[i,5])")
        delete!(loci_dict, mcmc_results[i,5])
        continue
      end

      loci_dict[mcmc_results[i,5]].snps = loci_dict[mcmc_results[i,5]].snps[idx1]
      loci_dict[mcmc_results[i,5]].ref = loci_dict[mcmc_results[i,5]].ref[idx1]
      loci_dict[mcmc_results[i,5]].alt = loci_dict[mcmc_results[i,5]].alt[idx1]
      if use_true
        loci_dict[mcmc_results[i,5]].est_h = [parse(Int64, s) for s in split(mcmc_results[i,7],"")][idx2]
        if loci_dict[mcmc_results[i,5]].est_h[1] == 1
          loci_dict[mcmc_results[i,5]].rho_mode = 1 - mcmc_results[i,13]
        else
          loci_dict[mcmc_results[i,5]].rho_mode = mcmc_results[i,13]
        end
      else
        loci_dict[mcmc_results[i,5]].est_h = [parse(Int64, s) for s in split(mcmc_results[i,8],"")][idx2]
        loci_dict[mcmc_results[i,5]].rho_mode = mcmc_results[i,13]
      end
    end
  end
end

function bin_reads{V<:AbstractString,W<:Integer}(loci_dict::Dict{V,IsoformPhaseEntry{V,W}}, 
  isoform_dict::Dict{V, Array{GPDEntry{V, W},1}}, 
  prefix::AbstractString,temp_dir::AbstractString; stitched::Bool=true, 
  verbose::Bool=false, verbose2::Bool=false, verbose3::Bool=false, use_true::Bool=false)
  seen_dict = Set{ASCIIString}()

  gpd_file_name = "$(temp_dir)/$(prefix)_gpd.bed"
  psl_file_name = "$(temp_dir)/$(prefix)_psl_phase.bed"
  if !stitched
    psl_file_name = "$(temp_dir)/$(prefix)_psl1.bed"
  end
  println(STDERR, "using psl file $psl_file_name")
  verboseOUT = open("verbose_out","w")
  for line in eachline(open(`bedtools intersect -wo -a $psl_file_name -b $gpd_file_name`)[1])
    fields = split(strip(line),"\t")
    block_sizes = [parse(Int64,s) for s in split(strip(fields[22],','),",")]
    q_starts = [parse(Int64,s) for s in split(strip(fields[23],','),",")]
    t_starts = [parse(Int64,s) for s in split(strip(fields[24],','),",")]
    num_blocks = length(block_sizes)

    read_name = fields[13]
    gene_name = fields[40]


    read = PSLEntry{ASCIIString,Int64}([parse(Int64,s) for s in fields[4:11]]...,fields[12],fields[13],[parse(Int64,s) for s in fields[14:16]]...,fields[17],[parse(Int64,s) for s in fields[18:21]]...,block_sizes,q_starts,t_starts)

    if haskey(loci_dict, gene_name)
      loci = loci_dict[gene_name]

      #for verbose3
      regions_on_transcript = hcat([1; cumsum(mapslices(diff,loci.exons,2))[1:(end-1)] + 1], cumsum(mapslices(diff,loci.exons,2)))

      num_regions = length(loci.effective_region_lengths)
      in_exons = Int64[]
      in_snps = Int64[]
      num_exons = size(loci.exons,1)
      num_snps = length(loci.snps)
      for k in 1:num_blocks
        for j in 1:num_exons
          if (loci.exons[j,1] < t_starts[k] < loci.exons[j,2]) ||
             (loci.exons[j,1] < t_starts[k] + block_sizes[k] < loci.exons[j,2]) ||
             ((loci.exons[j,1] >= t_starts[k]) && (loci.exons[j,2] <= t_starts[k] + block_sizes[k]))
             push!(in_exons,j)
           end
        end
        for j in 1:num_snps
          if (t_starts[k] + block_sizes[k] * .1) <= loci.snps[j] <= (t_starts[k] + block_sizes[k] * .9)
            push!(in_snps,j)
          end
        end
      end
      in_exons = unique(in_exons)
      in_snps = unique(in_snps)


      if length(loci.Y) == 0
        loci.Y = zeros(Int64, num_regions)
      end
      if length(in_snps) > 0 
        continue
      end
      if verbose
        println(STDERR, "in exons $in_exons")
      end
      if haskey(loci.regions, Set(in_exons))
        if read_name in seen_dict
          continue
        else
          push!(seen_dict, read_name)
        end
        loci.Y[loci.regions[Set(in_exons)]] += 1
        if verbose3
          println(verboseOUT, string(read.q_name, " maps to regions $in_exons which have reference coords ", join([join(loci.exons[i,:],",") for i in in_exons],"\t")," and query coords ", join([join(regions_on_transcript[i,:],",") for i in in_exons],"\t")," t_starts: ",join(t_starts,","),"t_ends: ",join(t_starts+block_sizes,",")))
        end
      elseif verbose || verbose2 || verbose3
        println(STDERR, "regions does not contain $in_exons, from read $read_name")
        println(verboseOUT, "regions does not contain $in_exons, from read $read_name")
      end
    end
  end
  for (gene_name, loci) in loci_dict
    for read_name in loci.read_names_order
      if read_name in seen_dict
        continue
      end
      if loci.read_counts[1] <= 0
        continue
      end
      read_idx  = findfirst(read_name .== [loci.reads[i].q_name for i in 1:loci.read_counts[1]])
      if read_idx == 0
        continue
      end
      read = loci.reads[read_idx]
      block_sizes = read.block_sizes
      q_starts = read.q_starts
      t_starts = read.t_starts
      num_blocks = length(block_sizes)

      regions_on_transcript = hcat([1; cumsum(mapslices(diff,loci.exons,2))[1:(end-1)] + 1], cumsum(mapslices(diff,loci.exons,2)))
      num_regions = length(loci.effective_region_lengths)
      num_exons, num_real_isoforms = size(loci.exons_in_isoforms)
      in_exons = Int64[]
      in_snps = Int64[]
      num_exons = size(loci.exons,1)
      num_snps = length(loci.snps)
      for k in 1:num_blocks
        for j in 1:num_exons
          if (loci.exons[j,1] < t_starts[k] < loci.exons[j,2]) ||
             (loci.exons[j,1] < t_starts[k] + block_sizes[k] < loci.exons[j,2]) ||
             ((loci.exons[j,1] >= t_starts[k]) && (loci.exons[j,2] <= t_starts[k] + block_sizes[k]))
             push!(in_exons,j)
           end
        end
        for j in 1:num_snps
          if (t_starts[k] + block_sizes[k] * .1) <= loci.snps[j] <= (t_starts[k] + block_sizes[k] * .9)
            push!(in_snps,j)
          end
        end
      end
      in_exons = unique(in_exons)
      in_snps = unique(in_snps)
      if verbose3
        println(verboseOUT, string(read.q_name, " maps to regions $in_exons and snps $in_snps which have reference coords ", join([join(loci.exons[i,:],",") for i in in_exons],"\t")," and query coords ", join([join(regions_on_transcript[i,:],",") for i in in_exons],"\t")," t_starts: ",join(t_starts,","),"t_ends: ",join(t_starts+block_sizes,",")))
      end
      if length(loci.Y) == 0
        loci.Y = zeros(Int64, num_regions)
      end
      if length(in_snps) > 0 
        in_snps = in_snps + num_exons
        in_exons = [in_exons; in_snps]
        k = findfirst((x)-> x == read_name, loci.read_names_order)
        if verbose
          println(STDERR, "$read_name: Found a read with a snp! It is read number $k")
        end
        if k == 0
          continue
        else
          #prob_h = haplo_lik(vec(loci.X[k,:]), loci.est_h, vec(loci.Q[k,:]), loci.rho_mode)
          #if verbose
          #  println(STDERR, "X is $(loci.X[k,:]) Q is $(loci.Q[k,:]) resulting in prob_h of $prob_h")
          #end
          diff1 = 0
          diff2 = 0
          idx = find((loci.X[k,:] .!= 3) & (loci.X[k,:] .!= 2))
          for j in idx
            if loci.X[k,j] != loci.est_h[j]
              diff1 += 1
            else
              diff2 += 1
            end
          end
          
          haplotype = loci.est_h
          #if rand() >= prob_h
          #if prob_h < 0.5
          if diff1 == diff2
            continue
          elseif diff1 > diff2
            haplotype = 1 - loci.est_h
          end
          for s in 1:length(haplotype)
            if haplotype[s] == 1
              idx = findfirst( (x) -> x == loci.snp_idx[s], in_exons)
              if idx == 0
                continue
              end
              in_exons[idx] = -in_exons[idx]
            end
          end
        end
      end
      if verbose
        println(STDERR, "in exons $in_exons")
      end
      if haskey(loci.regions, Set(in_exons))
        loci.Y[loci.regions[Set(in_exons)]] += 1
      elseif verbose || verbose2
        println(STDERR, "regions does not contain $in_exons, from read $read_name")
        println(verboseOUT, "regions does not contain $in_exons, $in_snps, from read $read_name")
      end
    end
  end
  close(verboseOUT)
end

function haplo_lik(x::AbstractVector, h::AbstractVector, q::AbstractVector, rho::Float64)
  idx = find((x) -> x != 3, x)
  value1 = 0.0
  value2 = 0.0
  for j in idx
    qj = max(min(q[j], 1.0),0.0)
    if x[j] == h[j]
      value1 += log(qj)
      value2 += log((1 - qj)/3)
    elseif x[j] == (1 - h[j])
      value1 += log((1 - qj)/3)
      value2 += log(qj)
    else
      value1 += log((1 - qj)/3)
      value2 += log((1 - qj)/3)
    end
  end
  #return rho * exp(value1) / (rho * exp(value1) + (1-rho) * exp(value2))
  return exp(value1) / (exp(value1) + exp(value2))
end

function add_phase_reads!{V<:AbstractString,W<:Integer}(loci_dict::Dict{V,IsoformPhaseEntry{V,W}},prefix::AbstractString,temp_dir::AbstractString, num_src_types::Int; whole_gpd::Bool=false)
  gpd_file_name = ""
  if whole_gpd
    gpd_file_name = "$(temp_dir)/$(prefix)_gpd.bed"
  else
    gpd_file_name = "$(temp_dir)/$(prefix)_gpd_with_snps.bed"
  end
  for i in 1:num_src_types
    for line in eachline(open(`bedtools intersect -wo -a $(temp_dir)/$(prefix)_psl$(i).bed -b $gpd_file_name`)[1])
      fields = split(strip(line),"\t")
      block_sizes = [parse(Int64,s) for s in split(strip(fields[22],','),",")]
      q_starts = [parse(Int64,s) for s in split(strip(fields[23],','),",")]
      t_starts = [parse(Int64,s) for s in split(strip(fields[24],','),",")]
      num_blocks = length(block_sizes)

      read_name = fields[13]
      gene_name = fields[40]
      if !haskey(loci_dict, gene_name)
        continue
      end
      #add psl entry
      if haskey(loci_dict[gene_name].read_names,read_name)
        continue
      else
        target_coords  = vcat([ collect((t_starts[k]+1):(t_starts[k]+block_sizes[k])) for k in 1:length(t_starts)]...)
        if whole_gpd || length(intersect(target_coords,loci_dict[gene_name].snps))>0
          loci_dict[gene_name].read_names[read_name] = gene_name
          push!(loci_dict[gene_name].reads, PSLEntry{ASCIIString,Int64}([parse(Int64,s) for s in fields[4:11]]...,fields[12],fields[13],[parse(Int64,s) for s in fields[14:16]]...,fields[17],[parse(Int64,s) for s in fields[18:21]]...,block_sizes,q_starts,t_starts))
        end
      end
    end
  end
end

function phase_isoform_sub(y::Vector{Int}, ls::Vector{Int}, lk::Vector{Int}, C::Matrix{Int})
  M = sum(y)
  K = length(lk)
  function negloglik(theta::Vector)
    value = 0.0
    lambda = ls .* C * exp(theta)
    return -sum(y .* log(lambda)) + sum(lambda)
  end
  function negscore!(theta::Vector, storage::Vector)
    for j in 1:K
      storage[j] = -sum((y .* C[:,j] * exp(theta[j])) ./ (C * exp(theta)) - ls .* C[:,j] * exp(theta[j]))
    end
  end
  function neghessian!(theta::Vector, storage::Matrix)
    for j in 1:K
      for l in 1:K
        storage[j,l] = sum((y .* C[:,j] .* C[:,l] * exp(theta[j]) * exp(theta[l])) ./ (C * exp(theta)) .^ 2)
        if j == l
          storage[j,l] += sum(ls .* C[:,j] * exp(theta[j]))
        end
      end
    end
  end
  function negloglik_and_score!(theta::Vector, storage)
    value = 0.0
    lambda = ls .* C * exp(theta)
    for j in 1:K
      storage[j] = -sum((y .* C[:,j] * exp(theta[j])) ./ (C[:,:] * exp(theta)) - ls .* C[:,j] * exp(theta[j]))
    end
    return -sum(y .* log(lambda)) + sum(lambda)
  end
  #d3 = TwiceDifferentiableFunction(negloglik, negscore!, neghessian!)
  init = log(M * ones(K)/K ./ lk)
  #optimize(d3, init, iterations=1000, method=:newton)
  optimize(negloglik, negscore!, neghessian!, init, Newton(), OptimizationOptions(iterations=1000))
end
function phase_isoform_sub_sim(y::Vector{Int}, ls::Vector{Int}, lk::Vector{Int}, C::Matrix{Int}, map::Vector{Int})
  M = sum(y)
  N = maximum(map)
  K = length(lk)
  theta_0 = rand(Gamma(1,1), K)
  tau_0 = zeros(Float64, K)
  for t in 1:N
    idx = find(map .== t)
    if length(idx) == 1
      tau_0[t] = 1
    elseif length(idx) == 2
      tau_0[idx[1]] = theta_0[idx[1]] / sum(theta_0[idx])
      tau_0[idx[2]] = theta_0[idx[2]] / sum(theta_0[idx])
    end
  end
  lambda_0 = ls .* C * theta_0
  sim_y = [rand(Poisson(s)) for s in lambda_0]
  function negloglik(theta::Vector)
    value = 0.0
    lambda = ls .* C * exp(theta)
    return -sum(y .* log(lambda)) + sum(lambda)
  end
  function negscore!(theta::Vector, storage::Vector)
    for j in 1:K
      storage[j] = -sum((y .* C[:,j] * exp(theta[j])) ./ (C * exp(theta)) - ls .* C[:,j] * exp(theta[j]))
    end
  end
  function neghessian!(theta::Vector, storage::Matrix)
    for j in 1:K
      for l in 1:K
        storage[j,l] = sum((y .* C[:,j] .* C[:,l] * exp(theta[j]) * exp(theta[l])) ./ (C * exp(theta)) .^ 2)
        if j == l
          storage[j,l] += sum(ls .* C[:,j] * exp(theta[j]))
        end
      end
    end
  end
  function negloglik_and_score!(theta::Vector, storage)
    value = 0.0
    lambda = ls .* C * exp(theta)
    for j in 1:K
      storage[j] = -sum((y .* C[:,j] * exp(theta[j])) ./ (C[:,:] * exp(theta)) - ls .* C[:,j] * exp(theta[j]))
    end
    return -sum(y .* log(lambda)) + sum(lambda)
  end
  #d3 = TwiceDifferentiableFunction(loglik, score!, hessian!)
  init = log(sum(sim_y) * ones(K)/K ./ lk)
  #return optimize(d3, init, iterations=1000, method=:newton), theta_0, tau_0, sim_y
  #d3 = TwiceDifferentiableFunction(negloglik, negscore!, neghessian!)
  #init = log(M * ones(K)/K ./ lk)
  #optimize(d3, init, iterations=1000, method=:newton)
  return optimize(negloglik, negscore!, neghessian!, init, Newton(), OptimizationOptions(iterations=1000)), theta_0, tau_0, sim_y
end
