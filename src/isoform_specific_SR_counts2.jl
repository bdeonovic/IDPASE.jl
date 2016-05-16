using Distributions
if length(ARGS) == 7
  gene_level_input_file = ARGS[1]
  gene_level_results_file = ARGS[2]
  isoform_level_input_file = ARGS[3]
  isoform_level_results_file = ARGS[4]
  regions_input_file = ARGS[5]
  use_true_haplotype = parse(Int64, ARGS[6]) == 1
  use_type = ARGS[7]

  gene_level_input = Dict{ASCIIString, ASCIIString}()
  gene_level_results = Dict{ASCIIString, Vector{Any}}()
  isoform_level_input = Dict{ASCIIString, ASCIIString}()
  isoform_level_results = Dict{ASCIIString, ASCIIString}()
  regions_input = Dict{ASCIIString, ASCIIString}()

  for line in readlines(open(gene_level_input_file, "r"))
    fields = split(strip(line),"\t")
    gene_name = fields[8]
    if !(gene_name in keys(gene_level_input))
      gene_level_input[gene_name] = line
    end
  end
  for line in readlines(open(gene_level_results_file, "r"))
    fields = split(strip(line),"\t")
    gene_name = fields[2]
    rho = fields[5] != "NaN" ? parse(Float64, fields[5]) : NaN
    h = "NA"
    if isnan(rho)
      continue
    else
      h = fields[4]
    end
    snps = [parse(Int64,s) for s in split(fields[3],",")]
    if !(gene_name in keys(gene_level_results))
      haplotype = [parse(Int64,s) for s in split(h,"")]
      if haplotype[1] == 1
        gene_level_results[gene_name] = [1-rho; snps; join(1-haplotype,"")]
      else
        gene_level_results[gene_name] = [rho; snps; h]
      end
    end
  end
  for line in readlines(open(isoform_level_input_file, "r"))
    fields = split(strip(line),"\t")
    gene_name = fields[2]
    if !(gene_name in keys(isoform_level_input))
      isoform_level_input[gene_name] = line #[convert(ASCIIString, s) for s in split(fields[10],",")]
    end
  end
  for line in readlines(open(isoform_level_results_file, "r"))
    fields = split(strip(line),"\t")
    isoform_name = fields[3]
    if !(isoform_name in keys(isoform_level_results))
      isoform_level_results[isoform_name] = line
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

  println(STDERR, "Processing complete...")
  #for (gene_name, isoform_names) in isoform_level_input
  for (gene_name, isoform_line) in isoform_level_input
    isoform_fields = split(strip(isoform_line),"\t")
    isoform_names = [convert(ASCIIString, s) for s in split(isoform_fields[10],",")]
    K = length(isoform_names)
    num_isoforms = parse(Int64, isoform_fields[3])
    num_regions = parse(Int64, isoform_fields[4])
    if num_regions == 0
      continue
    end
    regions_in_isoforms = reshape([parse(Int64, s) for s in split(isoform_fields[5],",")], num_regions, num_isoforms)
    isoform_map = [parse(Int64, s) for s in split(isoform_fields[7],",")]


    
    read_count_matrix = zeros(Int64,3,0,0)
    num_snps = 0
    num_reads = 0
    old_snps = Int64[]
    pseudo_h = Int64[]
    ref_counts = Int64[]
    alt_counts = Int64[]
    if gene_name in keys(gene_level_input)
      fields = split(strip(gene_level_input[gene_name]),"\t")
      num_snps = parse(Int64, fields[4])
      num_reads = parse(Int64, fields[5])
      read_counts = [parse(Int64, s) for s in split(fields[10],",")]
      old_snps = [parse(Int64, s) for s in split(fields[11],",")]
      X = reshape([parse(Int64, s) for s in split(fields[6],",")], num_snps, num_reads)
      pseudo_h = fill(2,num_snps)
      ref_counts = zeros(Int64,num_snps)
      alt_counts = zeros(Int64,num_snps)

      for i in 1:num_snps
        ref_counts[i] = sum(X[i,1:read_counts[1]] .==0)
        alt_counts[i] = sum(X[i,1:read_counts[1]] .==1)
        if ref_counts[i] > alt_counts[i]
          pseudo_h[i] = 0
        elseif ref_counts[i] < alt_counts[i]
          pseudo_h[i] = 1
        end
      end

      read_count_matrix = zeros(Int64, 3, num_snps, length(read_counts))
      n0 = 1
      for i in 1:length(read_counts)
        read_count_matrix[1, :, i] = mapslices(x-> sum(x .== 0), X[:, n0:sum(read_counts[1:i])],2)
        read_count_matrix[2, :, i] = mapslices(x-> sum(x .== 1), X[:, n0:sum(read_counts[1:i])],2)
        read_count_matrix[3, :, i] = mapslices(x-> sum(x .== 2), X[:, n0:sum(read_counts[1:i])],2)
        n0 = sum(read_counts[1:i]) + 1
      end
    else
      warn("$gene_name not in gene_level_input")
      continue
    end
    rho = 0.0
    snps = Int64[]
    h = ""
    if gene_name in keys(gene_level_results)
      result = gene_level_results[gene_name]
      rho = result[1]
      snps = result[2:(end-1)]
      h = result[end]
      #println(string(rho, "\t", snps,"\t", h))
    else
      warn("$gene_name not in gene_level_results")
      continue
    end

    idx = findin(old_snps, snps)
    pseudo_h = pseudo_h[idx]
    ref_counts = ref_counts[idx]
    alt_counts = alt_counts[idx] 
    haplotype = [parse(Int64, s) for s in split(h,"")]
    unreliable_pseudo = false
    if any(pseudo_h .== 2)
      pseudo_h = "NA"
    else
      if pseudo_h[1] == 1
        pseudo_h = 1 - pseudo_h
      end
      if sum(pseudo_h .== haplotype) != length(snps)
        unreliable_pseudo = true
      end
    end
    snps_in_isoforms = zeros(Int64, 0, K)
    reg_dict = Dict{Set{Int64}, Vector{Int64}}()
    regions = Array{Int64,1}[]
    if gene_name in keys(regions_input)
      #println(STDERR, gene_name)
      fields = split(strip(regions_input[gene_name]), "\t")
      snp_idx = [parse(Int64, s) for s in split(fields[8],",")]
      num_snps = length(snp_idx)
      snps_in_isoforms = reshape([parse(Int64, s) for s in split(fields[6],",")],num_snps, K) 
      regions = [ collect(eval(parse(string(strip(s,[')']),")"))))::Array{Int64,1} for s in split(fields[4],"),")]
      y = [parse(Int64, s) for s in split(fields[5],",")]
      #println(STDERR, regions)
      #println(STDERR, regions[1])
      #println(STDERR, regions[2])
      #println(STDERR, "++++++++++++")
      #snp_idx = sort(unique(abs(vcat(regions...)[find(x-> x<0, vcat(regions...))])))


      if length(snp_idx) != length(haplotype)
        warn("length of snp_idx $(length(snp_idx)) not equal to length of haplotype $(length(haplotype))")
        continue
      end
      for i in 1:length(regions)
        regions[i] = regions[i][sortperm(abs(regions[i]))]
        idx1 = find(x -> x in abs(regions[i]), snp_idx)
        idx2 = find(x -> x in snp_idx, abs(regions[i]))
        if length(idx1) > 0
          if all((regions[i][idx2] .< 0) .== (haplotype[idx1] .== 1))
            if Set(abs(regions[i])) in keys(reg_dict)
              reg_dict[Set(abs(regions[i]))] += [y[i], 0]
            else
              reg_dict[Set(abs(regions[i]))] = [y[i], 0]
            end
          elseif all((regions[i][idx2] .< 0) .== (haplotype[idx1] .== 0))
            if Set(abs(regions[i])) in keys(reg_dict)
              reg_dict[Set(abs(regions[i]))] += [0, y[i]]
            else
              reg_dict[Set(abs(regions[i]))] = [0, y[i]]
            end
          end
        end
      end
      #=
      for i in 1:length(regions)
        for j in 1:length(snp_idx)
          idx = findfirst(abs(regions[i]) .== snp_idx[j])
          if idx != 0
            if ((regions[i][idx] < 0) && (haplotype[j] == 1)) || ((regions[i][idx] > 0) && (haplotype[j] == 0))
              if Set(abs(regions[i])) in keys(reg_dict)
                reg_dict[Set(abs(regions[i]))] += [y[i], 0]
              else
                reg_dict[Set(abs(regions[i]))] = [y[i], 0]
              end
            else
              if Set(abs(regions[i])) in keys(reg_dict)
                reg_dict[Set(abs(regions[i]))] += [0, y[i]]
              else
                reg_dict[Set(abs(regions[i]))] = [0, y[i]]
              end
            end
          end
        end
      end
      =#
    else
      warn("$gene_name not in regions_input")
      continue
    end
    tau = zeros(Float64, K)
    theta = zeros(Float64, K)
    psi = zeros(Float64, K)
    for k in 1:K
      isoform_name = isoform_names[k]
      if isoform_name in keys(isoform_level_results)
        fields = split(strip(isoform_level_results[isoform_name]),"\t")
        tau[k] = fields[4] != "NaN" ? parse(Float64, fields[4]) : NaN
        theta[k] = parse(Float64, fields[5])
      else
        warn("$isoform_name not in isoform_level_results")
        continue
      end
    end
    for k in 1:K
      isoform_name = isoform_names[k]
      if !(isoform_name in keys(isoform_level_results))
        warn("$isoform_name not in isoform_level_results")
        continue
      end

      #=
      res = [0,0]
      for s in regions[find(mapslices(sum,regions_in_isoforms[:,find(isoform_map .== k)],2))]
        if Set(s) in keys(reg_dict)
          res += reg_dict[Set(s)]
        end
      end
      =#
      res = [0,0]
      regions0 = Array{Int64,1}[]
      ref = Float64[]
      alt = Float64[]
      ref0 = Int[]
      alt0 = Int[]
      idx_vec = Array{Int64,1}[]
      K_star = maximum(isoform_map)
      for i in 1:length(regions)
        s = regions[i]
        idx1 = find(x -> x in abs(s), snp_idx)
        if length(idx1) > 0
          #idx2 = find(mapslices(sum,snps_in_isoforms[idx1,:],1) .== length(idx1))
          #println(STDERR, 1:K_star)
          #println(STDERR, (K_star+1):(length(isoform_map)))
          #idx2 = find(regions_in_isoforms[i,1:K_star] + regions_in_isoforms[i,(K_star+1):end])
          idx2 = Int64[]
          for j in 1:K_star
            idx3 = find(isoform_map .== j)
            if sum(regions_in_isoforms[i,idx3]) > 0
              push!(idx2, j)
            end
          end
          if Set(s) in keys(reg_dict) && length(idx2) > 0 && all(s .> 0) && (k in idx2)
            res0 = reg_dict[Set(s)]
            push!(regions0, s)

            push!(ref, tau[k]/sum(tau[idx2]) * res0[1])
            push!(alt, (1 - tau[k])/sum(1 - tau[idx2]) * res0[2])

            push!(ref0, res0[1])
            push!(alt0, res0[2])

            push!(idx_vec, idx2)

            #=
            if !isfinite(tau[k]/sum(tau[idx2]) * res0[1])
              println(STDERR, res0[1], tau[k], tau[idx2])
            end
            if !isfinite(tau[k]/sum(tau[idx2]) * res0[2])
              println(STDERR, res0[2], tau[k], tau[idx2])
            end
            =#
          end
        end
      end
      res = [sum(ref), sum(alt)]
      #=
      res = [0,0]
      for s in regions[find(mapslices(sum,regions_in_isoforms[:,find(isoform_map .== k)],2))]
        idx1 = find(x -> x in abs(s), snp_idx)
        if length(idx1) == 0
          continue
        end
        idx2 = find(mapslices(sum,snps_in_isoforms[idx1,:],1) .== length(idx1))
        if length(idx2) == 0
          continue
        end
        if Set(s) in keys(reg_dict) && length(idx2) > 0
          #res += (tau[k] / sum(tau[idx2])) * reg_dict[Set(s)]
          res0 = reg_dict[Set(s)]
          res += [tau[k]/sum(tau[idx2]) * res0[1], (1 - tau[k])/sum(1 - tau[idx2]) * res0[2]]

          #=
          #a = zeros(K)
          #b = zeros(K)

          #a[idx2] = ones(K)[idx2]
          #b[idx2] = ((1 - tau)./tau)[idx2]

          a = ones(length(idx2))
          b = ((1 - tau[idx2])./tau[idx2])

          #println(STDERR,b)
          res0 = zeros(K)
          res0[idx2] = [a' ; b'] \ reg_dict[Set(s)]

          idx3 = find(res0 .< 0)
          idx4 = find(res0 .> 0)

          if length(idx3) > 0
            res0[idx4] -= sum(abs(res0[idx3]))/length(idx4)
            res0[idx3] = 0
          end
          res += [res0[k], ((1 - tau[k])/tau[k]) * res0[k]]
          =#
        end
      end
      =#
      n = sum(res)
      prob = [tau[k], 1 - tau[k]]
      stat = sum((res .- n * prob) .^2  ./ (n * prob))
      pval = 1 - cdf(Chisq(1), stat)
      psi[k] = res[1]/sum(res)

      output_string = "$k\t$isoform_name\t$gene_name\t$rho\t$(tau[k])\t$(theta[k])\t$h"
      output_string = string(output_string,"\t",join(regions0,","),"\t",join(idx_vec,","),"\t",join(ref,","),"\t",join(alt,","),"\t",join(ref0,","),"\t", join(alt0,","),"\t", res[1]/sum(res),"\t",stat,"\t",pval, "\t", unreliable_pseudo, "\t", join(pseudo_h,""),"\t",join(ref_counts,","),"\t",join(alt_counts,","))
      #=
      output_string = string(output_string, "\t", join(snps_in_isoforms[:,k],","))
      #for i in 1:length(read_counts)
      for i in 1:1
        #output_string = string(output_string,"\t",join(mapslices(x-> join(x,","),read_count_matrix[1:2,:,i],2),"\t"))
        #output_string = string(output_string,"\t",join(mapslices(x-> x[1]/sum(x),read_count_matrix[1:2,:,i],2),","))
        snp_ase = zeros(Float64, num_snps)
        for j in 1:num_snps
          snp_ase[j] = read_count_matrix[1,j,i] / sum(read_count_matrix[1:2,j,i])
        end
        output_string = string(output_string,"\t",join(snp_ase,","))
      end
      output_string = string(output_string,"\t",res[1]/sum(res),"\t",stat,"\t",pval)
      =#
      #=
      for i in 1:length(read_counts)
        theta_adj = zeros(Float64, num_snps)
        snp_ase = zeros(Float64, num_snps)
        for j in 1:num_snps
          idx = find(snps_in_isoforms[j,:])
          if snps_in_isoforms[j,k] == 1
            theta_adj[j] = theta[k] / sum(theta[idx])
          end
          snp_ase[j] = theta_adj[j] * read_count_matrix[1,j,i] / sum(read_count_matrix[1:2,j,i])
        end
        #output_string = string(output_string,"\t",join(mapslices(x-> join(x,","),round(Int64, theta * read_count_matrix[:,:,i]),2),"\t"))
        output_string = string(output_string,"\t",join(snp_ase,","))
      end
      =#
      println(STDOUT, output_string)
    end
  end
end
