extensions = quote
  import Distributions: length, size, DiscreteMultivariateDistribution

  #Allele Specific Haplotype Distribution
  immutable ASH <: DiscreteMultivariateDistribution
    n::Integer         #number of reads
    m::Integer         #number of SNPs
    Q::Matrix{Float64} #Quality matrix
    h::Vector{Int64}   #haplotype
    rho::Float64       

    function ASH(n::Integer, m::Integer, Q::Matrix{Float64}, h::Vector{Int64}, rho::Float64)
      if n <= 0
        throw(ArgumentError("n must be a positive integer."))
      end
      if m <= 0
        throw(ArgumentError("m must be a positive integer."))
      end
      if size(Q) != (m, n)
        throw(ArgumentError("Dimension of Q must match (m,n) = "string((m,n))"  but it is "string(size(Q))))
      end
      if size(h,1) != m
        throw(ArgumentError("Dimension of h must be m = $m."))
      end
      if !all( 0.0 .<= Q .<= 1.0 )
        throw(ArgumentError("All values of Q must be between 0 and 1."))
      end
      if !all( (h.==0) | (h.==1) )
        throw(ArgumentError("All values of h must be 0 or 1."))
      end
      if !(0.0 <= rho <= 1)
        throw(ArgumentError("rho = $rho must be between 0 and 1."))
      end
      new(round(Int64,n), round(Int64,m), Q, h, rho)
    end
  end

  function insupport{T<:Real}(d::ASH, X::AbstractVector{T}) 
    all( (X.==0) | (X.==1) | (X.==2) | (X.==3))
  end

  function logpdf{T<:Real}(d::ASH, X::AbstractVector{T}) 
    if !insupport(d,X)
      return -Inf
    end
    Xr = reshape(X,d.m,d.n)
    Q = max(min(d.Q, 1-eps()),eps())

    result = 0.0
    @inbounds for i in 1:d.n
      idx = find( (x) -> x != 3, Xr[:,i])

      #delta = [round(Int64,s) for s in (Xr[:,i] .== d.h)]

      sum1 = 0.0
      sum2 = 0.0
      @inbounds for j in idx
        #sum1 += (delta[j])*log(Q[j,i]) + (1-delta[j])*log(1-Q[j,i])
        #sum2 += (1-delta[j])*log(Q[j,i]) + (delta[j])*log(1-Q[j,i])
        if Xr[j,i] == d.h[j]
          sum1 += log(Q[j,i])
          sum2 += log((1 - Q[j,i])/3)
        elseif Xr[j,i] == (1 - d.h[j])
          sum1 += log((1 - Q[j,i])/3)
          sum2 += log(Q[j,i])
        else
          sum1 += log((1 - Q[j,i])/3)
          sum2 += log((1 - Q[j,i])/3)
        end
      end
      result += log( d.rho*exp(sum1) + (1-d.rho)*exp(sum2))
      #result += log(d.rho) + sum1
    end
    return result
  end

  pdf(d::ASH, X::AbstractVector) = exp(logpdf(d, X))
  length(d::ASH) = d.n*d.m
  export ASH
end
eval(Mamba, extensions)

function build_gamma_recursive!{V,E,W}(g::AbstractGraph{V,E},eweights::Vector{W},gamma::Vector{Vector{V}})
  if num_vertices(g) > 1
    parity, bestcut = min_cut(g,eweights)
    v1 = vertices(g)[parity]
    v2 = vertices(g)[!parity]
    push!(gamma,v1)
    push!(gamma,v2)

    g1 = inclist(v1,is_directed=false)
    g2 = inclist(v2,is_directed=false)

    edge_idx1 = Int64[]
    edge_idx2 = Int64[]

    covered_edges = Dict{Set{Int64},Bool}()
    for v in v1
      edges = out_edges(v,g)
      for e in edges
        s = source(e)
        t = target(e)
        if !haskey(covered_edges,Set([s,t]))
          if s in v1 && t in v1 
            add_edge!(g1,s,t)
            covered_edges[Set([s,t])] = true
            push!(edge_idx1,edge_index(e))
          end
        end
      end
    end
    build_gamma_recursive!(g1,eweights[edge_idx1],gamma)

    covered_edges = Dict{Set{Int64},Bool}()
    for v in v2
      edges = out_edges(v,g)
      for e in edges
        s = source(e)
        t = target(e)
        if !haskey(covered_edges,Set([s,t]))
          if s in v2 && t in v2 
            add_edge!(g2,s,t)
            covered_edges[Set([s,t])] = true
            push!(edge_idx2,edge_index(e))
          end
        end
      end
    end
    build_gamma_recursive!(g2,eweights[edge_idx2],gamma)
  elseif num_vertices(g) == 1
    push!(gamma,vertices(g))
  end
end

function build_gamma{W}(X::Matrix{W})
  m,n = size(X)
  g = simple_inclist(m, is_directed=false)
  wedges = Dict{Set{Int64},Int64}()
  for pair in combinations(1:m,2)
    wedges[Set([pair[1],pair[2]])] = sum((X[pair[1],:] .!= 2) & (X[pair[2],:] .!= 2))
  end
  w = length(wedges)
  eweights = Float64[]
  for (key,value) in wedges
    edge = collect(keys(key.dict))
    add_edge!(g, edge[1],edge[2])
    push!(eweights,value)
  end
  gamma = Array{Int64,1}[collect(1:m)]

  build_gamma_recursive!(g,eweights,gamma)
  unique(gamma)
end

function psl_to_bed(psl_file_name::AbstractString; smooth::Int64=68,chr::Vector{ASCIIString}=ASCIIString["chr$i" for i in [1:22;"X";"Y";"M"]],output_dir::AbstractString="",
  output_file_name::ASCIIString="psl.bed")

  num_reads = 0
  psl = open(psl_file_name,"r")
  out = [open(string(output_dir,"/",i,"_",output_file_name), "w") for i in chr]
  for line in eachline(psl)
    num_reads += 1
    fields = split(strip(line), "\t")

    if length(fields) < 21
      continue
    end

    out_idx = find((x)-> fields[14]==x,chr)
    if length(out_idx)!=1
      continue
    end

    #num_blocks = int(fields[18])
    block_sizes = [parse(Int64,s) for s in split(strip(fields[19],','),",")]
    block_starts = [parse(Int64,s) for s in split(strip(fields[21],','),",")]
    block_ends = block_starts + block_sizes - 1
    num_blocks = length(block_sizes)
    
    if num_blocks > 1
      idx = find(x-> x>smooth,block_starts[2:end] - block_ends[1:(end-1)]) 
      push!(idx,num_blocks)

      i = 1; j = 1
      while i < num_blocks
        print(out[out_idx[1]], string(fields[14],"\t",block_starts[i],"\t",block_ends[idx[j]],"\t",line))
        i = idx[j] + 1
        j = j + 1
      end
    elseif num_blocks == 1
      print(out[out_idx[1]], string(fields[14],"\t",block_starts[1],"\t",block_ends[1],"\t",line))
    end

  end
  close(psl)
  for i in 1:length(out)
    close(out[i])
  end
  return num_reads
end
function psl_to_bed2(psl_file_name::AbstractString; smooth::Int64=68,chr::Vector{ASCIIString}=ASCIIString["chr$i" for i in [1:22;"X";"Y";"M"]],output_dir::AbstractString="",
  output_file_name::ASCIIString="psl.bed")

  num_reads = 0
  psl = open(psl_file_name,"r")
  out = [open(string(output_dir,"/",i,"_",output_file_name), "w") for i in chr]
  for line in eachline(psl)
    num_reads += 1
    fields = split(strip(line), "\t")

    if length(fields) < 21
      continue
    end

    out_idx = find((x)-> fields[14]==x,chr)
    if length(out_idx)!=1
      continue
    end

    #num_blocks = int(fields[18])
    block_sizes = [parse(Int64,s) for s in split(strip(fields[19],','),",")]
    block_starts = [parse(Int64,s) for s in split(strip(fields[21],','),",")]
    block_ends = block_starts + block_sizes - 1
    num_blocks = length(block_sizes)
    
    for i in 1:num_blocks
      print(out[out_idx[1]], string(fields[14],"\t",block_starts[i],"\t",block_ends[i],"\t",line))
    end

  end
  close(psl)
  for i in 1:length(out)
    close(out[i])
  end
  return num_reads
end
function gpd_to_bed(gpd_file_name::AbstractString;chr::Vector{ASCIIString}=ASCIIString["chr$i" for i in [1:22;"X";"Y";"M"]],
                    output_dir::AbstractString="",to_keep::Dict{AbstractString,Int64}=Dict{AbstractString,Int64}(), output_file_name::ASCIIString="gpd.bed")
  gpd = open(gpd_file_name,"r")
  out = [open(string(output_dir,"/",i,"_",output_file_name), "w") for i in chr]
  to_keep_flag = (length(to_keep)>0)
  for line in eachline(gpd)
    if ismatch(r"^#",line)
      continue
    end
    fields = split(strip(line),"\t")
    block_ends   = [parse(Int64,s) for s in split(strip(fields[11],','),",")]
    block_starts = [parse(Int64,s) for s in split(strip(fields[10],','),",")]
    num_blocks = length(block_starts)
    out_idx = find((x)-> fields[3]==x,chr)
    if length(out_idx)!= 1
      continue
    end
    if to_keep_flag && !haskey(to_keep,fields[13])
      continue;
    end
    for i in 1:num_blocks
      println(out[out_idx[1]],string(fields[3],"\t",block_starts[i],"\t",block_ends[i],"\t",strip(line)))
    end

  end
  close(gpd)
  for i in 1:length(out)
    close(out[i])
  end
  return nothing
end
function vcf_to_bed(vcf_file_name::AbstractString;chr::Vector{ASCIIString}=ASCIIString["chr$i" for i in [1:22;"X";"Y";"M"]],
                      output_dir::AbstractString="",to_keep::Dict{AbstractString,Int64}=Dict{AbstractString,Int64}(), output_file_name::ASCIIString="vcf.bed")
  vcf = open(vcf_file_name,"r")
  out = [open(string(output_dir,"/",i,"_",output_file_name), "w") for i in chr]
  to_keep_flag = (length(to_keep)>0)
  for line in eachline(vcf)
    if ismatch(r"^#",line)
      continue
    end
    fields = split(strip(line),"\t")
    out_idx = find((x)-> fields[1]==x,chr)
    if length(out_idx)!= 1
      continue
    end
    if to_keep_flag && !haskey(to_keep,fields[1])
      continue;
    end

    println(out[out_idx[1]],"$(fields[1])\t$(parse(Int64,fields[2])-1)\t$(parse(Int64,fields[2]))\t"join(fields,"\t"))
  end
  close(vcf)
  for i in 1:length(out)
    close(out[i])
  end
  return nothing

end

function process_gpd(gpd_file_name::AbstractString)
  gpd_file = open(gpd_file_name,"r")
  isoforms = Dict{SubString{ASCIIString},Array{GPDEntry{SubString{ASCIIString},Int64},1}}()
  for line in eachline(gpd_file)
    if ismatch(r"^#",line)
      continue
    end
    fields = split(strip(line),"\t")
    block_ends   = [parse(Int64,s) for s in split(strip(fields[11],','),",")]
    block_starts = [parse(Int64,s) for s in split(strip(fields[10],','),",")]
    frames = [parse(Int64,s) for s in split(strip(fields[16],','),",")]
    isoform_length = sum(block_ends - block_starts)
    chr = fields[3]

    #m = match(r"chr[\dXYM]*",chr)
    #if m == nothing
    #  continue
    #end
    if haskey(isoforms,fields[13])
      #warn("Transcript $(fields[2]) already observed, skipping...")
      if isoforms[fields[13]][1].chr == chr
        push!(isoforms[fields[13]],GPDEntry(parse(Int64,fields[1]),fields[2],fields[3],fields[4],[parse(Int64,s)
        for s in
        fields[5:9]]...,block_starts,block_ends,parse(Int64,fields[12]),fields[13],fields[14],fields[15],frames))
      end
    else
      isoforms[fields[13]] =
      [GPDEntry(parse(Int64,fields[1]),fields[2],fields[3],fields[4],[parse(Int64,s) for s in fields[5:9]]...,block_starts,block_ends,parse(Int64,fields[12]),fields[13],fields[14],fields[15],frames)]
    end
  end
  close(gpd_file)
  return isoforms
end

function print_gpd(io::IO,entry::GPDEntry, bed_start::Integer=0, bed_end::Integer=0)
  println(io, string(entry.chr,"\t",
  entry.tx_start,"\t",
  entry.tx_end,"\t"),
  #bed_start,"\t",
  #bed_end,"\t"),
  string(entry.bin,"\t",
  entry.name,"\t",
  entry.chr,"\t",
  entry.strand,"\t",
  entry.tx_start,"\t",
  entry.tx_end,"\t",
  entry.cds_start,"\t",
  entry.cds_end,"\t",
  entry.exon_count,"\t",
  join(entry.exon_starts,","),"\t",
  join(entry.exon_ends,","),"\t",
  entry.score,"\t",
  entry.gene_name,"\t",
  entry.cds_start_stat,"\t",
  entry.cds_end_stat,"\t",
  join(entry.exon_frames,",")))
end


function initialize_loci(gpd_file_name::AbstractString, vcf_file_name::AbstractString, temp_dir::AbstractString, chr::AbstractString, num_src_types::Int,  fpkm_dict::Dict{ASCIIString, Array{Float64, 1}})
  gpd_file = open(gpd_file_name,"r")
  gpd_with_snps = open("$(temp_dir)/$(chr)_gpd_with_snps.bed","w")

  isoforms = Dict{ASCIIString,Array{GPDEntry{ASCIIString,Int64},1}}()
  loci_dict = Dict{ASCIIString,  PhaseEntry{ASCIIString, Int64}}()
  transcripts = Dict{ASCIIString,Bool}()

  for line in eachline(open(`bedtools intersect -wo -a $gpd_file_name -b $vcf_file_name`)[1])
    fields = split(strip(line),"\t")
    block_ends   = [parse(Int64,s) for s in split(strip(fields[14],','),",")]
    block_starts = [parse(Int64,s) for s in split(strip(fields[13],','),",")]
    frames = [parse(Int64,s) for s in split(strip(fields[19],','),",")]
    isoform_length = sum(block_ends - block_starts)

    if (length(keys(fpkm_dict)) > 0) 
      if (fields[5] in keys(fpkm_dict))
        fpkm = fpkm_dict[fields[5]]
        if fpkm[1] < 0.01 * fpkm[2]
          continue
        end
      else
        continue
      end
    end


    if haskey(transcripts,fields[5])
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
          error("Bad vcf")
          continue
        end
        #sort snps
        idx = sortperm(loci.snps)
        loci.ref = loci.ref[idx]
        loci.alt = loci.alt[idx]
        loci.snps = loci.snps[idx]
        loci.true_h = loci.true_h[idx]
      end
      loci_dict[fields[16]] = loci
    else
      isoform = GPDEntry{ASCIIString,Int64}(parse(Int64,fields[4]),fields[5],fields[6],fields[7],[parse(Int64,s) for s in fields[8:12]]...,block_starts,block_ends,parse(Int64,fields[15]),fields[16],fields[17],fields[18],frames)
      if haskey(isoforms,fields[16])
        #warn("Transcript $(fields[2]) already observed, skipping...")
        if isoforms[fields[16]][1].chr == chr
          push!(isoforms[fields[16]],isoform)
        end
      else
        isoforms[fields[16]]   = [isoform]
        loci =PhaseEntry{ASCIIString, Int64}(chr,Int64[],ASCIIString[],ASCIIString[],Dict{ASCIIString,Bool}(),Dict{ASCIIString,ASCIIString}(),Int64[],true,PSLEntry{ASCIIString,Int64}[],ASCIIString[],fill(2,0,0),fill(1.0,0,0),zeros(Int64,num_src_types))
        loci_dict[fields[16]]  = loci
      end
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
          continue
        end
      end
      print_gpd(gpd_with_snps,isoform)
      transcripts[fields[5]] = true
    end # in transcripts? 
  end #For loop
  close(gpd_file)
  close(gpd_with_snps)
  return isoforms, loci_dict
end


function add_reads!{V<:AbstractString,W<:Integer}(loci_dict::Dict{V,PhaseEntry{V,W}},prefix::AbstractString,temp_dir::AbstractString, num_src_types::Int; whole_gpd::Bool=false)
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
        target_coords  = vcat([ collect((t_starts[k]+round(Int64,block_sizes[k]*.1)):(t_starts[k]+round(Int64, block_sizes[k]*.90))) for k in 1:length(t_starts)]...)
        if whole_gpd || length(intersect(target_coords,loci_dict[gene_name].snps))>0
          loci_dict[gene_name].read_names[read_name] = gene_name
          push!(loci_dict[gene_name].reads, PSLEntry{ASCIIString,Int64}([parse(Int64,s) for s in fields[4:11]]...,fields[12],fields[13],[parse(Int64,s) for s in fields[14:16]]...,fields[17],[parse(Int64,s) for s in fields[18:21]]...,block_sizes,q_starts,t_starts))
        end
      end
    end
  end
end

# format = 1 is Phred+33, format = 2 is Phred+64
function extract_X_Q!(fastq_file_name::AbstractString, loci::Dict, reads::Dict, format::Integer=1, src::Int=1, verbose::Bool=false)
  fastq = open(fastq_file_name,"r")
  i = 1
  read_name = ""
  sequence = ""
  strand = ""
  quality = ""
  for line in eachline(fastq)
    if i == 4
      quality = strip(line)
      if haskey(reads,read_name)
        for loci_name in reads[read_name]
          psl_entries = loci[loci_name].reads
          for j in 1:length(psl_entries)
            if psl_entries[j].q_name == read_name
              loci[loci_name].read_counts[src] += 1
              push!(loci[loci_name].read_names_order, psl_entries[j].q_name)

              q_starts = psl_entries[j].q_starts +1
              q_ends = psl_entries[j].q_starts + psl_entries[j].block_sizes
              t_starts = psl_entries[j].t_starts + 1
              t_ends = psl_entries[j].t_starts + psl_entries[j].block_sizes
              blocks = length(q_starts)

              #1-based mapped coords
              query_coords  = vcat([ collect(q_starts[k]:q_ends[k]) for k in 1:blocks]...)
              target_coords = vcat([ collect(t_starts[k]:t_ends[k]) for k in 1:blocks]...)

              snps = loci[loci_name].snps  #1-based
              num_snps = length(snps)

              target_snp_pos = zeros(Int64,num_snps)

              for k in 1:num_snps
                target_snp_pos_k = find( (x)-> x == snps[k],target_coords)
                if length(target_snp_pos_k) == 1
                  target_snp_pos[k] = target_snp_pos_k[1]
                end
              end

              snps_detected = find( (x)-> x!=0,target_snp_pos)
              num_snps_detected = length(snps_detected)
              query_snp_pos = query_coords[target_snp_pos[snps_detected]]

              snp_calls = fill("_",num_snps)
              from_fasta = ASCIIString[]
              from_quality = ""
              if verbose
                println(loci_name)
                println(sequence)
                println(query_snp_pos)
              end
              if psl_entries[j].strand == "+"
                from_fasta = split(sequence[query_snp_pos],"")
                from_quality = quality[query_snp_pos]
              else
                from_fasta = split(convert(ASCIIString,reversecomplement(convert(ASCIIString,sequence),DNA_COMPLEMENT)[query_snp_pos]),"")
                from_quality = quality[reverse(length(quality) - query_snp_pos + 1)]
              end
              h = fill(3,num_snps)
              q = fill(0.0,num_snps)
              format_adj = 0
              if format == 1
                format_adj = 33
              elseif format == 2
                format_adj = 64
              else
                error("Improper format")
              end

              for k in 1:num_snps_detected
                if from_fasta[k] == loci[loci_name].ref[snps_detected[k]]
                  h[snps_detected[k]] = 0
                elseif from_fasta[k] == loci[loci_name].alt[snps_detected[k]]
                  h[snps_detected[k]] = 1
                else
                  h[snps_detected[k]] = 2
                end
                q[snps_detected[k]] = 1 - 10 ^ ((Int(from_quality[k]) - format_adj) / (-10))
              end
              if all(h.==3)
                loci[loci_name].read_counts[src] -= 1
              end
              if size(loci[loci_name].X) == (0,0)
                loci[loci_name].X = fill(2,0,num_snps)
              end
              if size(loci[loci_name].Q) == (0,0)
                loci[loci_name].Q = fill(0,0,num_snps)
              end
              loci[loci_name].X = vcat(loci[loci_name].X,h')
              loci[loci_name].Q = vcat(loci[loci_name].Q,q')
            end
          end
        end
      end
      i = 0
    elseif i == 3
      strand = strip(line)
    elseif i == 2
      sequence = strip(line)
    elseif i == 1
      m = match(r"@([\w\.\:\/]+)",strip(line))
      if m != nothing
        read_name = m.captures[1]
      else
        print(STDERR, line)
        error("Improper FASTQ file")
      end
    end
    i = i + 1
  end
  close(fastq)
end

function cleanup!(loci::Dict, max_snps=Inf)
  result = PHASEDataReal{keytype(loci),Float64,Int64}[]
  for (key,value) in loci
    idx = find((x) -> !all(value.X[x,:].== 3),1:size(value.X,1))
    value.X = value.X[idx,:]
    value.Q = value.Q[idx,:]
    value.read_names_order = value.read_names_order[idx]

    if size(value.X,1) == 0 || length(value.ref) > max_snps
      delete!(loci,key)
      continue
    end
    g = build_gamma(value.X')
    push!(result, PHASEDataReal{keytype(loci),Float64,Int64}(value.read_counts,value.snps,convert(SubString{ASCIIString},key),convert(SubString{ASCIIString},"NA"),value.chr,value.read_names_order,value.ref,value.alt,value.X',value.Q',g,value.true_h))
  end
  return result
end

function simulate_data{U<:AbstractString,V<:AbstractFloat,W<:Integer}(real_data::Vector{PHASEDataReal{U,V,W}},out_file_name::U, sim_rho::Bool=true)
  sim_data = PHASEDataSim{U,V,W}[]
  out = open(out_file_name,"w")
  for k in 1:length(real_data)
    m,n = size(real_data[k].X)
    true_h = [0,rand(0:1,m-1)]
    true_rho = sim_rho ? min(max(rand(Normal(0.5,0.12)),0.0),1.0) : 0.5
    true_delta = zeros(n)


    X = fill(3,m,n)
    for i in 1:n
      idx = find( (x)-> x != 3, real_data[k].X[:,i])
      if rand() < true_rho
        true_delta[i] = 1
        for j in idx
          X[j,i] = wsample([true_h[j],1-true_h[j],2], [real_data[k].Q[j,i], (1 - real_data[k].Q[j,i])/3, 2 * (1 - real_data[k].Q[j,i])/3])
        end
      else
        for j in idx
          X[j,i] = wsample([1-true_h[j],true_h[j],2], [real_data[k].Q[j,i], (1 - real_data[k].Q[j,i])/3, 2 * (1 - real_data[k].Q[j,i])/3])
        end
      end
    end
    idx = find(mapslices( (x)->!all(x.==3),X,1)[:])
    X = X[:,idx]
    Q = real_data[k].Q[:,idx]
    println(out,string(k,"\t",true_rho,"\t",join(true_h,""),"\t",size(X,1),"\t",size(X,2),"\t", join(X[:],","),"\t",join(Q[:],","),"\t",real_data[k].gene_name,"\t",real_data[k].isoform_name,"\t",join(real_data[k].read_counts,","),"\t",join(real_data[k].coords,","),"\t",join(real_data[k].reference_alleles,","),"\t",join(real_data[k].alternative_alleles,","),"\t",real_data[k].chr))
    push!(sim_data,PHASEDataSim{U,V,W}(real_data[k].read_counts, real_data[k].coords, real_data[k].gene_name, real_data[k].isoform_name, real_data[k].chr, real_data[k].read_names, real_data[k].reference_alleles, real_data[k].alternative_alleles, X, Q, real_data[k].gamma, true_h, true_rho, true_delta, Float64[]))
  end
  close(out)
  return sim_data
end

function read_nth_lines(stream, num)
  for i = 1:num-1
    readline(stream)
  end
  return readline(stream)
end

function phase_isoform(isoform_data_file::AbstractString, line_num::Integer, out_file_name::AbstractString, sim::Bool, nsim::Int=10)
        
  IN = open(isoform_data_file,"r")
  line = read_nth_lines(IN, line_num)
  close(IN)

  fields = split(strip(line),"\t")
  id = fields[1]
  gene_name = fields[2]
  num_isoforms = parse(Int64, fields[3])
  num_regions = parse(Int64, fields[4])
  C_vec = [parse(Int64, s) for s in split(fields[5],",")]
  C = reshape(C_vec, num_regions, num_isoforms)
  y = [parse(Int64, s) for s in split(fields[6],",")]
  isoform_map = [parse(Int64, s) for s in split(fields[7],",")]
  lk = [parse(Int64, s) for s in split(fields[8],",")]
  ls = [parse(Int64, s) for s in split(fields[9],",")]
  isoform_names = split(fields[10],",")

  OUT = open(out_file_name,"w")
  try 
    opt = phase_isoform_sub(y, ls, lk, C)
    conv = any([opt.iteration_converged, opt.f_converged, opt.gr_converged, opt.x_converged])
    theta = exp(opt.minimum)

    function get_isof_rho(isoform::Int64)
      idx = find(x-> x==isoform, isoform_map)
      return theta[idx[1]]/sum(theta[idx])
    end

    tau = map(get_isof_rho,1:maximum(isoform_map))
    rho = sum(theta[1:maximum(isoform_map)])/sum(theta)

    println(OUT,string(gene_name,"\t",rho,"\t",join(isoform_names,","),"\t",join(tau,","),"\t", join(theta,","),"\t",join(isoform_map,","),"\t",sum(y), "\t", join(lk,",") ))
  catch e
    showerror(STDERR, e, backtrace())
    println()
    println(STDERR, "Failed to optimize")
    println(OUT,string(gene_name,"\t","NA","\t",join(isoform_names,","),"\t","NA","\t", sum(y), "\t", join(lk,","),"\t","NA"))
  end
  close(OUT)

  if sim
    for i in 1:nsim
      SIMOUT_file = replace(isoform_data_file,r"(.*test_in/)(.*).txt$",Base.SubstitutionString("\\1simase/\\2::$gene_name::$i.txt"))
      SIMOUT = open(SIMOUT_file, "w")
      OUT = open(string(out_file_name,"_sim_$i"), "w")
      try
        opt, theta_0, tau_0,sim_y = phase_isoform_sub_sim(y, ls, lk, C, isoform_map)
        println(SIMOUT, join([fields[1:5]; join(sim_y,","); fields[7:end]],"\t"))

        conv = any([opt.iteration_converged, opt.f_converged, opt.gr_converged, opt.x_converged])
        theta = exp(opt.minimum)

        function get_isof_rho(isoform::Int64)
          idx = find(x-> x==isoform, isoform_map)
          return theta[idx[1]]/sum(theta[idx])
        end

        tau = map(get_isof_rho,1:maximum(isoform_map))
        rho = sum(theta[1:maximum(isoform_map)])/sum(theta)

        println(OUT,string(gene_name,"\t",rho,"\t",join(isoform_names,","),"\t",join(tau,","),"\t", join(theta,","),"\t",join(isoform_map,","),"\t",sum(y), "\t", join(lk,","),"\t",join(tau_0[1:maximum(isoform_map)::Int64],","),"\t",join(theta_0,",") ))
      catch e
        showerror(STDERR, e, backtrace())
        println()
        println(OUT,string(gene_name,"\t","NA","\t",join(isoform_names,","),"\t","NA","\t", sum(y), "\t", join(lk,","),"\t","NA","\t","NA","\t","NA" ))
        println(SIMOUT, join([fields[1:5]; join(["NA" for j in 1:length(y)],","); fields[7:end]],"\t"))
      end
      close(OUT)
      close(SIMOUT)
    end
  end
end

function read_data(true_data_file::AbstractString, read_file::AbstractString, line_num::Integer=1)
  IN = open(true_data_file,"r")
  line = read_nth_lines(IN, line_num)
  close(IN)

  fields = split(strip(line),"\t")
  id = fields[1]
  rho0 =  fields[2]
  if rho0 != "NA"
    rho0 = parse(Float64, fields[2])
  end
  h0 = [parse(Int64,s) for s in split(fields[3],"")]
  m = parse(Int64, fields[4])
  n = parse(Int64, fields[5])
  Xvec = [parse(Int64,s) for s in split(fields[6],",")]
  X = reshape(Xvec, m, n)
  Qvec = [parse(Float64,s) for s in split(fields[7],",")]
  Q = reshape(Qvec, m, n)
  gene_name = fields[8]
  isoform_name = fields[9]
  read_counts = [parse(Int64, s) for s in split(fields[10],",")]
  snps = [parse(Int64, s) for s in split(fields[11],",")]
  ref = split(fields[12],",")
  alt = split(fields[13],",")
  chr = fields[14]


  IN = open(read_file,"r")
  line = read_nth_lines(IN, line_num)
  close(IN)

  fields = split(strip(line),"\t")
  id = fields[1]
  read_names = split(fields[2],",")

  gamma = build_gamma(X)
  PHASEDataReal{ASCIIString,Float64,Int64}(read_counts, snps, convert(ASCIIString,gene_name), convert(ASCIIString,isoform_name), convert(ASCIIString,chr), convert(Array{ASCIIString,1},read_names), convert(Array{ASCIIString,1}, ref), convert(Array{ASCIIString,1}, alt), X, Q, gamma, h0)
end

function convert_gamma{W<:Integer}(gamma::Vector{Vector{W}},d::W)
  gamma_alt = Array{Int64,1}[]
  for set in gamma
    new_set = setdiff(collect(1:d),set)
    if 1 in set 
      if length(new_set) > 0
        push!(gamma_alt,new_set.-1)
      end
    else
      push!(gamma_alt,set.-1)
    end
  end
  return unique(gamma_alt)
end
function phase(data::PHASEData,max_iters::Int64,burnin::Int64,nchains::Int64,method::Int64, reads_to_use::Vector{Int64})
  PHASE = Dict{Symbol,Any}()

  PHASE[:m], PHASE[:n] = size(data.X)
  PHASE[:Q] = data.Q
  PHASE[:X] = data.X

  PHASE[:X] = PHASE[:X][:,reads_to_use]
  PHASE[:Q] = PHASE[:Q][:,reads_to_use]
  PHASE[:n] = length(reads_to_use)

  if PHASE[:m] == 1
    model = Model(
      X = Stochastic(1, (n,m,Q,hfull,rho) ->
            ASH(n,m,Q,convert(Array{Int64,1},hfull),convert(Float64,rho)),
          false
      ), 
      hfull = Logical(1, () -> Int64[0],true),
      rho = Stochastic(() -> Uniform(0.0,1.0), true)
    )
    inits = [ Dict{Symbol,Any}(:X => reshape(PHASE[:X],PHASE[:n]*PHASE[:m],1)[:], :rho => rand() ) for i in 1:nchains]
    scheme = [ Slice([:rho],[0.5],transform=true)]
    setsamplers!(model,scheme)
    tic()
    sim =  mcmc(model, PHASE, inits, burnin+500, burnin=burnin, thin=1, chains=nchains)
    run_time = toc()
    while gelmandiag(sim[:,"rho",:]).value[1,2,1] > 1.2 && sim.model.iter < max_iters && run_time < 3600
      tic()
      sim = mcmc(sim, 500)
      run_time += toc()
    end
    return sim, run_time
  else
    model = Model(
      X = Stochastic(1, (n,m,Q,hfull,rho) ->
            ASH(n,m,Q,convert(Array{Int64,1},hfull),convert(Float64,rho)),
          false
      ), 
      h = Stochastic(1, (m) ->
            UnivariateDistribution[DiscreteUniform(0,1) for i in 1:convert(Int64,m-1)],            
          false
      ),
      hfull = Logical(1, (h) ->
                [0 ; convert(Array{Int64,1},h)],
              true
      ),
      rho = Stochastic( () -> Uniform(0,1.0), true)
    )
    inits = Dict{Symbol,Any}[ Dict{Symbol,Any}(:X => reshape(PHASE[:X],PHASE[:n]*PHASE[:m],1)[:], :h => [rand(0:1) for j in 1:(PHASE[:m]-1)], :rho => rand() ) for i in 1:nchains]
    scheme = NaN
    if method == 1
      scheme = [ BMC3([:h], indexset = convert_gamma(data.gamma,PHASE[:m])), Slice([:rho],[0.5],transform=true)]
    elseif method == 2
      scheme = [ BMG([:h]), Slice([:rho],[0.5],transform=true)]
    elseif method == 3
      scheme = [ BHMC([:h], (2*PHASE[:m]+0.5)*pi), Slice([:rho],[0.5],transform=true)]
    else
      error("No such Method")
    end
    setsamplers!(model,scheme)
    tic()
    sim =  mcmc(model, PHASE, inits, burnin+500, burnin=burnin, thin=1, chains=nchains)
    run_time = toc()
    while gelmandiag(sim[:,"rho",:]).value[1,2,1] > 1.2 && sim.model.iter < max_iters && run_time < 18000
      tic()
      sim = mcmc(sim, 500)
      run_time += toc()
    end
    return sim, run_time
  end
end
function haplotype_map(sim::ModelChains)
  haplo_modes = Dict{Array{Int64,1},Array{Float64,1}}()
  num_iters, num_vars, num_chains = size(sim.value)
  max = 0
  max_h = Int64[]
  for i in 1:num_chains
    for j in 1:num_iters
      h = vec(round(Int64,sim[:,:hfull,:].value[j,:,i]))
      if h in keys(haplo_modes)
        haplo_modes[h] += [1; sim[:,:rho,:].value[j,:,i][1]]
      else
        haplo_modes[h] = [1; sim[:,:rho,:].value[j,:,i][1]]
      end
      if haplo_modes[h][1] > max
        max = haplo_modes[h][1]
        max_h = h[:]
      end
    end
  end
  return (max_h, haplo_modes[max_h][2]/haplo_modes[max_h][1])
end

function run_MCMC(args::Tuple{PHASEData,Int64,Int64,Int64,Int64,Vector{Bool},Vector{ASCIIString}})
  data, iters, burnin, nchains, method, use_src, type_names = args
  n0 = 1
  reads_to_use = Int64[]
  type_name = join(type_names[find(use_src)],"")
  for s in 1:length(use_src)
    if use_src[s]
      reads = n0:(n0 + data.read_counts[s] -1)
      num_lr = length(reads)
      append!(reads_to_use, reads)
    end
    n0 += data.read_counts[s]
  end

  if !any(use_src)
    error("must use atleast of one of the data types")
  end

  m, n = size(data.X)
  n = length(reads_to_use)

  if n > 0
    sim, run_time  = phase(data,iters,burnin, nchains, method, reads_to_use)

    d = summarystats(sim)
    q = quantile(sim)


    rhom = q.value[1,3]
    g = gelmandiag(sim[:,["rho"],:]).value[1,2,1]

    lth = mean(sim.value[:,1,:] .< 0.5)
    gth = mean(sim.value[:,1,:] .> 0.5)

    h_0 = [round(Int64,s) for s in d.value[2:(m+1),1]]

    kdrho = kde(vec(sim[:,["rho"],:].value))
    rho_mode = kdrho.x[findmax(kdrho.density)[2]]

    h_mode, cond_rho = haplotype_map(sim)

    result1 = [length(data.read_counts),join(data.read_counts,","),n,m,data.gene_name,data.isoform_name, join(data.true_h,""),join(h_0,""), join(h_mode,""), cond_rho, rho_mode, q.value[1,3], d.value[1,1], q.value[1,1], q.value[1,5],hpd(sim).value[1,1],hpd(sim).value[1,2],d.value[1,5],run_time,g,lth,gth,sim.model.iter, sim.model.burnin, join(data.coords,","), type_name]  
    return result1
  else
    result1 = [length(data.read_counts); join(data.read_counts,","); n; m; data.gene_name; data.isoform_name; join(data.true_h,""); ["NA" for i in 1:18]; type_name]  
    return result1
  end
end
function make_X_Q(args::Tuple{AbstractString, Dict, Bool, Bool})
  chr, parsed_args, sim, isoform = args
  num_src_types = length(parsed_args["psl"])
  read_dict = Dict{ASCIIString,Set{ASCIIString}}()
  fpkm_dict = Dict{ASCIIString,Array{Float64,1}}()
  if parsed_args["fpkm"] != nothing
    @time fpkm_dict = get_fpkm(parsed_args["fpkm"])
  end
  if isoform
    println(STDERR, "Getting isoforms with SNPS...")
    @time isoform_dict, loci_dict, transcript_dict = init_loci("$(parsed_args["temp"])/$(chr)_gpd.bed", chr, num_src_types, fpkm_dict); 
    @time get_snps!(loci_dict, "$(parsed_args["temp"])/$(chr)_gpd.bed", "$(parsed_args["temp"])/$(chr)_vcf.bed", parsed_args["temp"], chr, fpkm_dict)
    @time get_gene_level_results(loci_dict, parsed_args["results"], parsed_args["type"], !parsed_args["estimate"])

    if false
      REGIONS = open(string(parsed_args["out"],"/",parsed_args["prefix"],"_inter","_$(chr).txt"),"w")
      i = 1
      for (gene_name, loci) in loci_dict
        println(REGIONS, string(i, "\t", gene_name,"\t", join(loci.snps, ","),"\t", join(loci.est_h,""),"\t", loci.rho_mode ))
        i += 1
      end
      close(REGIONS)
    else
      @time assemble_loci!(loci_dict,isoform_dict, read_length = parsed_args["read_length"], min_overlap_length = 1, mean_insert_length = 0);
      @time add_phase_reads!(loci_dict, chr, parsed_args["temp"], num_src_types)

      for (key,value) in loci_dict
        for (read_name, gene_name) in value.read_names
          if read_name in keys(read_dict)
            push!(read_dict[read_name], gene_name)
          else
            read_dict[read_name] = Set([gene_name])
          end
        end
      end
      println(STDERR,"Processing FASTQ files...")
      for i in 1:num_src_types
        @time extract_X_Q!(parsed_args["fastq"][i], loci_dict, read_dict, parsed_args["format"][i], i, false)
      end
      println(STDERR,"Cleaning up...")
      all_data = cleanup!(loci_dict);

      #@time bin_reads(loci_dict, isoform_dict, chr, parsed_args["temp"], stitched = (parsed_args["phase_psl"] != nothing))
      @time bin_reads(loci_dict, isoform_dict, chr, parsed_args["temp"], stitched = false)
      EXTRA = open(string(parsed_args["out"],"/",parsed_args["prefix"],"_extra_SPECIAL","_$(chr).txt"),"w")
      REGIONS = open(string(parsed_args["out"],"/",parsed_args["prefix"],"_regions","_$(chr).txt"),"w")
      i = 1
      for (gene_name, loci) in loci_dict
        println(EXTRA,string(i,"\t",gene_name,"\t",length(loci.effective_isoform_lengths),"\t",length(loci.effective_region_lengths),"\t",join(loci.regions_in_isoforms,","),"\t",join(loci.Y,","),"\t",join(loci.isoform_map,","),"\t",join(loci.effective_isoform_lengths,","),"\t",join(loci.effective_region_lengths,","),"\t",join(loci.isoform_names,","))) 
        println(REGIONS, string(i, "\t", gene_name,"\t", join(values(loci.regions),","), "\t", join(keys(loci.regions),","),"\t", join(loci.Y,","),"\t",join(loci.snps_in_isoforms,","),"\t",join(loci.exons,","),"\t",join(loci.snp_idx,",")))
        i += 1
      end
      close(EXTRA)
      close(REGIONS)
    end
  else
    println(STDERR, "Getting isoforms with SNPS...")
    @time isoform_dict,loci_dict = initialize_loci("$(parsed_args["temp"])/$(chr)_gpd.bed","$(parsed_args["temp"])/$(chr)_vcf.bed",parsed_args["temp"],chr, num_src_types, fpkm_dict)
    println(STDERR, "Adding reads...")
    @time add_reads!(loci_dict,chr,parsed_args["temp"], num_src_types)

    for (key,value) in loci_dict
      for (read_name, gene_name) in value.read_names
        if read_name in keys(read_dict)
          push!(read_dict[read_name], gene_name)
        else
          read_dict[read_name] = Set([gene_name])
        end
      end
    end

    println(STDERR,"Processing FASTQ files...")
    for i in 1:num_src_types
      @time extract_X_Q!(parsed_args["fastq"][i], loci_dict, read_dict, parsed_args["format"][i], i, false)
    end
    println(STDERR,"Cleaning up...")
    all_data = cleanup!(loci_dict)
    if !parsed_args["only_sim"]
      TRUE = open(string(parsed_args["out"],"/",parsed_args["prefix"],"_true","_$(chr).txt"),"w")
      READS = open(string(parsed_args["out"],"/",parsed_args["prefix"],"_reads","_$(chr).txt"),"w")
      for i in 1:length(all_data)
        loci = loci_dict[all_data[i].gene_name]
        println(TRUE,string(i,"\t","NA","\t",join(all_data[i].true_h,""),"\t",size(all_data[i].X,1),"\t",size(all_data[i].X,2),"\t",join(vec(all_data[i].X),","),"\t",join(vec(all_data[i].Q),","),"\t",all_data[i].gene_name,"\t","NA","\t",join(all_data[i].read_counts,","),"\t",join(loci.snps,","),"\t",join(loci.ref,","),"\t",join(loci.alt,","),"\t",chr))
        println(READS, string(i,"\t",join(all_data[i].read_names,",")))
      end
      close(TRUE)
      close(READS)
    end
    if sim
      sim_data = simulate_data(all_data, string(parsed_args["out"],"/",parsed_args["prefix"],"_sim","_$(chr).txt"))
    end
  end
end
