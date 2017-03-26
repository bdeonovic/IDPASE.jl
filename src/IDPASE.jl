module IDPASE

  using Distributions
  using StatsBase
  using DataStructures
  using Mamba
  using Graphs
  using KernelDensity
  using Optim
  using Combinatorics
  
  immutable PSLEntry{V <: AbstractString, W <: Integer}
    matches::W
    mismatches::W
    repmatches::W
    n_count::W
    q_num_insert::W
    q_base_insert::W
    t_num_insert::W
    t_base_insert::W
    strand::V
    q_name::V
    q_size::W
    q_start::W
    q_end::W
    t_name::V
    t_size::W
    t_start::W
    t_end::W
    block_count::W
    block_sizes::Vector{W}
    q_starts::Vector{W}
    t_starts::Vector{W}
  end
  
  immutable GPDEntry{V <: AbstractString, W <: Integer}
    bin::W
    name::V
    chr::V
    strand::V
    tx_start::W
    tx_end::W
    cds_start::W
    cds_end::W
    exon_count::W
    exon_starts::Vector{W}
    exon_ends::Vector{W}
    score::W
    gene_name::V
    cds_start_stat::V
    cds_end_stat::V
    exon_frames::Vector{W}
  end
  type IsoformPhaseEntry{U <:AbstractString, V <: Integer}
    chr::U
    snps::Vector{V}
    ref::Vector{U}
    alt::Vector{U}

    snp_names::Dict{U,Bool}
    read_names::Dict{U,U}
    true_h::Vector{V}
    est_h::Vector{V}
    rho_mode::Float64
    phased::Bool

    reads::Vector{PSLEntry{U,V}}
    read_names_order::Vector{U}
    isoform_names::Vector{U}
    X::Matrix{V}
    Q::Matrix{Float64}

    read_counts::Vector{V}
    exons::Matrix{V}
  
    snps_in_isoforms::Matrix{V}
    exons_in_isoforms::Matrix{V}
    regions_in_isoforms::Matrix{V}

    effective_region_lengths::Vector{V}
    effective_isoform_lengths::Vector{V}
    regions::DataStructures.OrderedDict{Set{V},V}
    snp_idx::Vector{V}
    isoform_map::Vector{V}
    Y::Vector{V}
  end
  type LociEntry{U <: AbstractString, V <: Integer}
    chr::U
    snps::Vector{V}
    ref::Vector{U}
    alt::Vector{U}

    snp_names::Dict{U,Bool}
    read_names::Dict{U,Bool}
    true_h::Vector{V}
    phased::Bool

    exons::Matrix{V}
  
    exons_in_isoforms::Matrix{V}
    regions_in_isoforms::Matrix{V}

    effective_region_lengths::Vector{V}
    effective_isoform_lengths::Vector{V}
  
    Y1::Vector{V}
    X1::Matrix{V}
  
    Y2::Vector{V}
    X2::Matrix{V}

    X3::Matrix{V}
  end
  type PhaseEntry{U <: AbstractString, V <: Integer}
    chr::U
    snps::Vector{V}
    ref::Vector{U}
    alt::Vector{U}

    snp_names::Dict{U,Bool}
    read_names::Dict{U,U}
    true_h::Vector{V}
    phased::Bool

    reads::Vector{PSLEntry{U,V}}
    read_names_order::Vector{U}
    X::Matrix{V}
    Q::Matrix{Float64}

    read_counts::Vector{V}
  end


  abstract PHASEData
  type PHASEDataReal{T <: AbstractString, U <: AbstractFloat, W <: Integer} <: PHASEData
    read_counts::Vector{W}
    coords::Vector{W}
    gene_name::T
    isoform_name::T
    chr::T
    read_names::Vector{T}

    reference_alleles::Vector{T}
    alternative_alleles::Vector{T}

    X::Matrix{W}
    Q::Matrix{U}

    gamma::Vector{Vector{W}}
    true_h::Vector{W}
  end

  type PHASEDataSim{T <: AbstractString, U <: AbstractFloat, W <: Integer} <: PHASEData
    read_counts::Vector{W}
    coords::Vector{W}
    gene_name::T
    isoform_name::T
    chr::T
    read_names::Vector{T}

    reference_alleles::Vector{T}
    alternative_alleles::Vector{T}

    X::Matrix{W}
    Q::Matrix{U}

    gamma::Vector{Vector{W}}

    true_h::Vector{W}
    rho::U
    delta::Vector{W}

    rel_abundance::Vector{U}
  end
  coords(d::PHASEData) = d.coords
  isoform_lengths(d::PHASEData) = d.isoform_lengths
  isoforms(d::PHASEData) = d.isoforms
  reference_alleles(d::PHASEData) = d.reference_alleles
  alternative_alleles(d::PHASEData) = d.alternative_alleles
  X(d::PHASEData) = d.X
  Q(d::PHASEData) = d.Q
  gamma(d::PHASEData) = d.gamma

  const DNA_COMPLEMENT = zeros(Integer,256)
  DNA_COMPLEMENT[convert(Int64,'.')] = '.'
  DNA_COMPLEMENT[convert(Int64,'-')] = '-'
  DNA_COMPLEMENT[convert(Int64,'A')] = 'T'
  DNA_COMPLEMENT[convert(Int64,'B')] = 'V'
  DNA_COMPLEMENT[convert(Int64,'C')] = 'G'
  DNA_COMPLEMENT[convert(Int64,'D')] = 'H'
  DNA_COMPLEMENT[convert(Int64,'G')] = 'C'
  DNA_COMPLEMENT[convert(Int64,'H')] = 'D'
  DNA_COMPLEMENT[convert(Int64,'K')] = 'M'
  DNA_COMPLEMENT[convert(Int64,'M')] = 'K'
  DNA_COMPLEMENT[convert(Int64,'T')] = 'A'
  DNA_COMPLEMENT[convert(Int64,'V')] = 'B'
  DNA_COMPLEMENT[convert(Int64,'S')] = 'S'
  DNA_COMPLEMENT[convert(Int64,'W')] = 'W'
  DNA_COMPLEMENT[convert(Int64,'R')] = 'Y'
  DNA_COMPLEMENT[convert(Int64,'Y')] = 'R'
  DNA_COMPLEMENT[convert(Int64,'N')] = 'N'
  DNA_COMPLEMENT[convert(Int64,'X')] = 'X'
  DNA_COMPLEMENT[convert(Int64,'a')] = 't'
  DNA_COMPLEMENT[convert(Int64,'b')] = 'v'
  DNA_COMPLEMENT[convert(Int64,'c')] = 'g'
  DNA_COMPLEMENT[convert(Int64,'d')] = 'h'
  DNA_COMPLEMENT[convert(Int64,'g')] = 'c'
  DNA_COMPLEMENT[convert(Int64,'h')] = 'd'
  DNA_COMPLEMENT[convert(Int64,'k')] = 'm'
  DNA_COMPLEMENT[convert(Int64,'m')] = 'k'
  DNA_COMPLEMENT[convert(Int64,'t')] = 'a'
  DNA_COMPLEMENT[convert(Int64,'v')] = 'b'
  DNA_COMPLEMENT[convert(Int64,'s')] = 's'
  DNA_COMPLEMENT[convert(Int64,'w')] = 'w'
  DNA_COMPLEMENT[convert(Int64,'r')] = 'y'
  DNA_COMPLEMENT[convert(Int64,'y')] = 'r'
  DNA_COMPLEMENT[convert(Int64,'x')] = 'x'
  DNA_COMPLEMENT[convert(Int64,'n')] = 'n'

  function reversecomplement!(s::AbstractArray,y::Array)
    i=1
    j=length(s)
    while i<=j
      I=convert(Int64,s[i])
      J=convert(Int64,s[j])
      s[i] = convert(Char,y[J])
      s[j] = convert(Char,y[I])
      i += 1
      j -= 1
    end
    s
  end
  reversecomplement!(s::String,y::Array) = reversecomplement!(Char[s[i] for i in 1:length(s)],y)

  reversecomplement(s::AbstractArray,y::Array) = reversecomplement!(identity(s),y)
  reversecomplement(s::String,y::Array) = reversecomplement!(identity(s),y)


  include("utilities.jl")
  include("isoform_utilities.jl")
  include("chisqdiag.jl")

  export
  PSLEntry,
  GPDEntry,
  LociEntry,
  PhaseEntry,
  IsoformPhaseEntry,
  DNA_COMPLEMENT

  export
    bin_reads,
    haplo_lik,
    phase_isoform,
    add_phase_reads!,
    assemble_loci!,
    read_data,
    make_X_Q,
    get_snps!,
    get_gene_level_results,
    init_loci,
    reversecomplement,
    reversecomplement!,
    psl_to_bed,
    psl_to_bed2,
    process_gpd,
    process_gpd_with_vcf,
    extract_X_Q!,
    cleanup!,
    gpd_to_bed,
    vcf_to_bed,
    get_psl_entries_on_snps,
    get_psl_entries_on_genes,
    initialize_loci,
    initialize_loci_with_vcf,
    exon_composition,
    get_regions,
    run_model,
    run_model2,
    run_model3,
    run_model4,
    run_model5,
    run_model6,
    run_model7,
    run_model10,
    run_model11,
    rho_conditional_mean,
    phase,
    add_reads!,
    add_short_reads!,
    add_long_reads!,
    run_MCMC,
    simulate_data
end # module
