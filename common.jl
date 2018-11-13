using LinearAlgebra
using Random
using SparseArrays
using StatsBase

const SpIntMat = SparseMatrixCSC{Int64,Int64}
const Edges = Vector{NTuple{2,Int64}}

struct Comparisons
    edges_formed::Edges
    edges_not_formed::Edges
end

# Get prediction scores
function prediction_scores(A_cf::SpIntMat, comps::Comparisons)
    # Set up score functions
    B = A_cf' * A_cf
    cdegs = vec(sum(A_cf, dims=1))
    common_nbrs_score(w::Int64, z::Int64) = B[w, z]
    jaccard_score(w::Int64, z::Int64) = B[w, z] / (cdegs[w] + cdegs[z] - B[w, z])
    
    # Comparisons based on common neighbors and Jaccard similarity
    com_correct = 0.0
    jac_correct = 0.0
    total = 0
    for (formed_edge, not_formed_edge) in zip(comps.edges_formed, comps.edges_not_formed)
        new_com_score = common_nbrs_score(formed_edge...)
        not_com_score = common_nbrs_score(not_formed_edge...)
        new_jac_score = jaccard_score(formed_edge...)
        not_jac_score = jaccard_score(not_formed_edge...)
        if new_com_score  > not_com_score; com_correct += 1.0; end
        if new_com_score == not_com_score; com_correct += 0.5; end
        if new_jac_score  > not_jac_score; jac_correct += 1.0; end
        if new_jac_score == not_jac_score; jac_correct += 0.5; end
        total += 1
    end
    return (com_correct / total, jac_correct / total, sum(B) - sum(diag(B)))
end

# Get comparisons
function get_comparisons(A_core::SpIntMat, new_edges::Set{NTuple{2,Int64}},
                         core::Set{Int64}, sample_factor::Int64)
    total = length(new_edges) * sample_factor
    new_edge_samples = rand(collect(new_edges), total)
    core_vec = collect(core)
    no_edge_samples = Edges()
    for _ = 1:total
        while true
            i, j = rand(core_vec, 2)
            if i != j && A_core[i, j] == 0
                push!(no_edge_samples, (min(i, j), max(i, j)))
                break
            end
        end
    end
    return Comparisons(new_edge_samples, no_edge_samples)
end

# Make a n x n sparse matrix with all nonzeros equal to 1 and no diagonal
function make_simple_sym_sparse(I::Vector{Int64}, J::Vector{Int64}, n::Int64)
    A = convert(SpIntMat, sparse(I, J, 1, n, n))
    A = A + A'
    A -= Diagonal(A)
    LinearAlgebra.fillstored!(A, 1)
    return A
end

# Make a m x n sparse matrix with all nonzeros equal to 1
function make_simple_sparse(I::Vector{Int64}, J::Vector{Int64}, m::Int64, n::Int64)
    A = convert(SpIntMat, sparse(I, J, 1, m, n))
    LinearAlgebra.fillstored!(A, 1)
    return A
end

# Add a key to a dictionary which maps keys to [1, ..., n]
function add_and_get_val(dict::Dict{T,Int64}, key::T) where T
    if !haskey(dict, key)
        curr = length(dict) + 1
        dict[key] = curr
        return curr
    end
    return dict[key]
end

# Write out value --> key for a Dict
function write_data_map(data_map::Dict, outfile::AbstractString)
    sorted_data_map = sort([(value, key) for (key, value) in data_map])
    open(outfile, "w") do f
        for (value, key) in sorted_data_map
            write(f, "$(value) $(key)\n")
        end
    end
end

function read_core(dataset::AbstractString)
    core = Set{Int64}()
    open("data/$dataset/core-$dataset.txt") do f
        for line in eachline(f); push!(core, parse(Int64, line)); end
    end
    return core
end

function read_data(dataset::AbstractString)
    data = Vector{NTuple{3,Int64}}()
    datasetdir = "data/$dataset"
    open("$(datasetdir)/$dataset.txt") do f
        for line in eachline(f)
            i, j, t = [parse(Int64, x) for x in split(line)]
            if i == j; continue; end
            # undirected edges
            push!(data, (min(i, j), max(i, j), t))
        end
    end
    # remove duplicate in time (only consider first timestamp)
    sort!(data)
    I, J, times = Int64[], Int64[], Int64[]
    curr_i, curr_j = -1, -1
    for (i, j, time) in data
        if i != curr_i || j != curr_j
            curr_i = i
            curr_j = j
            push!(I, curr_i)
            push!(J, curr_j)
            push!(times, time)
        end
    end

    core = read_core(dataset)
    node_map = Dict{Int64, Int64}()
    
    # Core goes first
    core_new = Set{Int64}()
    for c in core; push!(core_new, add_and_get_val(node_map, c)); end
    I_new = [add_and_get_val(node_map, i) for i in I]
    J_new = [add_and_get_val(node_map, j) for j in J]    
    return I_new, J_new, times, core_new
end
