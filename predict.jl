include("common.jl")

function make_predictions(dataset::AbstractString, trial::Int64, rand_ordering::Bool)
    print(stdout, "trial = $trial\n")
    print(stdout, "random ordering = $rand_ordering\n")
    Random.seed!(1234 * trial)  # for reproducibility
    
    I, J, times, core = read_data(dataset)
    core_times = Int64[]
    I_core, J_core = Int64[], Int64[]
    for (i, j, t) in zip(I, J, times)
        if (i in core) && (j in core)
            push!(core_times, t)
            push!(I_core, i)
            push!(J_core, j)
        end
    end
    cutoff_time = percentile(core_times, 80)
    ncore = maximum(core)

    # Get the ordering
    adjacencies = Vector{Set{Int64}}()
    new_edges = Set{NTuple{2,Int64}}()
    nnodes = max(maximum(I), maximum(J))
    for i in 1:nnodes; push!(adjacencies, Set{Int64}()); end
    for (i, j, t) in zip(I, J, times)
        if t <= cutoff_time
            i_in = i in core
            j_in = j in core
            if i_in + j_in == 1
                if !i_in; push!(adjacencies[i], j); end
                if !j_in; push!(adjacencies[j], i); end
            end
        else
            if (i in core) && (j in core)
                push!(new_edges, (i, j))
            end
        end
    end

    A_core = make_simple_sym_sparse(I_core, J_core, ncore)
    comparisons = get_comparisons(A_core, new_edges, core, 10)
    counts = [length(x) for x in adjacencies]
    sp = sortperm(counts, rev=true)
    ordering = Int64[]
    for (node, count) in zip(sp, counts[sp])
        if count == 0; break; end
        push!(ordering, node)
    end
    if rand_ordering; shuffle!(ordering); end
    
    active_nodes = zeros(Int64, nnodes)
    for c in core; active_nodes[c] = 1; end
    common_scores    = Float64[]
    jaccard_scores   = Float64[]
    total_num_common = Int64[]
    for ind = 1:(length(ordering) + 1)
        I_curr, J_curr = Int64[], Int64[]
        for (i, j, time) in zip(I, J, times)
            if time > cutoff_time; continue; end
            if active_nodes[i] == 1 && active_nodes[j] == 1
                if (i, j) in new_edges || (j, i) in new_edges; continue; end                
                if i in core
                    push!(I_curr, j)
                    push!(J_curr, i)
                end
                if j in core
                    push!(I_curr, i)
                    push!(J_curr, j)
                end
            end
        end
        A_curr = make_simple_sparse(I_curr, J_curr, maximum(I_curr), ncore)
        cn_score, j_score, num_common = prediction_scores(A_curr, comparisons)
        push!(common_scores,  cn_score)
        push!(jaccard_scores, j_score)
        push!(total_num_common, num_common)
        print(stdout, "$ind of $(length(ordering)) \r")
        if ind <= length(ordering)
            active_nodes[ordering[ind]] = 1
        end
    end
    print(stdout, "\n")
    flush(stdout)

    outfilename = "scores-$dataset-$trial.txt"
    outdir = "output/$dataset"
    if !isdir(outdir); mkdir(outdir); end
    if rand_ordering; outfilename = "rand-$(outfilename)"; end
    open("$outdir/$outfilename", "w") do f
        for (cs, js, nc) in zip(common_scores, jaccard_scores, total_num_common)
            write(f, "$(cs) $(js) $(nc)\n")
        end
    end
end
;
