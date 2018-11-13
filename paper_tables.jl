include("common.jl")

# List the summary statistics of a table
function summary_stats(dataset::AbstractString)
    I, J, times, core = read_data(dataset)
    ncc = 0
    ncf = 0
    for (i, j) in zip(I, J)
        num_in_core = (i in core) + (j in core)
        if num_in_core == 2; ncc += 1; end        
        if num_in_core == 1; ncf += 1; end
    end
    nc = length(core)
    nf = length(unique([I; J])) - nc
    println("$nc & $nf & $ncc & $ncf")
end
;
