include("common.jl")

# Function used to split a dataset into pieces (used for Enron)
function split_dataset(dataset::AbstractString, nsplit::Int64)
    core_vec = collect(read_core(dataset))
    Random.seed!(78910)  # for reproducibility
    shuffle!(core_vec)
    ncore = length(core_vec)
    chunk_size = convert(Int64, floor(ncore / nsplit))

    for chunk= 1:convert(Int64, ceil(ncore / chunk_size))
        start_ind = (chunk - 1) * chunk_size + 1
        end_ind = min(chunk * chunk_size, ncore)
        chunk_core = Set{Int64}(core_vec[start_ind:end_ind])
        chunk_data = Vector{NTuple{3,Int64}}()        
        open("data/$dataset/$dataset.txt") do f
            for line in eachline(f)
                i, j, t = [parse(Int64, x) for x in split(line)]
                if i in chunk_core || j in chunk_core
                    push!(chunk_data, (i, j, t))
                end
            end
        end
        datadir = "data/$dataset-$chunk"
        if !isdir(datadir); mkdir(datadir); end
        open("$datadir/core-$dataset-$chunk.txt", "w") do f
            for c in chunk_core; write(f, "$c\n"); end
        end
        open("$datadir/$dataset-$chunk.txt", "w") do f
            for (i, j, t) in chunk_data; write(f, "$i $j $t\n"); end
        end
    end
end
#split("email-Enron", 4)
