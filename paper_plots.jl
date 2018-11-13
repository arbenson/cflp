using PyPlot
using StatsBase

function sbm_snr(p::Float64, q::Float64, r::Float64, s::Float64,
                 nc::Int64, nf::Int64)
    # expected number of common neighbors from core
    pq = p * q
    pp = p * p
    rr = r * r
    ss = s * s
    rs = r * s

    # Same
    exp_same, var_same = 0.0, 0.0
    # block 1 --> block 1
    exp_same += (nc - 2) * pp
    var_same += (nc - 2) * pp * (1 - pp)
    # block 2 --> block 1
    exp_same += nc * pp
    var_same += nc * pp * (1 - pp)
    # block 3 --> block 1
    exp_same += nf * rr
    var_same += nf * rr * (1 - rr)
    # block 4 --> block 1
    exp_same += nf * ss
    var_same += nf * ss * (1 - ss)
    
    # Diff
    exp_diff, var_diff = 0.0, 0.0
    # block 1 + 2
    exp_diff += 2 * (nc - 1) * pq
    var_diff += 2 * (nc - 1) * pq * (1 - pq)
    # block 3 + 4
    exp_diff += 2 * nf * rs
    var_diff += 2 * nf * rs * (1 - rs)

    diff = (exp_same - exp_diff)
    var = var_same + var_diff
    return diff / sqrt(var)
end

function sbm_plot(p::Float64, q::Float64, r::Float64, ss::Vector{Float64}, nc::Int64, nfs::UnitRange{Int64},
                  outfilename::AbstractString)
    close()
    linestyles = ["-", "--", ":"]
    for (ind, s) in enumerate(ss)
        snrs = [sbm_snr(p, q, r, s, nc, nf) for nf in nfs]
        plot(collect(nfs), snrs, label="s=$s", linestyle=linestyles[(ind - 1) % 3 + 1])
    end
    fsz = 22
    ax = gca()
    ax[:tick_params]("both", labelsize=fsz-5, length=5, width=2)
    legend(fontsize=fsz-2)
    xlabel("Fringe size (d)", fontsize=fsz)
    ylabel(L"SNR(Z_d)", fontsize=fsz)
    title("SBM (p=$p, q=$q, r=$r)", fontsize=fsz)
    savefig(outfilename, bbox_inches="tight")
    show()
end
# sbm_plot(0.5, 0.3, 0.2, [0.0, 0.2], 10, 0:100, "SBM1.pdf")
# sbm_plot(0.5, 0.3, 0.2, [0.1], 10, 0:300, "SBM2.pdf")

function lattice_snr(a::Int64, b::Int64, c::Int64, d::Int64)
    avec = Float64[]
    bvec = Float64[]
    u = -c
    v = u + a
    z = c
    w = z - b
    for s = -(c + d):(c + d)
        if s != u && s != v; push!(avec, 1 / abs(s - u) / abs(s - v)); end
        if s != w && s != z; push!(bvec, 1 / abs(s - w) / abs(s - z)); end
    end
    exp = sum(avec) - sum(bvec)
    var = sum(avec) - sum(avec .^ 2) + sum(bvec) - sum(bvec .^ 2)
    return exp / sqrt(var)
end

function lattice_plot(as::UnitRange{Int64}, b::Int64, c::Int64, ds::UnitRange{Int64})
    close()
    inds = collect(ds)
    for a in as
        snrs = [lattice_snr(a, b, c, d) for d in inds]
        maxval, ind = findmax(snrs)
        plot(inds, snrs, marker=".", markersize=12, label="v-u = $a", lw=3)
        plot(inds[ind], maxval, marker="o", markersize=15, color="black", fillstyle="none")
    end
    fsz = 22
    ax = gca()
    ax[:tick_params]("both", labelsize=fsz-5, length=5, width=2)
    legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.,
           frameon=false, fontsize=fsz-2)
    xlabel("Fringe size (d)", fontsize=fsz)
    ylabel(L"SNR(Z_d)", fontsize=fsz)    
    title("Lattice model (c = $c, z-w = $b)", fontsize=fsz)

    savefig("lattice-$b-$c.pdf", bbox_inches="tight")
    show()
end
# lattice_plot(2:6, 9, 10, 0:12)

function μ_σ_scores(scores::Vector{Vector{Float64}})
    μs = Float64[]
    σs = Float64[]
    for i = 1:length(scores[1])
        vals = [score_vec[i] for score_vec in scores]
        μ, σ = mean_and_std(vals)
        push!(μs, μ)
        push!(σs, σ)
    end
    return μs, σs
end

function read_scores(datafile::AbstractString, score_func_ind::Int64)
    scores = Vector{Float64}()
    open(datafile) do f
        for line in eachline(f)
            push!(scores, parse(Float64, split(line)[score_func_ind]))
        end
    end
    return scores
end

function empirical_plot(dataset::AbstractString, score_func::AbstractString, ntrials::Int64=10)
    close()
    conn_scores, rand_scores = Vector{Vector{Float64}}(), Vector{Vector{Float64}}()
    score_func_ind = 1
    if score_func == "Jaccard"; score_func_ind = 2; end
    for trial = 1:ntrials
        push!(conn_scores, read_scores("output/$dataset/scores-$dataset-$trial.txt", score_func_ind))
        push!(rand_scores, read_scores("output/$dataset/rand-scores-$dataset-$trial.txt", score_func_ind))
    end
    
    μ_conn_scores = μ_σ_scores(conn_scores)[1]
    μ_rand_scores = μ_σ_scores(rand_scores)[1]
    inds = vec([0.8; collect(1:(length(μ_conn_scores)-1))])
    colors = ["#e41a1c", "#377eb8"]
    plot(inds, μ_conn_scores, color=colors[1], lw=2.5, ls="-",  label="most conn.")
    plot(inds, μ_rand_scores, color=colors[2], lw=1.5, ls="--", label="random")
    for (ind, scores) in enumerate([μ_conn_scores, μ_rand_scores])
        maxval = maximum(scores)
        maxind = inds[minimum([i for (i, s) in enumerate(scores) if s >= maxval])]
        plot([maxind], [maxval], marker="o", markersize=15, color=colors[ind], fillstyle="none")
    end
    fsz = 22
    ax = gca()
    ax[:set_xscale]("log")
    ax[:tick_params]("both", labelsize=fsz-5, length=5, width=1)
    legend(fontsize=fsz-2)
    xlabel("Number of fringe nodes (d)", fontsize=fsz)
    ylabel("Fraction correct",           fontsize=fsz)
    title("$dataset ($score_func)",      fontsize=fsz)
    savefig(replace("$dataset-$(score_func).pdf", " " => "_"), bbox_inches="tight")
    show()
end
;
