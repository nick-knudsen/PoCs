using PyCall
using Random
using StatsBase
using Plots
using GLM

plt = pyimport("matplotlib.pyplot")
display = pyimport("IPython.display")
ndimage = pyimport("scipy.ndimage")
stats = pyimport("scipy.stats")

# probability of a spark at a given site
get_site_prob(row, col, l) = ℯ^(-1*(row/l))*ℯ^(-1*(col/l))

function calc_avg_trees_burned(probs::Array{Float64, 2}, world::Array{Int, 2}, L::Int)
    conn_comps, num_conn_comps = ndimage.label(world,[[0,1,0],[1,1,1],[0,1,0]])
    avg_trees_burned = 0
    for i=1:L
        for j=1:L
            if world[i,j] == 0
                continue
            end
            avg_trees_burned += probs[i,j]*sum(conn_comps .== conn_comps[i, j])
        end
    end
    return avg_trees_burned
end

function calc_avg_yield(world::Array{Int, 2}, probs::Array{Float64, 2}, L::Int)
    num_sites = length(world)
    num_trees = sum(world)
    avg_trees_burned = calc_avg_trees_burned(probs, world, L)

    density = num_trees/num_sites
    cost = avg_trees_burned/num_sites

    return density - cost
end

# populate the world with probabilities
function make_site_probs(L)
    probs = Array{Float64,2}(undef, L, L)
    l = L/10
    for i=1:L
        for j=1:L
            probs[i, j] = get_site_prob(i, j, l)
        end
    end
    return probs ./ sum(probs)
end

function choose_tree(d::Int64, world::Array{Int}, probs::Array{Float64, 2}, L)
    avg_yields = Vector{Float64}()
    open_sites = findall(x -> x==0, world)
    shuffle!(open_sites)
    used_spots = Vector{eltype(open_sites)}()
    forest_sizes = Vector{Int}
    for i=1:d
        if isempty(open_sites)
            continue
        end
        world_copy = deepcopy(world)
        site = pop!(open_sites)
        world_copy[site] = 1
        avg_yield = calc_avg_yield(world_copy, probs, L)

        push!(avg_yields, avg_yield)
        push!(used_spots, site)

    end

    max_yield, max_avg_yield_index = findmax(avg_yields)
    world[used_spots[max_avg_yield_index]] = 1
    
    return max_yield, deepcopy(world)
end

function make_world(L)
    # a matrix of sites
    world = Array{Float64, 2}(undef, L, L)
    probs = make_site_probs(L)
    return world, probs
end

function draw_world(world::Array{Int,2})
    plt.cla()
    plt.imshow(world, cmap="hot")
end

function main()
    
    # dimension
    Ls = [128]
    
    # number of randomly chosen placements of the next tree
    for L in Ls
        D = [1,2,16] 
        probs = make_site_probs(L)
        avg_yields_mat = Array{Float64}(undef, length(D), L^2)
        max_matrices = []

        forest_sizes = Vector{Vector{Int64}}()

        @time for (i, d) in pairs(D)
            println("L = "*string(L)*"  D = "*string(d))
            avg_yields = Vector{Float64}()
            max_yield = 0
            max_matrix = nothing
            world = zeros(Int, L,L)
            tree_num = 0
            while zero(eltype(world)) in world
                avg_yield, sample_world = choose_tree(d, world, probs, L)
                tree_num += 1
                if avg_yield > max_yield
                    max_yield = avg_yield
                    max_matrix = deepcopy(sample_world)
                end
                push!(avg_yields, avg_yield) 
                print("\rGrowing tree "*string(tree_num)*"/"*string(L^2))
            end
            push!(max_matrices, max_matrix)
            avg_yields_mat[i,:] = avg_yields
            conn_comps, num_conn_comps = ndimage.label(max_matrix,[[0,1,0],[1,1,1],[0,1,0]])
            forest_sizes_for_d = Vector{Int64}()
            
            for i=1:num_conn_comps
                forest_size = sum(conn_comps .== i)
                push!(forest_sizes_for_d, forest_size)
                
            end
            
            push!(forest_sizes,forest_sizes_for_d)

            println()
            println("Max Yield: "*string(max_yield))
            println("Density: "*string(sum(max_matrix)/L^2))
        end
        

        folder = pwd()*"/graphs"
        if !isdir(folder)
            mkdir(folder)
        end

        
        for i=1:length(D)
            forest_sizes_for_d = forest_sizes[i]
            sort!(forest_sizes_for_d, rev=true)
            log_ranks = log10.(1:length(forest_sizes_for_d))
            log_forest_sizes = log10.(forest_sizes_for_d)
            slope, intercept, rvalue, pvalue, se = stats.linregress(log_ranks, log_forest_sizes)

            plt.figure()
            fig = plt.scatter(log_ranks, log_forest_sizes)
            plt.plot(log_ranks, intercept .+ slope*log_ranks)
            plt.title("L: "*string(L)*" D: "*string(D[i])*" α: "*string(round(-slope, digits=3)))
            plt.xlabel("log_10 rank")
            plt.ylabel("log_10 forest size")
            plt.savefig(folder*"/zipf_L_"*string(L)*"_d_"*string(D[i])*".png")
            plt.clf()

            fig, ax = plt.subplots(figsize=(10, 10))
            ax.set_ylim(0,L)
            ax.set_xlim(0,L)
            ax.set_aspect("equal")
            draw_world(max_matrices[i])
            plt.title("D="*string(D[i])*" L="*string(L), fontsize=20)
            plt.savefig(folder*"/max_yield_"*"L_"*string(L)*"_d_"*string(D[i])*".png")
            plt.clf()


        end 
        densities = [i/L^2 for i=1:L^2]
        plot(densities, avg_yields_mat', labels=reshape("D=".*string.(D), 1, :), legend=:topleft, fmt = png, size=(700,500))
        title!("Yield Density Curves")
        xlabel!("Density")
        ylabel!("Yield")
        savefig(folder*"/yield_density_"*"L_"*string(L)*"_d_"*".png")
    end
end

main()