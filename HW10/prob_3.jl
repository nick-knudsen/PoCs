using PyCall
ndimage = PyCall("scipy.ndimage")

# probability of a spark at a given site
get_site_prob(row, col) = ℯ^(row/l)*ℯ^(col/l)

function calc_avg_trees_burned(probs::Array{Float64, 2}, world:Array{Int, 2}, L::int)
    conn_comps, num_conn_comps = ndimage.label(world,[[0,1,0],[1,1,1],[0,1,0]])
    avg_trees_burned = 0
    for i=1:L
        for j=1:L
            if world(i, j) == 0
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
    cost = num_trees_burned/num_sites

    return density - cost
end

# populate the world with probabilities
function make_site_probs(L)
    probs = Array{Float64,2}(undef, L, L)
    for i=1:L
        for j=1:L
            probs[i, j] = get_site_prob(i, j)
        end
    end
    return probs ./ sum(probs)
end

function choose_tree(d::Int64, world::Array{Int}, probs::Array{Int})
    avg_yields = Vector{Float64}()
    open_sites = findall(x -> x==0, world)
    shuffle!(open_sites)
    all_open_spots = deepcopy(open_spots)

    for i=1:d
        world_copy = deepcopy(world)
        site = pop!(open_sites)
        world_copy[site] = 1
        avg_yield = calc_avg_yield(world_copy, probs, L)

        push!(avg_yield, avg_yields)
    end

    max_avg_yield_index = findmax(avg_yields)
    world[a]
end

function make_world(L)
    # a matrix of sites
    world = Array{Float64, 2}(undef, L, L)
    probs = make_site_probs(L)
    return world, probs
end

# dimension
L = 32
l= L/10
# number of randomly chosen placements of the next tree
D = [1, 2, L, L^2]

for d in D
    world = zeros(Int, L,L)
    while zero(eltype(world)) in world
        choose_tree(d, world)
    end

end
