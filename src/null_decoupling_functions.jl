#--------------------------------------------------------------#
#-Null models I, I-A, II and II-B: decoupling the two layers---#
#--------------------------------------------------------------#

# These null models test whether the coupling of pollination and seed dispersal
# by dual-role birds shapes the extinctions that hurricanes cause in the island
# networks (Extended Data Fig. 1). Every null network decouples the pollination
# and seed-dispersal layers of an empirical network while holding its remaining
# topological properties fixed: species richness, the number of interactions in
# each layer, and the identity and degree of the plants in each layer. Two
# classes of null models decouple the layers.
#
# Class I decouples the layers at the level of bird species. Each decoupled
# connector (a bird that both pollinates and disperses seeds) retains all of
# its interactions in one layer, and its interactions in the other layer are
# rewired to single-role birds of that layer.
#   Null model I.    All connectors are decoupled, so the layers separate.
#   Null model I-A.  Half of the connectors, drawn at random, are decoupled.
#
# Class II decouples the layers at the level of interactions. Interactions of
# dual-role birds are drawn one at a time, regardless of which bird holds them,
# and rewired to single-role birds of the same layer.
#   Null model II.   Rewiring continues until no bird interacts in both layers.
#   Null model II-B. Rewiring stops once dual-role birds retain half of the
#                    interactions they held in the empirical network.
#
# In all null models, a rewired interaction keeps its plant and its layer, and
# the receiving bird (1) interacts exclusively in that layer, so that no new
# connector appears; (2) does not yet interact with that plant; and (3) is not
# a hummingbird receiving a fruit. In a few cases, a plant already interacts
# with every bird that could receive the interaction, which is then rewired to
# another plant of the same layer.
#
# Inputs:  the empirical layers and the hummingbird list, read with the
#          helpers of src/null_io_functions.jl (included before this file)
# Outputs: the null matrices and the manifest of src/null_<model>/, written
#          by build_null_ensemble
# Helpers: src/null_decoupling_checks.jl (included after this file) checks
#          every null network before it is written
# Used by: scripts/01_network_preparation/build_null_decoupling.jl, which
#          includes the three files in this order

using Random

#--------------------------------------------------------------#
#-Layers, dual-role birds and receiving birds------------------#
#--------------------------------------------------------------#

# The two layers of a network are a vector of binary incidence matrices
# (plants x birds) sharing the same rows and columns: layers[POLLINATION]
# holds the flower visits and layers[DISPERSAL] the fruit consumption.
const POLLINATION = 1
const DISPERSAL = 2
const LAYER_NAMES = ["pollination", "dispersal"]

# The layer in which a connector's interactions are rewired once it retains
# the other one. Arguments: retained_layer is POLLINATION or DISPERSAL.
rewired_layer(; retained_layer::Int64) = 3 - retained_layer

# The layer of each bird: POLLINATION or DISPERSAL for single-role birds, and
# 0 for dual-role birds, which interact in both layers.
# Arguments: layers holds the two matrices. Returns one Int64 per bird.
function bird_layer(; layers::Vector{Matrix{Int64}})
    in_pol = vec(sum(layers[POLLINATION]; dims=1)) .> 0
    in_disp = vec(sum(layers[DISPERSAL]; dims=1)) .> 0
    layer = zeros(Int64, length(in_pol))
    layer[in_pol .& .!in_disp] .= POLLINATION
    layer[in_disp .& .!in_pol] .= DISPERSAL
    return layer
end

# The dual-role birds of a network: the connectors in the empirical network.
# Arguments: layers holds the two matrices. Returns their column indices.
function dual_role_birds(; layers::Vector{Matrix{Int64}})
    in_pol = vec(sum(layers[POLLINATION]; dims=1)) .> 0
    in_disp = vec(sum(layers[DISPERSAL]; dims=1)) .> 0
    return findall(in_pol .& in_disp)
end

# Number of interactions held by dual-role birds, summed over both layers.
# Class II reduces this number. Arguments: layers holds the two matrices.
function count_dual_role_interactions(; layers::Vector{Matrix{Int64}})
    dual = dual_role_birds(layers=layers)
    return sum(layers[POLLINATION][:, dual]) + sum(layers[DISPERSAL][:, dual])
end

# Every interaction held by a dual-role bird, as (layer, plant, bird).
# Arguments: layers holds the two matrices. Returns a vector of tuples.
function dual_role_interactions(; layers::Vector{Matrix{Int64}})
    dual = dual_role_birds(layers=layers)
    n_plants = size(layers[POLLINATION], 1)
    return [(layer, plant, bird)
            for layer in (POLLINATION, DISPERSAL), plant in 1:n_plants, bird in dual
            if layers[layer][plant, bird] == 1]
end

# Whether a bird can receive the interaction of a plant in a layer. The bird
# must interact exclusively in that layer, must not yet interact with that
# plant, and must not be a hummingbird receiving a fruit.
# Arguments: layers holds the current matrices; layer, plant and bird locate
# the candidate interaction; receiver_layer gives the layer in which each bird
# can receive interactions (0 for none); hummingbirds marks the hummingbirds.
function can_receive(; layers::Vector{Matrix{Int64}}, layer::Int64, plant::Int64,
                       bird::Int64, receiver_layer::Vector{Int64},
                       hummingbirds::Vector{Bool})
    return receiver_layer[bird] == layer &&
           layers[layer][plant, bird] == 0 &&
           !(layer == DISPERSAL && hummingbirds[bird])
end

# The birds that can receive the interaction of a plant in a layer.
# Arguments: as in can_receive, without bird. Returns their column indices.
function receiving_birds(; layers::Vector{Matrix{Int64}}, layer::Int64, plant::Int64,
                           receiver_layer::Vector{Int64}, hummingbirds::Vector{Bool})
    return filter(bird -> can_receive(layers=layers, layer=layer, plant=plant, bird=bird,
                                      receiver_layer=receiver_layer,
                                      hummingbirds=hummingbirds),
                  1:size(layers[layer], 2))
end

# The (plant, bird) cells of other plants of the same layer that can receive
# an interaction. Only plants that already belong to the layer in the
# empirical network can receive, so that no plant changes layers.
# Arguments: as in receiving_birds; plant is the plant to exclude;
# plants_in_layer marks the plants of the layer in the empirical network.
function receiving_cells_at_other_plants(; layers::Vector{Matrix{Int64}}, layer::Int64,
                                           plant::Int64, plants_in_layer::BitVector,
                                           receiver_layer::Vector{Int64},
                                           hummingbirds::Vector{Bool})
    n_plants, n_birds = size(layers[layer])
    return [(other, bird) for other in 1:n_plants, bird in 1:n_birds
            if other != plant && plants_in_layer[other] &&
               can_receive(layers=layers, layer=layer, plant=other, bird=bird,
                           receiver_layer=receiver_layer, hummingbirds=hummingbirds)]
end

#--------------------------------------------------------------#
#-Rewiring one interaction-------------------------------------#
#--------------------------------------------------------------#

# Rewire the interaction between plant and bird to a randomly chosen receiving
# bird at the same plant. The plant keeps its degree in the layer.
# Arguments: as in receiving_birds, plus bird (the current holder) and rng.
function rewire_at_same_plant!(; layers::Vector{Matrix{Int64}}, layer::Int64,
                                 plant::Int64, bird::Int64, receiver_layer::Vector{Int64},
                                 hummingbirds::Vector{Bool}, rng::Random.Xoshiro)
    receivers = receiving_birds(layers=layers, layer=layer, plant=plant,
                                receiver_layer=receiver_layer, hummingbirds=hummingbirds)
    layers[layer][plant, bird] = 0
    layers[layer][plant, rand(rng, receivers)] = 1
    return nothing
end

# Rewire the interaction between plant and bird to a randomly chosen receiving
# bird at another plant of the same layer. Two plants change their degree by
# one interaction. Arguments: as in receiving_cells_at_other_plants, plus bird
# (the current holder) and rng.
function rewire_to_other_plant!(; layers::Vector{Matrix{Int64}}, layer::Int64,
                                  plant::Int64, bird::Int64, plants_in_layer::BitVector,
                                  receiver_layer::Vector{Int64}, hummingbirds::Vector{Bool},
                                  rng::Random.Xoshiro)
    cells = receiving_cells_at_other_plants(layers=layers, layer=layer, plant=plant,
                                            plants_in_layer=plants_in_layer,
                                            receiver_layer=receiver_layer,
                                            hummingbirds=hummingbirds)
    isempty(cells) &&
        error("No bird can receive a $(LAYER_NAMES[layer]) interaction at another plant.")
    layers[layer][plant, bird] = 0
    other_plant, other_bird = rand(rng, cells)
    layers[layer][other_plant, other_bird] = 1
    return nothing
end

#--------------------------------------------------------------#
#-Class I: decoupling at the level of bird species-------------#
#--------------------------------------------------------------#

# Count the interactions that would be rewired to another plant under one
# assignment of retained layers, and whether the assignment is feasible. At
# each plant, a layer has as many interactions to rewire as decoupled
# connectors that give it up there, and as many receiving birds as birds of
# that layer that do not yet interact with the plant. Interactions in excess
# have no receiver at their plant. Receivers in excess are the places where
# such interactions can go, so each layer needs at least as many of them.
# Arguments: layers holds the empirical matrices; decoupled lists the
# connectors to decouple; retained_layers gives the layer each one retains;
# hummingbirds marks the hummingbirds. Returns the count and the feasibility.
function count_rewired_to_other_plants(; layers::Vector{Matrix{Int64}},
                                         decoupled::Vector{Int64},
                                         retained_layers::Vector{Int64},
                                         hummingbirds::Vector{Bool})
    receiver_layer = bird_layer(layers=layers)
    receiver_layer[decoupled] .= retained_layers
    n_without_receiver = zeros(Int64, 2)   # per layer
    n_spare_receivers = zeros(Int64, 2)    # per layer

    for layer in (POLLINATION, DISPERSAL), plant in 1:size(layers[layer], 1)
        any(==(1), layers[layer][plant, :]) || continue   # plant not in this layer
        n_to_rewire = count(i -> retained_layers[i] != layer &&
                                 layers[layer][plant, decoupled[i]] == 1,
                            eachindex(decoupled))
        n_receivers = length(receiving_birds(layers=layers, layer=layer, plant=plant,
                                             receiver_layer=receiver_layer,
                                             hummingbirds=hummingbirds))
        n_without_receiver[layer] += max(n_to_rewire - n_receivers, 0)
        n_spare_receivers[layer] += max(n_receivers - n_to_rewire, 0)
    end

    return (n=sum(n_without_receiver), feasible=all(n_without_receiver .<= n_spare_receivers))
end

# The allowed assignments of retained layers for a set of decoupled connectors:
# among all combinations (each connector retains pollination or dispersal),
# the feasible ones that rewire the fewest interactions to another plant.
# Arguments: as in count_rewired_to_other_plants, without retained_layers.
# Returns a vector of assignments, each a vector of retained layers in the
# order of decoupled. The vector is empty when no combination is feasible.
function allowed_assignments(; layers::Vector{Matrix{Int64}}, decoupled::Vector{Int64},
                               hummingbirds::Vector{Bool})
    combinations = vec([collect(Int64, retained) for retained in
                        Iterators.product(fill((POLLINATION, DISPERSAL), length(decoupled))...)])
    costs = [count_rewired_to_other_plants(layers=layers, decoupled=decoupled,
                                           retained_layers=retained, hummingbirds=hummingbirds)
             for retained in combinations]
    feasible = [cost.feasible for cost in costs]
    any(feasible) || return Vector{Int64}[]
    fewest = minimum(cost.n for cost in costs[feasible])
    return combinations[feasible .& [cost.n == fewest for cost in costs]]
end

# Draw the connectors that null model I-A decouples: half of them, or, with an
# odd number of connectors, the number just below or above half with equal
# probability. A drawn set without any feasible assignment is redrawn.
# Arguments: layers holds the empirical matrices; connectors lists all
# connectors; share is the share of connectors to decouple; hummingbirds marks
# the hummingbirds; cache stores the allowed assignments of the sets already
# drawn; rng is the generator of the replicate.
# Returns the decoupled connectors and their allowed assignments.
function draw_decoupled_connectors(; layers::Vector{Matrix{Int64}},
                                     connectors::Vector{Int64}, share::Float64,
                                     hummingbirds::Vector{Bool},
                                     cache::Dict{Vector{Int64}, Vector{Vector{Int64}}},
                                     rng::Random.Xoshiro)
    expected = share * length(connectors)
    below = floor(Int64, expected)
    n_decoupled = rand(rng) < expected - below ? below + 1 : below

    for attempt in 1:1000
        decoupled = sort(shuffle(rng, connectors)[1:n_decoupled])
        allowed = get!(cache, decoupled) do
            allowed_assignments(layers=layers, decoupled=decoupled, hummingbirds=hummingbirds)
        end
        isempty(allowed) || return (decoupled=decoupled, allowed=allowed)
    end
    error("No set of $n_decoupled connectors can be decoupled.")
end

# Build one null network of class I. Each decoupled connector retains all of
# its interactions in its retained layer. Its interactions in the other layer
# are rewired, one by one, to randomly chosen receiving birds at the same
# plant. The interactions without a receiving bird at their plant are then
# rewired to another plant of the same layer. The connectors left untouched
# keep both layers and never receive interactions.
# Arguments: empirical_layers holds the matrices to decouple; decoupled lists
# the connectors to decouple and retained_layers the layer each one retains;
# hummingbirds marks the hummingbirds; rng is the generator of the replicate.
# Returns the null layers and the number of interactions rewired to another
# plant.
function build_class_I_network(; empirical_layers::Vector{Matrix{Int64}},
                                 decoupled::Vector{Int64}, retained_layers::Vector{Int64},
                                 hummingbirds::Vector{Bool}, rng::Random.Xoshiro)
    layers = [copy(empirical_layers[POLLINATION]), copy(empirical_layers[DISPERSAL])]
    plants_in_layer = [vec(sum(empirical_layers[layer]; dims=2)) .> 0
                       for layer in (POLLINATION, DISPERSAL)]

    # Single-role birds receive in their layer; decoupled connectors in the
    # layer they retain; untouched connectors never receive
    receiver_layer = bird_layer(layers=empirical_layers)
    receiver_layer[decoupled] .= retained_layers

    # Rewire the interactions of each decoupled connector in its other layer to
    # receiving birds at the same plant. Interactions without a receiving bird
    # wait for the second pass.
    without_receiver = Tuple{Int64, Int64, Int64}[]   # (layer, plant, bird)
    for (i, connector) in enumerate(decoupled)
        layer = rewired_layer(retained_layer=retained_layers[i])
        for plant in findall(==(1), empirical_layers[layer][:, connector])
            receivers = receiving_birds(layers=layers, layer=layer, plant=plant,
                                        receiver_layer=receiver_layer,
                                        hummingbirds=hummingbirds)
            if isempty(receivers)
                push!(without_receiver, (layer, plant, connector))
            else
                rewire_at_same_plant!(layers=layers, layer=layer, plant=plant,
                                      bird=connector, receiver_layer=receiver_layer,
                                      hummingbirds=hummingbirds, rng=rng)
            end
        end
    end

    # Rewire the remaining interactions to another plant of the same layer
    for (layer, plant, connector) in without_receiver
        rewire_to_other_plant!(layers=layers, layer=layer, plant=plant, bird=connector,
                               plants_in_layer=plants_in_layer[layer],
                               receiver_layer=receiver_layer,
                               hummingbirds=hummingbirds, rng=rng)
    end

    check_class_I_network(layers=layers, empirical_layers=empirical_layers,
                          decoupled=decoupled, retained_layers=retained_layers,
                          hummingbirds=hummingbirds,
                          n_rewired_to_other_plants=length(without_receiver))

    return (layers=layers, n_rewired_to_other_plants=length(without_receiver))
end

#--------------------------------------------------------------#
#-Class II: decoupling at the level of interactions------------#
#--------------------------------------------------------------#

# Whether rewiring the interaction of a bird in a layer would separate the two
# layers: the bird is the last dual-role bird and this is its last interaction
# in that layer. Null model II-B never makes this move.
# Arguments: layers holds the current matrices; layer and bird locate the
# interaction. Returns true when the move would separate the layers.
function separates_layers(; layers::Vector{Matrix{Int64}}, layer::Int64, bird::Int64)
    return length(dual_role_birds(layers=layers)) == 1 && sum(layers[layer][:, bird]) == 1
end

# Build one null network of class II. Interactions of dual-role birds are drawn
# one at a time, with uniform probability and regardless of which bird holds
# them, and rewired to a randomly chosen receiving bird at the same plant. A
# connector that loses all of its interactions in one layer becomes a
# single-role bird and can receive from then on. Rewiring continues until
# dual-role birds retain at most (1 - share) of the interactions they held in
# the empirical network. When no interaction can be rewired at its own plant,
# one is rewired to another plant of the same layer, provided its plant keeps
# at least one partner in that layer.
#
# Two rules cover the cases in which nothing can be rewired. Null model II-B
# never separates the two layers, so it excludes the moves that would, and
# stops when no other move is left (Antigua, whose single connector is the only
# frugivorous bird, and St. Eustatius, where every frugivorous bird also visits
# flowers). Null model II then finishes the separation with the rule of class
# I: the remaining dual-role birds each retain one layer (St. Eustatius).
#
# Arguments: empirical_layers holds the matrices to decouple; hummingbirds
# marks the hummingbirds; share is the share of the dual-role birds'
# interactions to rewire (1.0 in null model II, 0.5 in null model II-B); rng
# is the generator of the replicate. Returns the null layers, the number of
# interactions rewired to another plant and the rule that stopped the rewiring.
function build_class_II_network(; empirical_layers::Vector{Matrix{Int64}},
                                  hummingbirds::Vector{Bool}, share::Float64,
                                  rng::Random.Xoshiro)
    layers = [copy(empirical_layers[POLLINATION]), copy(empirical_layers[DISPERSAL])]
    plants_in_layer = [vec(sum(empirical_layers[layer]; dims=2)) .> 0
                       for layer in (POLLINATION, DISPERSAL)]
    target = (1.0 - share) * count_dual_role_interactions(layers=empirical_layers)
    keep_coupled = share < 1.0   # null model II-B never separates the layers
    n_rewired_to_other_plants = 0
    stop_rule = "target reached"

    while count_dual_role_interactions(layers=layers) > target
        receiver_layer = bird_layer(layers=layers)
        candidates = dual_role_interactions(layers=layers)
        if keep_coupled
            filter!(x -> !separates_layers(layers=layers, layer=x[1], bird=x[3]), candidates)
        end

        # Draw one interaction with a receiving bird at its plant and rewire it
        rewirable = filter(x -> !isempty(receiving_birds(layers=layers, layer=x[1], plant=x[2],
                                                         receiver_layer=receiver_layer,
                                                         hummingbirds=hummingbirds)),
                           candidates)
        if !isempty(rewirable)
            layer, plant, bird = rand(rng, rewirable)
            rewire_at_same_plant!(layers=layers, layer=layer, plant=plant, bird=bird,
                                  receiver_layer=receiver_layer, hummingbirds=hummingbirds,
                                  rng=rng)
            continue
        end

        # No interaction can be rewired at its own plant: rewire one to another
        # plant of the same layer, keeping every plant in its layer
        rewirable = filter(x -> sum(layers[x[1]][x[2], :]) > 1 &&
                                !isempty(receiving_cells_at_other_plants(
                                    layers=layers, layer=x[1], plant=x[2],
                                    plants_in_layer=plants_in_layer[x[1]],
                                    receiver_layer=receiver_layer, hummingbirds=hummingbirds)),
                           candidates)
        if !isempty(rewirable)
            layer, plant, bird = rand(rng, rewirable)
            rewire_to_other_plant!(layers=layers, layer=layer, plant=plant, bird=bird,
                                   plants_in_layer=plants_in_layer[layer],
                                   receiver_layer=receiver_layer, hummingbirds=hummingbirds,
                                   rng=rng)
            n_rewired_to_other_plants += 1
            continue
        end

        # Nothing can be rewired. Null model II-B stops here.
        if keep_coupled
            stop_rule = "stopped: no interaction could be rewired without separating the layers"
            break
        end

        # Null model II finishes the separation with the rule of class I
        connectors = dual_role_birds(layers=layers)
        allowed = allowed_assignments(layers=layers, decoupled=connectors,
                                      hummingbirds=hummingbirds)
        isempty(allowed) && error("The rule of class I cannot finish the separation.")
        finished = build_class_I_network(empirical_layers=layers, decoupled=connectors,
                                         retained_layers=rand(rng, allowed),
                                         hummingbirds=hummingbirds, rng=rng)
        layers = finished.layers
        n_rewired_to_other_plants += finished.n_rewired_to_other_plants
        stop_rule = "finished with the rule of null model I"
        break
    end

    check_class_II_network(layers=layers, empirical_layers=empirical_layers,
                           hummingbirds=hummingbirds, share=share,
                           stopped_early=startswith(stop_rule, "stopped"),
                           n_rewired_to_other_plants=n_rewired_to_other_plants)

    return (layers=layers, n_rewired_to_other_plants=n_rewired_to_other_plants,
            stop_rule=stop_rule)
end

#--------------------------------------------------------------#
#-Building one null network of either class--------------------#
#--------------------------------------------------------------#

# Build one null network of the requested class from an empirical network.
# Class I decouples all connectors (share 1.0) or draws the connectors to
# decouple (share 0.5), draws one allowed assignment of retained layers at
# random, and rewires. Class II rewires interactions one at a time.
# Arguments: empirical_layers holds the matrices; connectors lists the
# empirical connectors; class is "species" (class I) or "interactions"
# (class II); share is the share of connectors or of dual-role interactions
# decoupled; hummingbirds marks the hummingbirds; cache stores the allowed
# assignments already computed for this island; rng is the generator of the
# replicate. Returns a NamedTuple with the same fields for both classes: the
# null layers, the decoupled connectors and their retained layers (class I),
# the number of allowed assignments (class I), the number of interactions
# rewired to another plant, and the rule that stopped the rewiring (class II).
function build_null_network(; empirical_layers::Vector{Matrix{Int64}},
                              connectors::Vector{Int64}, class::String, share::Float64,
                              hummingbirds::Vector{Bool},
                              cache::Dict{Vector{Int64}, Vector{Vector{Int64}}},
                              rng::Random.Xoshiro)
    if class == "species"
        if share == 1.0
            decoupled = connectors
            allowed = get!(cache, decoupled) do
                allowed_assignments(layers=empirical_layers, decoupled=decoupled,
                                    hummingbirds=hummingbirds)
            end
            isempty(allowed) && error("No assignment of retained layers separates the layers.")
        else
            decoupled, allowed = draw_decoupled_connectors(layers=empirical_layers,
                                                           connectors=connectors, share=share,
                                                           hummingbirds=hummingbirds,
                                                           cache=cache, rng=rng)
        end
        retained_layers = rand(rng, allowed)
        null = build_class_I_network(empirical_layers=empirical_layers, decoupled=decoupled,
                                     retained_layers=retained_layers,
                                     hummingbirds=hummingbirds, rng=rng)
        return (layers=null.layers, decoupled=decoupled, retained_layers=retained_layers,
                n_allowed_assignments=length(allowed),
                n_rewired_to_other_plants=null.n_rewired_to_other_plants,
                stop_rule=missing)
    else
        null = build_class_II_network(empirical_layers=empirical_layers,
                                      hummingbirds=hummingbirds, share=share, rng=rng)
        return (layers=null.layers, decoupled=Int64[], retained_layers=Int64[],
                n_allowed_assignments=missing,
                n_rewired_to_other_plants=null.n_rewired_to_other_plants,
                stop_rule=null.stop_rule)
    end
end

#--------------------------------------------------------------#
#-Building and writing the ensemble of one null model----------#
#--------------------------------------------------------------#

# Describe the layer retained by each decoupled connector, for example
# "Bananaquit=pollination; Caribbean Elaenia=dispersal".
# Arguments: birds holds the bird labels; decoupled and retained_layers come
# from build_null_network. Returns a String (empty for class II).
function describe_retained_layers(; birds::Vector{String}, decoupled::Vector{Int64},
                                    retained_layers::Vector{Int64})
    return join(["$(birds[decoupled[i]])=$(LAYER_NAMES[retained_layers[i]])"
                 for i in eachindex(decoupled)], "; ")
end

# Describe the layers in which each empirical connector interacts in a null
# network, for example "Bananaquit=pollination; Caribbean Elaenia=both".
# Arguments: birds holds the bird labels; connectors lists the empirical
# connectors; layers holds the null matrices. Returns a String.
function describe_connector_states(; birds::Vector{String}, connectors::Vector{Int64},
                                     layers::Vector{Matrix{Int64}})
    layer = bird_layer(layers=layers)
    state(connector) = layer[connector] == 0 ? "both" : LAYER_NAMES[layer[connector]]
    return join(["$(birds[connector])=$(state(connector))" for connector in connectors], "; ")
end

# Build the ensemble of one null model for every island and write it to disk.
# Arguments: source_dir is src/, with the empirical matrices and the
# hummingbird list; model names the null model in the file names ("I",
# "I-A", "II" or "II-B"); class is "species" or "interactions"; share is the
# share of connectors (class I) or of dual-role interactions (class II)
# decoupled; seed is the base seed of the null model; n_null is the number of
# null networks per island; output_dir receives the matrices and must be new
# or empty. Writes output_dir/multi_pol/B_<island>_null_<model>_<id>.csv,
# output_dir/multi_sd/S_..., and output_dir/network_manifest.csv (one row per
# null network). Returns the manifest as a DataFrame.
function build_null_ensemble(; source_dir::String, model::String, class::String,
                               share::Float64, seed::Int64, n_null::Int64,
                               output_dir::String)
    class in ("species", "interactions") || error("class must be \"species\" or \"interactions\".")
    0.0 < share <= 1.0 || error("share must lie in (0, 1].")
    ispath(output_dir) && !isempty(readdir(output_dir)) &&
        error("Output folder is not empty: $output_dir. Move or delete it first.")
    mkpath(joinpath(output_dir, "multi_pol"))
    mkpath(joinpath(output_dir, "multi_sd"))

    networks = read_empirical_networks(source_dir=source_dir)
    hummingbird_path = joinpath(source_dir, "hummingbirds.csv")
    hummingbird_names = read_hummingbirds(path=hummingbird_path)
    manifest = NamedTuple[]

    for net in networks
        empirical_layers = [net.B, net.S]
        hummingbirds = hummingbird_mask(birds=net.birds, hummingbirds=hummingbird_names)
        any(net.S[:, hummingbirds] .> 0) &&
            error("Empirical hummingbird-fruit interaction in $(net.network).")
        connectors = dual_role_birds(layers=empirical_layers)
        cache = Dict{Vector{Int64}, Vector{Vector{Int64}}}()

        # Every null network starts from the empirical network with its own seed
        for null_id in 1:n_null
            network_seed = null_network_seed(base_seed=seed, network=net.network, null_id=null_id)
            null = build_null_network(empirical_layers=empirical_layers, connectors=connectors,
                                      class=class, share=share, hummingbirds=hummingbirds,
                                      cache=cache, rng=Xoshiro(network_seed))

            # Write the two layers with the empirical labels and order
            name = "$(net.network)_null_$(model)_$(lpad(string(null_id), 4, '0'))"
            pol_file = "multi_pol/B_$(name).csv"
            sd_file = "multi_sd/S_$(name).csv"
            write_null_matrix(matrix=null.layers[POLLINATION], plants=net.plants,
                              birds=net.birds, path=joinpath(output_dir, pol_file))
            write_null_matrix(matrix=null.layers[DISPERSAL], plants=net.plants,
                              birds=net.birds, path=joinpath(output_dir, sd_file))

            push!(manifest, (
                network = net.network,
                model = model,
                null_id = null_id,
                seed = network_seed,
                base_seed = seed,
                n_connectors = length(connectors),
                n_decoupled = count(bird_layer(layers=null.layers)[connectors] .> 0),
                n_allowed_assignments = null.n_allowed_assignments,
                retained_layers = describe_retained_layers(birds=net.birds, decoupled=null.decoupled,
                                                           retained_layers=null.retained_layers),
                connector_states = describe_connector_states(birds=net.birds, connectors=connectors,
                                                             layers=null.layers),
                dual_role_interactions = count_dual_role_interactions(layers=null.layers),
                dual_role_interactions_empirical = count_dual_role_interactions(layers=empirical_layers),
                n_rewired_to_other_plants = null.n_rewired_to_other_plants,
                stop_rule = null.stop_rule,
                pol_file = pol_file,
                sd_file = sd_file,
                pol_sha256 = net.pol_sha256,
                sd_sha256 = net.sd_sha256,
                hummingbird_sha256 = bytes2hex(sha256(read(hummingbird_path))),
                julia_version = string(VERSION)
            ))
        end
    end

    manifest_df = DataFrame(manifest)
    CSV.write(joinpath(output_dir, "network_manifest.csv"), manifest_df)
    return manifest_df
end
