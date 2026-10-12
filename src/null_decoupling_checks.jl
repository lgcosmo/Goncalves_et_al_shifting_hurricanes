#--------------------------------------------------------------#
#-Checks of every null network of the decoupling null models---#
#--------------------------------------------------------------#

# Each null network of null models I to II-B (Extended Data Fig. 1) is
# checked before it is written to disk, and a failed check stops the run. The
# checks verify the properties that every null model holds fixed, and the
# properties specific to each class. They use the helpers of
# src/null_decoupling_functions.jl, which must be included before this file.
# build_class_I_network and build_class_II_network call the checks on every
# null network they build, so no unchecked network reaches build_null_ensemble.
#
# Inputs:  the null and empirical layers passed by the building functions
# Outputs: none. Each check returns nothing or stops the run with an error.
# Used by: src/null_decoupling_functions.jl, included after it by
#          scripts/01_network_preparation/build_null_decoupling.jl

# Properties shared by all null models: the matrices stay binary; species
# richness, the identity of the plants in each layer and the number of
# interactions in each layer are unchanged; no single-role bird becomes a
# connector or loses an interaction; hummingbirds receive no fruits; and plant
# degrees change only through the interactions rewired to another plant.
# Arguments: layers holds the null matrices and empirical_layers the empirical
# ones; hummingbirds marks the hummingbirds; n_rewired_to_other_plants is the
# number of interactions rewired to another plant. Returns nothing.
function check_shared_properties(; layers::Vector{Matrix{Int64}},
                                   empirical_layers::Vector{Matrix{Int64}},
                                   hummingbirds::Vector{Bool},
                                   n_rewired_to_other_plants::Int64)
    B, S = layers
    B_emp, S_emp = empirical_layers

    all(x -> x in (0, 1), B) && all(x -> x in (0, 1), S) ||
        error("A null network has a nonbinary entry.")

    # Every bird keeps at least one interaction
    birds_active = vec(sum(B; dims=1) .+ sum(S; dims=1)) .> 0
    birds_active_emp = vec(sum(B_emp; dims=1) .+ sum(S_emp; dims=1)) .> 0
    birds_active == birds_active_emp || error("A bird lost all of its interactions.")

    # Single-role birds stay in their layer and keep their interactions
    single = bird_layer(layers=empirical_layers) .> 0
    all(bird_layer(layers=layers)[single] .== bird_layer(layers=empirical_layers)[single]) ||
        error("A single-role bird changed its layer or became a connector.")
    for layer in (POLLINATION, DISPERSAL)
        all(layers[layer][:, single] .>= empirical_layers[layer][:, single]) ||
            error("A single-role bird lost an interaction.")
    end

    # Every plant stays in the layers it belonged to
    for layer in (POLLINATION, DISPERSAL)
        (vec(sum(layers[layer]; dims=2)) .> 0) == (vec(sum(empirical_layers[layer]; dims=2)) .> 0) ||
            error("A plant changed the layers it belongs to.")
    end

    sum(B) == sum(B_emp) && sum(S) == sum(S_emp) ||
        error("A layer changed its number of interactions.")

    any(S[:, hummingbirds] .> 0) && error("A hummingbird received a fruit interaction.")

    # Each interaction rewired to another plant changes two plants by one
    degree_change = sum(sum(abs.(vec(sum(layers[layer]; dims=2)) .-
                                 vec(sum(empirical_layers[layer]; dims=2))))
                        for layer in (POLLINATION, DISPERSAL))
    degree_change <= 2 * n_rewired_to_other_plants ||
        error("Plant degrees changed beyond the interactions rewired to another plant.")

    return nothing
end

# Properties of class I: each decoupled connector keeps every interaction of
# its retained layer and none of the other layer, untouched connectors keep
# exactly their empirical interactions, and only untouched connectors still
# interact in both layers. Arguments: as in build_class_I_network, plus layers
# (the null matrices) and n_rewired_to_other_plants. Returns nothing.
function check_class_I_network(; layers::Vector{Matrix{Int64}},
                                 empirical_layers::Vector{Matrix{Int64}},
                                 decoupled::Vector{Int64}, retained_layers::Vector{Int64},
                                 hummingbirds::Vector{Bool}, n_rewired_to_other_plants::Int64)
    check_shared_properties(layers=layers, empirical_layers=empirical_layers,
                            hummingbirds=hummingbirds,
                            n_rewired_to_other_plants=n_rewired_to_other_plants)

    for (i, connector) in enumerate(decoupled)
        retained = retained_layers[i]
        rewired = rewired_layer(retained_layer=retained)
        all(layers[retained][:, connector] .>= empirical_layers[retained][:, connector]) ||
            error("A decoupled connector lost an interaction in its retained layer.")
        sum(layers[rewired][:, connector]) == 0 ||
            error("A decoupled connector kept an interaction in its rewired layer.")
    end

    untouched = setdiff(dual_role_birds(layers=empirical_layers), decoupled)
    for connector in untouched, layer in (POLLINATION, DISPERSAL)
        layers[layer][:, connector] == empirical_layers[layer][:, connector] ||
            error("An untouched connector changed its interactions.")
    end
    dual_role_birds(layers=layers) == untouched ||
        error("The birds interacting in both layers are not the untouched connectors.")

    return nothing
end

# Properties of class II: dual-role birds retain at most (1 - share) of their
# empirical interactions, unless the rewiring stopped early, and null model
# II-B keeps the two layers coupled. Arguments: as in build_class_II_network,
# plus layers (the null matrices), stopped_early (true when null model II-B
# stopped before reaching its target) and n_rewired_to_other_plants.
# Returns nothing.
function check_class_II_network(; layers::Vector{Matrix{Int64}},
                                  empirical_layers::Vector{Matrix{Int64}},
                                  hummingbirds::Vector{Bool}, share::Float64,
                                  stopped_early::Bool, n_rewired_to_other_plants::Int64)
    check_shared_properties(layers=layers, empirical_layers=empirical_layers,
                            hummingbirds=hummingbirds,
                            n_rewired_to_other_plants=n_rewired_to_other_plants)

    stopped_early ||
        count_dual_role_interactions(layers=layers) <=
        (1.0 - share) * count_dual_role_interactions(layers=empirical_layers) ||
        error("Dual-role birds retained too many interactions.")

    share == 1.0 || count_dual_role_interactions(layers=layers) > 0 ||
        error("Null model II-B separated the two layers.")

    return nothing
end
