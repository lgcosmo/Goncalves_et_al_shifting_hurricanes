#--------------------------------------------------------------#
#-Building the decoupling null models for the Lesser Antilles--#
#--------------------------------------------------------------#

# Build the null networks of null models I, I-A, II and II-B, which Extended
# Data Fig. 1 compares with the empirical networks. Each null network
# decouples the pollination and seed-dispersal layers of one island network
# while holding its remaining topological properties fixed, and so tests
# whether the coupling of the two layers by dual-role birds (the connectors)
# shapes the extinctions that hurricanes cause. Class I decouples bird
# species. All connectors (null model I) or half of them (null model I-A)
# retain a single layer. Class II decouples interactions. The cross-layer
# interactions of dual-role birds are rewired until the layers separate (null
# model II) or until dual-role birds retain half of their interactions (null
# model II-B). Every island receives 100 null networks per null model.
#
# Inputs:  src/multi_pol/B_<island>.csv and src/multi_sd/S_<island>.csv
#          (build_networks_multilayer.R)
#          src/hummingbirds.csv (hummingbird list)
# Outputs: src/null_<model>/multi_pol/B_<island>_null_<model>_<id>.csv
#          src/null_<model>/multi_sd/S_<island>_null_<model>_<id>.csv
#          src/null_<model>/network_manifest.csv (one row per null network)
# Helpers: src/null_io_functions.jl reads and writes the matrices,
#          src/null_decoupling_functions.jl rewires, and
#          src/null_decoupling_checks.jl checks every null network
# Used by: main_simulations_null_decoupling.jl simulates the null networks
#          under the empirical hurricane regime of their island

#--------------------------------------------------------------#
#-Loading packages and functions-------------------------------#
#--------------------------------------------------------------#

using DrWatson
@quickactivate "Goncalves_et_al_shifting_hurricanes"

include(srcdir("null_io_functions.jl"))
include(srcdir("null_decoupling_functions.jl"))
include(srcdir("null_decoupling_checks.jl"))

#--------------------------------------------------------------#
#-Defining the null models-------------------------------------#
#--------------------------------------------------------------#

# The four null models of the Methods. share is the share of connectors
# (class I) or of dual-role interactions (class II) that is decoupled. Each
# null network derives its seed from the base seed of its null model, its
# island and its replicate number.
n_null = 100   # null networks per island and null model

null_models = [
    (model="I",    class="species",      share=1.0, seed=20261003),
    (model="I-A",  class="species",      share=0.5, seed=20261005),
    (model="II",   class="interactions", share=1.0, seed=20261006),
    (model="II-B", class="interactions", share=0.5, seed=20261007)
]

#--------------------------------------------------------------#
#-Building the ensembles---------------------------------------#
#--------------------------------------------------------------#

# Each null model is written to src/null_<model>/, which must be new or empty
for null_model in null_models
    println("Building null model $(null_model.model)...")
    build_null_ensemble(source_dir=srcdir(), model=null_model.model,
                        class=null_model.class, share=null_model.share,
                        seed=null_model.seed, n_null=n_null,
                        output_dir=srcdir("null_$(null_model.model)"))
end
