#--------------------------------------------------------------#
#-Building the decoupling null models for the Lesser Antilles--#
#--------------------------------------------------------------#

# Build the null networks of null models I, I-A, II and II-B. Each one
# decouples the pollination and seed-dispersal layers of an island network
# and keeps its other topological properties fixed. Class I decouples bird
# species: all connectors (the dual-role birds, null model I) or half of
# them (I-A) retain a single layer. Class II decouples interactions: those
# of dual-role birds are rewired until the layers separate (II) or until
# dual-role birds retain half of their interactions (II-B).
#
# Inputs:  src/multi_pol/B_<island>.csv
#          src/multi_sd/S_<island>.csv
#          src/hummingbirds.csv
# Outputs: src/null_<model>/multi_pol/B_<island>_null_<model>_<id>.csv
#          src/null_<model>/multi_sd/S_<island>_null_<model>_<id>.csv
#          src/null_<model>/network_manifest.csv

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

# share is the share of connectors (class I) or of dual-role interactions
# (class II) to decouple. seed is the base seed of the null model. The seed
# of each null network comes from it, the island and the replicate number.
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

# Write each null model to src/null_<model>/. The folder must be new or empty.
for null_model in null_models
    println("Building null model $(null_model.model)...")
    build_null_ensemble(source_dir=srcdir(), model=null_model.model,
                        class=null_model.class, share=null_model.share,
                        seed=null_model.seed, n_null=n_null,
                        output_dir=srcdir("null_$(null_model.model)"))
end
