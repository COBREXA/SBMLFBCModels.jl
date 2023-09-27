"""
$(TYPEDEF)

Thin wrapper around the model from SBML.jl library. Allows easy conversion from
SBML to any other model format.

# Fields
$(TYPEDFIELDS)
"""
struct SBMLFBCModel <: A.AbstractFBCModel
    sbml::SBML.Model
    reaction_ids::Vector{String}
    reaction_idx::Dict{String,Int}
    metabolite_ids::Vector{String}
    metabolite_idx::Dict{String,Int}
    gene_ids::Vector{String}
    active_objective::String
end

"""
$(TYPEDEF)

Construct the SBML model and add the necessary cached indexes, possibly choosing an active objective.
"""
function SBMLFBCModel(sbml::SBML.Model, active_objective::String = "")
    rxns = sort(collect(keys(sbml.reactions)))
    mets = sort(collect(keys(sbml.species)))
    genes = sort(collect(keys(sbml.gene_products)))

    SBMLFBCModel(
        sbml,
        rxns,
        Dict(rxns .=> eachindex(rxns)),
        mets,
        Dict(mets .=> eachindex(mets)),
        genes,
        active_objective,
    )
end

const Maybe = A.Maybe
