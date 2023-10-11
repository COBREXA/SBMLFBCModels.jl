
import SparseArrays: spzeros, sparse

A.reactions(model::SBMLModel)::Vector{String} = model.reaction_ids
A.metabolites(model::SBMLModel)::Vector{String} = model.metabolite_ids

A.n_reactions(model::SBMLModel)::Int = length(model.reaction_ids)
A.n_metabolites(model::SBMLModel)::Int = length(model.metabolite_ids)

function A.stoichiometry(model::SBMLModel)::A.SparseMat

    # find the vector size for preallocation
    nnz = 0
    for (_, r) in model.sbml.reactions
        for _ in r.reactants
            nnz += 1
        end
        for _ in r.products
            nnz += 1
        end
    end

    Rows = Int[]
    Cols = Int[]
    Vals = Float64[]
    sizehint!(Rows, nnz)
    sizehint!(Cols, nnz)
    sizehint!(Vals, nnz)

    row_idx = Dict(k => i for (i, k) in enumerate(model.metabolite_ids))
    for (ridx, rid) in enumerate(model.reaction_ids)
        r = model.sbml.reactions[rid]
        for sr in r.reactants
            push!(Rows, model.metabolite_idx[sr.species])
            push!(Cols, ridx)
            push!(Vals, isnothing(sr.stoichiometry) ? -1.0 : -sr.stoichiometry)
        end
        for sr in r.products
            push!(Rows, model.metabolite_idx[sr.species])
            push!(Cols, ridx)
            push!(Vals, isnothing(sr.stoichiometry) ? 1.0 : sr.stoichiometry)
        end
    end
    return sparse(Rows, Cols, Vals, n_metabolites(model), n_reactions(model))
end

"""
$(TYPEDSIGNATURES)

Get the lower and upper flux bounds of model [`SBMLModel`](@ref). Throws `DomainError` in
case if the SBML contains mismatching units.
"""
function A.bounds(model::SBMLModel)::Tuple{Vector{Float64},Vector{Float64}}
    # There are multiple ways in SBML to specify a lower/upper bound. There are
    # the "global" model bounds that we completely ignore now because no one
    # uses them. In reaction, you can specify the bounds using "LOWER_BOUND"
    # and "UPPER_BOUND" parameters, but also there may be a FBC plugged-in
    # parameter name that refers to the parameters.  We extract these, using
    # the units from the parameters. For unbounded reactions we use -Inf or Inf
    # as a default.

    common_unit = ""

    function get_bound(rid, fld, param, default)
        rxn = model.sbml.reactions[rid]
        param_name = SBML.mayfirst(getfield(rxn, fld), param)
        param = get(
            rxn.kinetic_parameters,
            param_name,
            get(model.sbml.parameters, param_name, default),
        )
        unit = SBML.mayfirst(param.units, "")
        if unit != ""
            if common_unit != ""
                if unit != common_unit
                    throw(
                        DomainError(
                            unit,
                            "The SBML file uses multiple units; loading would need conversion",
                        ),
                    )
                end
            else
                common_unit = unit
            end
        end
        return param.value
    end

    return (
        get_bound.(
            model.reaction_ids,
            :lower_bound,
            "LOWER_BOUND",
            Ref(SBML.Parameter(value = -Inf)),
        ),
        get_bound.(
            model.reaction_ids,
            :upper_bound,
            "UPPER_BOUND",
            Ref(SBML.Parameter(value = Inf)),
        ),
    )
end

"""
$(TYPEDSIGNATURES)

Balance vector of a [`SBMLModel`](@ref). This is always zero.
"""
A.balance(model::SBMLModel)::A.SparseVec = spzeros(Float64, n_metabolites(model))

"""
$(TYPEDSIGNATURES)

Objective of the [`SBMLModel`](@ref). Notably, SBML may have multiple stored
objectives; this function primarily takes the one specified by the
`active_objective`. If no objectives are specified, the old-style objective
specification via `OBJECTIVE_COEFFICIENT` is tried too.
"""
function A.objective(model::SBMLModel)::A.SparseVec
    res = spzeros(n_reactions(model))

    objective = get(model.sbml.objectives, model.active_objective, nothing)
    if isnothing(objective) && length(model.sbml.objectives) == 1
        objective = first(values(model.sbml.objectives))
    end
    if !isnothing(objective)
        direction = objective.type == "maximize" ? 1.0 : -1.0
        for (rid, coef) in objective.flux_objectives
            res[model.reaction_idx[rid]] = float(direction * coef)
        end
    else
        # old-style objectives
        for (rid, r) in model.sbml.reactions
            oc = get(r.kinetic_parameters, "OBJECTIVE_COEFFICIENT", nothing)
            isnothing(oc) || (res[model.reaction_idx[rid]] = float(oc.value))
        end
    end
    return res
end

A.genes(model::SBMLModel)::Vector{String} = model.gene_ids
A.n_genes(model::SBMLModel)::Int = length(model.gene_ids)
A.reaction_gene_association_dnf(
    model::SBMLModel,
    rid::String,
)::Maybe{A.GeneAssociationDNF} =
    maybemap(parse_grr, model.sbml.reactions[rid].gene_product_association)

function A.reaction_gene_products_available(
    model::SBMLModel,
    rid::String,
    available::Function,
)
    ev(x::SBML.GPAAnd) = all(ev, x.terms)
    ev(x::SBML.GPAOr) = any(ev, x.terms)
    ev(x::SBML.GPARef) = available(x.gene_product)

    maybemap(ev, model.reaction[rid].gene_product_association)
end

A.metabolite_formula(model::SBMLModel, mid::String)::Maybe{A.MetaboliteFormula} =
    maybemap(parse_formula, model.sbml.species[mid].formula)
A.metabolite_compartment(model::SBMLModel, mid::String) =
    model.sbml.species[mid].compartment
A.metabolite_charge(model::SBMLModel, mid::String)::Maybe{Int} =
    model.sbml.species[mid].charge

function A.reaction_stoichiometry(m::SBMLModel, rid::String)::Dict{String,Float64}
    s = Dict{String,Float64}()
    default1(x) = isnothing(x) ? 1 : x
    for sr in m.sbml.reactions[rid].reactants
        s[sr.species] = get(s, sr.species, 0.0) - default1(sr.stoichiometry)
    end
    for sr in m.sbml.reactions[rid].products
        s[sr.species] = get(s, sr.species, 0.0) + default1(sr.stoichiometry)
    end
    return s
end

A.reaction_name(model::SBMLModel, rid::String) = model.sbml.reactions[rid].name
A.metabolite_name(model::SBMLModel, mid::String) = model.sbml.species[mid].name
A.gene_name(model::SBMLModel, gid::String) = model.sbml.gene_products[gid].name
A.reaction_annotations(model::SBMLModel, rid::String) =
    sbml_import_cvterms(model.sbml.reactions[rid].sbo, model.sbml.reactions[rid].cv_terms)
A.metabolite_annotations(model::SBMLModel, mid::String) =
    sbml_import_cvterms(model.sbml.species[mid].sbo, model.sbml.species[mid].cv_terms)
A.gene_annotations(model::SBMLModel, gid::String) = sbml_import_cvterms(
    model.sbml.gene_products[gid].sbo,
    model.sbml.gene_products[gid].cv_terms,
)
A.reaction_notes(model::SBMLModel, rid::String) =
    sbml_import_notes(model.sbml.reactions[rid].notes)
A.metabolite_notes(model::SBMLModel, mid::String) =
    sbml_import_notes(model.sbml.species[mid].notes)
A.gene_notes(model::SBMLModel, gid::String) =
    sbml_import_notes(model.sbml.gene_products[gid].notes)
