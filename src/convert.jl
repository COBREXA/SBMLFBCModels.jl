
"""
$(TYPEDSIGNATURES)

Convert any metabolic model to [`SBMLFBCModel`](@ref).
"""
function Base.convert(::Type{SBMLFBCModel}, mm::A.AbstractFBCModel)
    if typeof(mm) == SBMLFBCModel
        return mm
    end

    function default(d::T, x::Maybe{T})::T where {T}
        isnothing(x) ? d : x
    end

    mets = A.metabolites(mm)
    rxns = A.reactions(mm)
    stoi = A.stoichiometry(mm)
    (lbs, ubs) = A.bounds(mm)
    comps = default.("compartment", A.metabolite_compartment.(Ref(mm), mets))
    compss = Set(comps)

    #TODO re-think this
    metid(x) = startswith(x, "M_") ? x : "M_$x"
    rxnid(x) = startswith(x, "R_") ? x : "R_$x"
    gprid(x) = startswith(x, "G_") ? x : "G_$x"

    bound_params = Dict{Float64,String}()
    for b in [lbs; ubs]
        haskey(bound_params, b) && continue
        bound_params[b] = "bound_value_$(length(bound_params)+1)"
    end

    return SBMLFBCModel(
        SBML.Model(
            compartments = Dict(
                comp => SBML.Compartment(constant = true) for comp in compss
            ),
            species = Dict(
                metid(mid) => SBML.Species(
                    name = A.metabolite_name(mm, mid),
                    compartment = default("compartment", comps[mi]),
                    formula = maybemap(unparse_formula, A.metabolite_formula(mm, mid)),
                    charge = A.metabolite_charge(mm, mid),
                    constant = false,
                    boundary_condition = false,
                    only_substance_units = false,
                    sbo = sbml_export_sbo(A.metabolite_annotations(mm, mid)),
                    notes = sbml_export_notes(A.metabolite_notes(mm, mid)),
                    metaid = metid(mid),
                    cv_terms = sbml_export_cvterms(A.metabolite_annotations(mm, mid)),
                ) for (mi, mid) in enumerate(mets)
            ),
            parameters = Dict(p => SBML.Parameter(value = v) for (v,p) in bound_params),
            reactions = Dict(
                rxnid(rid) => SBML.Reaction(
                    name = A.reaction_name(mm, rid),
                    reactants = [
                        SBML.SpeciesReference(
                            species = metid(mets[i]),
                            stoichiometry = -stoi[i, ri],
                            constant = true,
                        ) for
                        i in SparseArrays.nonzeroinds(stoi[:, ri]) if stoi[i, ri] <= 0
                    ],
                    products = [
                        SBML.SpeciesReference(
                            species = metid(mets[i]),
                            stoichiometry = stoi[i, ri],
                            constant = true,
                        ) for
                        i in SparseArrays.nonzeroinds(stoi[:, ri]) if stoi[i, ri] > 0
                    ],
                    lower_bound = bound_params[lbs[ri]],
                    upper_bound = bound_params[ubs[ri]],
                    gene_product_association = maybemap(
                        unparse_grr,
                        map(gs -> map(gprid, gs), A.reaction_gene_association_dnf(mm, rid)),
                    ),
                    reversible = true,
                    sbo = sbml_export_sbo(A.reaction_annotations(mm, rid)),
                    notes = sbml_export_notes(A.reaction_notes(mm, rid)),
                    metaid = rxnid(rid),
                    cv_terms = sbml_export_cvterms(A.reaction_annotations(mm, rid)),
                ) for (ri, rid) in enumerate(rxns)
            ),
            gene_products = Dict(
                gprid(gid) => SBML.GeneProduct(
                    label = gid,
                    name = A.gene_name(mm, gid),
                    sbo = sbml_export_sbo(A.gene_annotations(mm, gid)),
                    notes = sbml_export_notes(A.gene_notes(mm, gid)),
                    metaid = gprid(gid),
                    cv_terms = sbml_export_cvterms(A.gene_annotations(mm, gid)),
                ) for gid in A.genes(mm)
            ),
            active_objective = "objective",
            objectives = Dict(
                "objective" => SBML.Objective(
                    "maximize",
                    Dict(rid => oc for (rid, oc) in zip(rxns, A.objective(mm)) if oc != 0),
                ),
            ),
        ),
    )
end
