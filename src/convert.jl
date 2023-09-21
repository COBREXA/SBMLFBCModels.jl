
"""
$(TYPEDSIGNATURES)

Convert any metabolic model to [`SBMLModel`](@ref).
"""
function Base.convert(::Type{SBMLModel}, mm::MetabolicModel)
    if typeof(mm) == SBMLModel
        return mm
    end

    mets = metabolites(mm)
    rxns = reactions(mm)
    stoi = stoichiometry(mm)
    (lbs, ubs) = bounds(mm)
    comps = _default.("compartment", metabolite_compartment.(Ref(mm), mets))
    compss = Set(comps)

    metid(x) = startswith(x, "M_") ? x : "M_$x"
    rxnid(x) = startswith(x, "R_") ? x : "R_$x"
    gprid(x) = startswith(x, "G_") ? x : "G_$x"

    return SBMLModel(
        SBML.Model(
            compartments = Dict(
                comp => SBML.Compartment(constant = true) for comp in compss
            ),
            species = Dict(
                metid(mid) => SBML.Species(
                    name = metabolite_name(mm, mid),
                    compartment = _default("compartment", comps[mi]),
                    formula = _maybemap(_unparse_formula, metabolite_formula(mm, mid)),
                    charge = metabolite_charge(mm, mid),
                    constant = false,
                    boundary_condition = false,
                    only_substance_units = false,
                    sbo = _sbml_export_sbo(metabolite_annotations(mm, mid)),
                    notes = _sbml_export_notes(metabolite_notes(mm, mid)),
                    metaid = metid(mid),
                    cv_terms = _sbml_export_cvterms(metabolite_annotations(mm, mid)),
                ) for (mi, mid) in enumerate(mets)
            ),
            reactions = Dict(
                rxnid(rid) => SBML.Reaction(
                    name = reaction_name(mm, rid),
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
                    kinetic_parameters = Dict(
                        "LOWER_BOUND" => SBML.Parameter(value = lbs[ri]),
                        "UPPER_BOUND" => SBML.Parameter(value = ubs[ri]),
                    ),
                    lower_bound = "LOWER_BOUND",
                    upper_bound = "UPPER_BOUND",
                    gene_product_association = _maybemap(
                        x -> _unparse_grr(SBML.GeneProductAssociation, x),
                        reaction_gene_association(mm, rid),
                    ),
                    reversible = true,
                    sbo = _sbml_export_sbo(reaction_annotations(mm, rid)),
                    notes = _sbml_export_notes(reaction_notes(mm, rid)),
                    metaid = rxnid(rid),
                    cv_terms = _sbml_export_cvterms(reaction_annotations(mm, rid)),
                ) for (ri, rid) in enumerate(rxns)
            ),
            gene_products = Dict(
                gprid(gid) => SBML.GeneProduct(
                    label = gid,
                    name = gene_name(mm, gid),
                    sbo = _sbml_export_sbo(gene_annotations(mm, gid)),
                    notes = _sbml_export_notes(gene_notes(mm, gid)),
                    metaid = gprid(gid),
                    cv_terms = _sbml_export_cvterms(gene_annotations(mm, gid)),
                ) for gid in genes(mm)
            ),
            active_objective = "objective",
            objectives = Dict(
                "objective" => SBML.Objective(
                    "maximize",
                    Dict(rid => oc for (rid, oc) in zip(rxns, objective(mm)) if oc != 0),
                ),
            ),
        ),
    )
end
