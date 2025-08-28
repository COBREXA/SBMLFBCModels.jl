
"""
$(TYPEDSIGNATURES)

A helper for handling `nothing`s.
"""
maybemap(fn, x) = isnothing(x) ? nothing : fn(x)

"""
$(TYPEDSIGNATURES)

A helper for producing predictable unique sequences. Might be faster if
compacting would be done directly in sort().
"""
function sortunique(x)
    o = collect(x)
    sort!(o)
    put = prevind(o, firstindex(o))
    for i in eachindex(o)
        if put >= firstindex(o) && o[i] == o[put]
            # we already have this one
            continue
        else
            put = nextind(o, put)
            if put != i
                o[put] = o[i]
            end
        end
    end
    o[begin:put]
end

"""
$(TYPEDSIGNATURES)

Parse `SBML.GeneProductAssociation` structure and convert it to a strictly
positive DNF. Negation (`SBML.GPANot`) is not supported.
"""
function parse_grr(gpa::SBML.GeneProductAssociation)::A.GeneAssociationDNF

    function fold_and(dnfs::Vector{Vector{Vector{String}}})::Vector{Vector{String}}
        if isempty(dnfs)
            [String[]]
        else
            sortunique(
                sortunique(String[l; r]) for l in dnfs[1] for r in fold_and(dnfs[2:end])
            )
        end
    end

    dnf(x::SBML.GPARef) = [[x.gene_product]]
    dnf(x::SBML.GPAOr) = sortunique(vcat(dnf.(x.terms)...))
    dnf(x::SBML.GPAAnd) = fold_and(dnf.(x.terms))

    return dnf(gpa)
end

"""
$(TYPEDSIGNATURES)

Convert a gene association DNF to the corresponding `SBML.jl` structure.
"""
function unparse_grr(x::A.GeneAssociationDNF)::SBML.GeneProductAssociation
    or(x) = length(x) == 1 ? first(x) : SBML.GPAOr(x)
    and(x) = length(x) == 1 ? first(x) : SBML.GPAAnd(x)
    or([and([SBML.GPARef(j) for j in i]) for i in x])
end

"""
$(TYPEDSIGNATURES)

Parse a formula in format `C2H6O` into a `MetaboliteFormula`, which is
basically a dictionary of atom counts in the molecule.
"""
function parse_formula(f::String)::A.MetaboliteFormula
    res = Dict{String,Int}()
    pattern = @r_str "([A-Z][a-z]*)([1-9][0-9]*)?"

    for m in eachmatch(pattern, f)
        res[m.captures[1]] = isnothing(m.captures[2]) ? 1 : parse(Int, m.captures[2])
    end

    return res
end

"""
$(TYPEDSIGNATURES)

Format `MetaboliteFormula` to `String`.
"""
function unparse_formula(f::A.MetaboliteFormula)::String
    return join(["$elem$n" for (elem, n) in f])
end

function parse_sbml_identifiers_org_uri(uri::String)::Tuple{String,String}
    m = match(r"^http://identifiers.org/([^/]+)/(.*)$", uri)
    isnothing(m) ? ("RESOURCE_URI", uri) : (m[1], m[2])
end

function sbml_import_cvterms(sbo::Maybe{String}, cvs::Vector{SBML.CVTerm})::A.Annotations
    res = A.Annotations()
    isnothing(sbo) || (res["sbo"] = [sbo])
    for cv in cvs
        cv.biological_qualifier == :is || continue
        for (id, val) in parse_sbml_identifiers_org_uri.(cv.resource_uris)
            push!(get!(res, id, []), val)
        end
    end
    return res
end

function sbml_export_cvterms(annotations::A.Annotations)::Vector{SBML.CVTerm}
    isempty(annotations) && return []
    length(annotations) == 1 && haskey(annotations, "sbo") && return []
    [
        SBML.CVTerm(
            biological_qualifier = :is,
            resource_uris = [
                id == "RESOURCE_URI" ? val : "http://identifiers.org/$id/$val" for
                (id, vals) in annotations if id != "sbo" for val in vals
            ],
        ),
    ]
end

function sbml_export_sbo(annotations::A.Annotations)::Maybe{String}
    haskey(annotations, "sbo") || return nothing
    if length(annotations["sbo"]) != 1
        throw(
            DomainError(
                annotations["sbo"],
                "data loss: SBO term is not unique for SBML export",
            ),
        )
    end
    return annotations["sbo"][1]
end

sbml_import_notes(notes::Maybe{String})::A.Notes =
    isnothing(notes) ? A.Notes() : A.Notes("" => [notes])

sbml_export_notes(notes::A.Notes)::Maybe{String} =
    isempty(notes) ? nothing :
    collect(keys(notes)) == [""] && length(notes[""]) == 1 ? notes[""][1] :
    throw(DomainError(notes, "data loss: structured notes can not exported to SBML"))
