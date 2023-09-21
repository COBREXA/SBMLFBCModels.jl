
"""
$(TYPEDSIGNATURES)

A helper for producing predictable unique sequences. Might be faster if
compacting would be done directly in sort().
"""
function _sortunique(x)
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
positive DNF [`GeneAssociation`](@ref). Negation (`SBML.GPANot`) is not
supported.
"""
function _parse_grr(gpa::SBML.GeneProductAssociation)::GeneAssociation

    function fold_and(dnfs::Vector{Vector{Vector{String}}})::Vector{Vector{String}}
        if isempty(dnfs)
            [String[]]
        else
            _sortunique(
                _sortunique(String[l; r]) for l in dnfs[1] for r in fold_and(dnfs[2:end])
            )
        end
    end

    dnf(x::SBML.GPARef) = [[x.gene_product]]
    dnf(x::SBML.GPAOr) = _sortunique(vcat(dnf.(x.terms)...))
    dnf(x::SBML.GPAAnd) = fold_and(dnf.(x.terms))
    dnf(x) = throw(
        DomainError(
            x,
            "unsupported gene product association contents of type $(typeof(x))",
        ),
    )
    return dnf(gpa)
end

"""
$(TYPEDSIGNATURES)

Convert a GeneAssociation to the corresponding `SBML.jl` structure.
"""
function _unparse_grr(
    ::Type{SBML.GeneProductAssociation},
    x::GeneAssociation,
)::SBML.GeneProductAssociation
    SBML.GPAOr([SBML.GPAAnd([SBML.GPARef(j) for j in i]) for i in x])
end


function _parse_sbml_identifiers_org_uri(uri::String)::Tuple{String,String}
    m = match(r"^http://identifiers.org/([^/]+)/(.*)$", uri)
    isnothing(m) ? ("RESOURCE_URI", uri) : (m[1], m[2])
end

function _sbml_import_cvterms(sbo::Maybe{String}, cvs::Vector{SBML.CVTerm})::Annotations
    res = Annotations()
    isnothing(sbo) || (res["sbo"] = [sbo])
    for cv in cvs
        cv.biological_qualifier == :is || continue
        for (id, val) in _parse_sbml_identifiers_org_uri.(cv.resource_uris)
            push!(get!(res, id, []), val)
        end
    end
    return res
end

function _sbml_export_cvterms(annotations::Annotations)::Vector{SBML.CVTerm}
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

function _sbml_export_sbo(annotations::Annotations)::Maybe{String}
    haskey(annotations, "sbo") || return nothing
    if length(annotations["sbo"]) != 1
        @_io_log @error "Data loss: SBO term is not unique for SBML export" annotations["sbo"]
        return
    end
    return annotations["sbo"][1]
end

function _sbml_import_notes(notes::Maybe{String})::Notes
    isnothing(notes) ? Notes() : Notes("" => [notes])
end

function _sbml_export_notes(notes::Notes)::Maybe{String}
    isempty(notes) || @_io_log @error "Data loss: notes not exported to SBML" notes
    nothing
end
