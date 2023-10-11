
"""
$(TYPEDSIGNATURES)

Load and return a SBML XML model in `file_name`.
"""
load(::Type{SBMLModel}, file_name::String)::SBMLModel = SBMLModel(SBML.readSBML(file_name))

"""
$(TYPEDSIGNATURES)

Write a given SBML model to `file_name`.
"""
save(model::SBMLModel, file_name::String) = SBML.writeSBML(m.sbml, file_name)
