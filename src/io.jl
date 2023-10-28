
A.load(::Type{SBMLFBCModel}, file_name::String)::SBMLFBCModel =
    SBMLFBCModel(SBML.readSBML(file_name))

A.save(model::SBMLFBCModel, file_name::String) = SBML.writeSBML(model.sbml, file_name)

A.filename_extensions(::Type{SBMLFBCModel}) = ["xml"]
