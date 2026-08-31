
A.load(::Type{SBMLFBCModel}, file_name::String)::SBMLFBCModel =
    SBMLFBCModel(SBML.readSBML(file_name))

A.load_from_string(::Type{SBMLFBCModel}, data::String)::SBMLFBCModel =
    SBMLFBCModel(SBML.readSBMLFromString(data))

A.save(model::SBMLFBCModel, file_name::String) = SBML.writeSBML(model.sbml, file_name)

A.save_to_string(model::SBMLFBCModel) = SBML.writeSBML(model.sbml)

A.filename_extensions(::Type{SBMLFBCModel}) = ["xml", "sbml"]
