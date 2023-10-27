
import AbstractFBCModels as A
import SBMLFBCModels: SBMLFBCModel

using Test

@testset "SBMLFBCModels tests" begin
    A.run_fbcmodel_type_tests(SBMLFBCModel)

    modeldir = joinpath(@__DIR__, "test-models")
    mkpath(modeldir)

    for (name, hash) in ["e_coli_core" => "xxx", "iJO1366" => "yyy"]
        path = joinpath(modeldir, "$name.xml")
        A.download_data_file("http://bigg.ucsd.edu/static/models/$name.xml", path, hash)
        A.run_fbcmodel_file_tests(SBMLFBCModel, path; name)
    end
end
