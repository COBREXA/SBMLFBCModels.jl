
import AbstractFBCModels as A
import SBMLFBCModels as S
import SBMLFBCModels: SBMLFBCModel

using Test

@testset "SBMLFBCModels tests" begin
    A.run_fbcmodel_type_tests(SBMLFBCModel)

    modeldir = joinpath(@__DIR__, "test-models")
    mkpath(modeldir)

    for (name, url, hash, ts) in [
        (
            "e_coli_core",
            "http://bigg.ucsd.edu/static/models/e_coli_core.xml",
            "b4db506aeed0e434c1f5f1fdd35feda0dfe5d82badcfda0e9d1342335ab31116",
            true,
        ),
        (
            "iJO1366",
            "http://bigg.ucsd.edu/static/models/iJO1366.xml",
            "d6d9ec61ef6f155db5bb2f49549119dc13b96f6098b403ef82ea4240b27232eb",
            true,
        ),
        (
            "ecoli_core_model",
            "http://systemsbiology.ucsd.edu/sites/systemsbiology.ucsd.edu/files/Attachments/Images/downloads/Ecoli_core/ecoli_core_model.xml",
            "78692f8509fb36534f4f9b6ade23b23552044f3ecd8b48d84d484636922ae907",
            false,
        ),
    ]
        path = joinpath(modeldir, "$name.xml")
        A.download_data_file(url, path, hash)
        A.run_fbcmodel_file_tests(SBMLFBCModel, path; name, test_save = ts)
    end

    include("misc.jl")
end
