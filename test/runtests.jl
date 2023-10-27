
import AbstractFBCModels as A
import SBMLFBCModels as S
using Test

@testset "SBMLFBCModels tests" begin
    A.run_fbcmodel_type_tests(S.SBMLModel)
end
