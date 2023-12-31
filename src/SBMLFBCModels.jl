"""
$(README)
"""
module SBMLFBCModels

using DocStringExtensions

import AbstractFBCModels as A
import SBML
import SparseArrays

include("types.jl")
include("interface.jl")
include("convert.jl")
include("io.jl")
include("misc.jl")

end # module SBMLFBCModels
