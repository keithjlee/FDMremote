module FDMremote

using LinearAlgebra
using Statistics
using SparseArrays
using Optimization
using OptimizationNLopt
using Zygote
using HTTP
using JSON

include("FDM.jl")
include("types.jl")
include("optimization.jl")
include("analysis.jl")
include("communication.jl")

export FDMsolve!

end # module FDMremote
