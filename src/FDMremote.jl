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
include("Optim.jl")
include("Solve.jl")

export FDMsolve!
export FDMecho!

end # module FDMremote
