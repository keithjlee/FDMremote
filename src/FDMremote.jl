module FDMremote

using LinearAlgebra
using Statistics
using SparseArrays
using Optimization
using OptimizationNLopt
using Zygote
using HTTP.WebSockets

include("FDM.jl")
include("Optim.jl")

export FDMsolve!

end # module FDMremote
