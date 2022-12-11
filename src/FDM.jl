#######
# contains functions for analyzing an FDM network
#######

```
Explicit solver function
```
function solve_explicit(q::Union{Vector{Float64}, Vector{Int64}}, #vector of force densities
        Cn::SparseMatrixCSC{Int64, Int64}, #index matrix of free nodes
        Cf::SparseMatrixCSC{Int64, Int64}, #index matrix of fixed nodes
        Pn::Union{Matrix{Float64}, Matrix{Int64}}, #Matrix of free node loads
        Nf::Union{Matrix{Float64}, Matrix{Int64}}) #Matrix of fixed node positions

    Q = sparse(diagm(q)) # build diagonal force density Matrix
    return (Cn' * Q * Cn) \ (Pn - Cn' * Q * Cf * Nf)
end

"""
Returns full xyz matrix including fixed points
"""
function fullXYZ(xyzcurrent::Matrix{Float64}, xyzf::Matrix{Float64}, N::Vector{Int64}, F::Vector{Int64})
    xyzunsorted = [xyzcurrent; xyzf]
    i = sortperm([N; F])

    return xyzunsorted[i, :]
end