"""
Penalizes values in vector that are below a threshold 
"""
function minPenalty(x::Vector{Float64}, l::Float64)
    thresh = l .- x
    return sum(thresh .+ abs.(thresh))
end

"""
Penalizes values in vector that are above a threshold
"""
function maxPenalty(x::Vector{Float64}, l::Float64)
    thresh = x .- l
    return sum(thresh .+ abs.(thresh))
end

#####
#Composite loss functions: this is ugly but differentiable :)
#####
function lossFunc(xyz::Matrix{Float64}, lengths::Vector{Float64}, forces::Vector{Float64}, receiver::Receiver)
    loss = 0.0
    for (id, w) in zip(receiver.OBJids, receiver.OBJweights)
        if id == -1
            loss += 0.0
        elseif id == 0 #TARGET OBJ
            loss += w * norm(xyz - receiver.XYZtarget)
        elseif id == 1 #LENGTH VARIATION OBJ
            loss += w * -reduce(-, extrema(lengths))
        elseif id == 2 #FORCE VARIATION OBJ
            loss += w * -reduce(-, extrema(forces))
        elseif id == 3 #∑FL
            loss += w * dot(lengths, forces)
        elseif id == 4 #minimum length
            loss += w * minPenalty(lengths, receiver.MinLength)
        elseif id == 5 #maximum length
            loss += w * maxPenalty(lengths, receiver.MaxLength)
        elseif id == 6 # minimum force
            loss += w * minPenalty(forces, receiver.MinForce)
        elseif id == 7 # maximum force
            loss += w * maxPenalty(forces, receiver.MaxForce)
        end
    end
    
    return loss
end
