struct Receiver

    #FORCE DENSITY
    Q::Vector{Float64}

    #NETWORK INFORMATION
    N::Vector{Int64}
    F::Vector{Int64}
    
    #OBJECTIVES
    OBJids::Vector{Int64}
    OBJweights::Vector{Float64}

    #PARAMETERS
    MinForce::Float64
    MaxForce::Float64
    MinLength::Float64
    MaxLength::Float64
    Tol::Float64
    Freq::Int64
    MaxIter::Int64
    Show::Bool

    #DERIVED VALUES
    LB::Vector{Float64}
    UB::Vector{Float64}

    XYZf::Matrix{Float64}
    XYZtarget::Matrix{Float64}
    Pn::Matrix{Float64}

    C::SparseMatrixCSC{Int64, Int64}
    Cn::SparseMatrixCSC{Int64, Int64}
    Cf::SparseMatrixCSC{Int64, Int64}

    #constructor
    function Receiver(problem::Dict)
        # point geometry
        x = Float64.(problem["X"])
        z = Float64.(problem["Z"])
        y = Float64.(problem["Y"])

        # initial force densities
        q = Float64.(problem["Q"])

        # global info
        ne = Int(problem["Ne"])
        nn = Int(problem["Nn"])

        # free/fixed
        N = Int.(problem["Njulia"])
        F = Int.(problem["Fjulia"])

        # loads
        px = Float64.(problem["Px"])
        py = Float64.(problem["Py"])
        pz = Float64.(problem["Pz"])

        # connectivity
        i = Int.(problem["Ijulia"])
        j = Int.(problem["Jjulia"])
        v = Int.(problem["V"])
        C = sparse(i, j, v, ne, nn)
        Cn = C[:, N]
        Cf = C[:, F]

        # in matrix form
        xyz = hcat(x, y, z)
        xyzf = xyz[F, :]
        xyz_target = xyz[N, :]
        Pn = hcat(px, py, pz)

        # objective functions
        objids = Int64.(problem["OBJids"])
        objweights = Float64.(problem["OBJweights"])

        # check for null objective
        if length(objids) == 0
            objids = [-1]
            objweights = [1.0]
        end

        # values for min/max force/length
        minforce = Float64(problem["MinForce"])
        maxforce = Float64(problem["MaxForce"])
        minlength = Float64(problem["MinLength"])
        maxlength = Float64(problem["MaxLength"])

        #optimization parameters
        lb = repeat([problem["LowerBound"]], ne)
        ub = repeat([problem["UpperBound"]], ne)
        abstol = problem["AbsTolerance"]
        freq = problem["UpdateFrequency"]
        maxiter = problem["MaxIterations"]
        show = problem["ShowIterations"]

        receiver = new(q, 
            N, 
            F, 
            objids, 
            objweights, 
            minforce, 
            maxforce, 
            minlength, 
            maxlength, 
            abstol, 
            freq, 
            maxiter, 
            show, 
            lb, 
            ub, 
            xyzf, 
            xyz_target, 
            Pn, 
            C, 
            Cn, 
            Cf)

        return receiver
    end

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

### optimiztaion
function FDMoptim!(receiver::Receiver, ws::WebSockets.WebSocket)
    # objective function
    if length(receiver.OBJids) == 1 && receiver.OBJids[1] == -1

        println("Single Solve")

        xyznew = solve_explicit(receiver.q, receiver.Cn, receiver.Cf, receiver.Pn, receiver.XYZf)

        xyz = fullXYZ(xyznew, receiver.XYZf, receiver.N, receiver.F)

        msgout = Dict("Finished" => true,
                "Iter" => 1, 
                "Loss" => 0.,
                "Q" => receiver.q, 
                "X" => xyz[:,1], 
                "Y" => xyz[:,2], 
                "Z" => xyz[:,3],
                "Losstrace" => [0.])
        
        for _ = 1:5
            WebSockets.send(ws, json(msgout))
        end

        println("message sent")
        
    else

        println("OPTIMIZING")
        function obj(q::Vector{Float64}, p)

            xyznew = solve_explicit(q, receiver.Cn, receiver.Cf, receiver.Pn, receiver.XYZf)

            xyzfull = fullXYZ(xyznew, receiver.XYZf, receiver.N, receiver.F)
            
            lengths = norm.(eachrow(receiver.C * xyzfull))
            forces = q .* lengths

            loss = lossFunc(xyznew, lengths, forces, receiver)

            return loss, xyzfull
        end

        #trace
        i = 0
        iters = Vector{Vector{Float64}}()
        losses = Vector{Float64}()

        #callback function
        function cb(q::Vector{Float64}, loss::Float64, xyz::Matrix{Float64})
            if receiver.Show && i % receiver.Freq == 0
                push!(iters, deepcopy(q))
                push!(losses, loss)

                #send intermediate message
                msgout = Dict("Finished" => false,
                    "Iter" => i, 
                    "Loss" => loss,
                    "Q" => q, 
                    "X" => xyz[:,1], 
                    "Y" => xyz[:,2], 
                    "Z" => xyz[:,3],
                    "Losstrace" => losses)
                    
                WebSockets.send(ws, json(msgout))
                i += 1
                println("Iteration $i")
                return false
            else
                i += 1
                println("Iteration $i")
                return false
            end
        end

        # OPTIMIZATION
        opf = Optimization.OptimizationFunction(obj, Optimization.AutoZygote())
        opp = Optimization.OptimizationProblem(opf, receiver.Q,
            p = SciMLBase.NullParameters(),
            lb = receiver.LB,
            ub = receiver.UB)

        sol = Optimization.solve(opp, NLopt.LD_LBFGS(),
            abstol = receiver.Tol,
            maxiters = receiver.MaxIter,
            callback = cb)

        println("SOLUTION FOUND")
        # PARSING SOLUTION
        xyz_final = solve_explicit(sol.u, receiver.Cn, receiver.Cf, receiver.Pn, receiver.XYZf)
        xyz_full_final = fullXYZ(xyz_final, receiver.XYZf, receiver.N, receiver.F)


        msgout = Dict("Finished" => true,
            "Iter" => i,
            "Loss" => sol.minimum,
            "Q" => sol.u,
            "X" => xyz_full_final[:, 1],
            "Y" => xyz_full_final[:, 2],
            "Z" => xyz_full_final[:, 3],
            "Losstrace" => losses)

        

        for _ = 1:5
            WebSockets.send(ws, json(msgout))
        end

        println("Final message sent")
    end
end