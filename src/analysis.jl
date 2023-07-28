
### optimiztaion
function FDMoptim!(receiver::Receiver, ws::WebSockets.WebSocket)
    # objective function
    if length(receiver.OBJids) == 1 && receiver.OBJids[1] == -1

        println("SOLVING")

        xyznew = solve_explicit(receiver.Q, receiver.Cn, receiver.Cf, receiver.Pn, receiver.XYZf)

        xyz = fullXYZ(xyznew, receiver.XYZf, receiver.N, receiver.F)

        msgout = Dict("Finished" => true,
                "Iter" => 1, 
                "Loss" => 0.,
                "Q" => receiver.Q, 
                "X" => xyz[:,1], 
                "Y" => xyz[:,2], 
                "Z" => xyz[:,3],
                "Losstrace" => [0.])

        WebSockets.send(ws, json(msgout))
        
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

        for _ = 1:3
            WebSockets.send(ws, json(msgout))
        end

    end
end