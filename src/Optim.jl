function FDMsolve!(;host = "127.0.0.1", port = 2000)
    #start server
    println("SERVER OPENED--FIRST RUN MIGHT TAKE A WHILE :) ")

    ## PERSISTENT LOOP
    server = WebSockets.listen(host, port) do ws
        # FOR EACH MESSAGE SENT FROM CLIENT
        for msg in ws
            println("MSG RECEIVED")
            problem = JSON.parsefile(msg)

            # CLOSE INSTRUCTION
            if problem == "close"
                HTTP.WebSockets.send(ws, "CONNECTION ENDED")
                close(server)
                println("SERVER ENDED BY CLIENT")
                return
            end

            # MAIN ALGORITHM
            if problem["Valid"] == true #check all required info is provided

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

                # objective function
                function obj(q::Vector{Float64}, p)
                    xyznew = solve_explicit(q, Cn, Cf, Pn, xyzf)
                    xyzfull = fullXYZ(xyznew, xyzf, N, F)
                    lengths = norm.(eachrow(C * xyzfull))
                    forces = q .* lengths

                    #loss
                    loss = 0.0
                    for (id, w) in zip(objids, objweights)
                        if id == 0 #TARGET OBJ
                            loss += w * norm(xyznew - xyz_target)
                        elseif id == 1 #LENGTH VARIATION OBJ
                            loss += w * -reduce(-, extrema(lengths))
                        elseif id == 2 #FORCE VARIATION OBJ
                            loss += w * -reduce(-, extrema(forces))
                        elseif id == 3 #∑FL
                            loss += dot(lengths, forces)
                        end
                    end

                    return loss, xyzfull
                end

                #optimization parameters
                lb = problem["LowerBound"]
                ub = problem["UpperBound"]
                abstol = problem["AbsTolerance"]
                freq = problem["UpdateFrequency"]
                maxiter = problem["MaxIterations"]
                show = problem["ShowIterations"]

                #trace
                if show
                    i = 0
                    iters = Vector{Vector{Float64}}()
                    losses = Vector{Float64}()
                end

                #callback function
                function cb(q::Vector{Float64}, loss::Float64, xyz::Matrix{Float64})
                    if show && i % freq == 0
                        push!(iters, depcopy(q))
                        push!(losses, loss)

                        #send intermediate message
                        msgout = json(Dict("Finished" => false,
                            "Iter" => i, 
                            "Loss" => loss,
                            "Q" => q, 
                            "X" => xyz[:,1], 
                            "Y" => xyz[:,2], 
                            "Z" => xyz[:,3],
                            "Qtrace" => iters,
                            "Losstrace" => losses))
                            
                        HTTP.WebSockets.send(ws, msgout)
                        i += 1
                        return false
                    else
                        i += 1
                        return false
                    end
                end

                # OPTIMIZATION
                opf = Optimization.OptimizationFunction(obj, Optimization.AutoZygote())
                opp = Optimization.OptimizationProblem(obf, q,
                    p = SciMLBase.NullParameters(),
                    lb = lb,
                    ub = ub)

                sol = Optimization.solve(opp, NLopt.LD_LBFGS(),
                    abstol = abstol,
                    maxiters = maxiter,
                    callback = cb)

                # PARSING SOLUTION
                xyz_final = solve_explicit(sol.u, Cn, Cf, Pn, xyzf)
                xyz_full_final = fullXYZ(xyz_final, xyzf, N, F)

                msgout = Dict("Finished" => true,
                    "Iter" => i,
                    "Loss" => sol.minimum,
                    "Q" => sol.u,
                    "X" => xyz_full_final[:, 1],
                    "Y" => xyz_full_final[:, 2],
                    "Z" => xyz_full_final[:, 3],
                    "Qtrace" => iters,
                    "Losstrace" => losses)

                HTTP.WebSockets.send(ws, json(msgout))
            end

        end
    end

    println("SERVER CLOSED")
    close(server)
    return
end
