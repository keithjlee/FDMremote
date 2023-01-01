function FDMecho!(;host = "127.0.0.1", port = 2000)
    #start server
    println("SERVER OPENED--FIRST RUN MIGHT TAKE A WHILE :) ")

    ## PERSISTENT LOOP
    server = WebSockets.listen!(host, port) do ws
        # FOR EACH MESSAGE SENT FROM CLIENT
        for msg in ws
            println("MSG RECEIVED")

            # CLOSE INSTRUCTION
            if msg == "CLOSE"
                send(ws, "CONNECTION ENDED")
                close(server)
                println("SERVER ENDED BY CLIENT")
                return
            end

            send(ws, msg)
        end
    end
end

function FDMsolve!(;host = "127.0.0.1", port = 2000)
    #start server
    println("###############################################")
    println("SERVER OPENED")
    println("###############################################")

    ## initialize variable
    msgout = Dict

    ## min/max values
    minforce = 1.0
    maxforce = 10000.
    minlength = 1.0
    maxlength = 1000.

    counter = 0

    ## PERSISTENT LOOP
    server = WebSockets.listen!(host, port) do ws
        # FOR EACH MESSAGE SENT FROM CLIENT
        for msg in ws
            println("MSG RECEIVED")

            # CLOSE INSTRUCTION
            if msg == "CLOSE"
                WebSockets.send(ws, "CONNECTION ENDED")
                close(server)
                println("SERVER ENDED BY CLIENT")
                return
            end

            if msg == "init" || msg == "Hello World"
                println("CONNECTION INITIALIZED")
                continue
            end

            try
                problem = JSON.parse(msg)
                # MAIN ALGORITHM
                if haskey(problem, "Valid") && problem["Valid"]
                    problem = JSON.parse(msg)
                    println("READING DATA")
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

                    if length(objids) == 0
                        objids = [-1]
                        objweights = [1.0]
                    end

                    # values for min/max force/length
                    minforce = Float64(problem["MinForce"])
                    maxforce = Float64(problem["MaxForce"])
                    minlength = Float64(problem["MinLength"])
                    maxlength = Float64(problem["MaxLength"])

                    # solving
                    println("OPTIMIZING")
                    if counter == 0
                        println("First run will take a while! :--)")
                        counter += 1
                    end
                    # objective function
                    function obj(q::Vector{Float64}, p)
                        xyznew = solve_explicit(q, Cn, Cf, Pn, xyzf)
                        xyzfull = fullXYZ(xyznew, xyzf, N, F)
                        lengths = norm.(eachrow(C * xyzfull))
                        forces = q .* lengths

                        #####
                        #Composite loss functions: this is ugly but differentiable :)
                        #####
                        loss = 0.0
                        for (id, w) in zip(objids, objweights)
                            if id == -1
                                loss += 0.0
                            elseif id == 0 #TARGET OBJ
                                loss += w * norm(xyznew - xyz_target)
                            elseif id == 1 #LENGTH VARIATION OBJ
                                loss += w * -reduce(-, extrema(lengths))
                            elseif id == 2 #FORCE VARIATION OBJ
                                loss += w * -reduce(-, extrema(forces))
                            elseif id == 3 #∑FL
                                loss += w * dot(lengths, forces)
                            elseif id == 4 #minimum length
                                loss += w * minPenalty(lengths, minlength)
                            elseif id == 5 #maximum length
                                loss += w * maxPenalty(lengths, maxlength)
                            elseif id == 6 # minimum force
                                loss += w * minPenalty(forces, minforce)
                            elseif id == 7 # maximum force
                                loss += w * maxPenalty(forces, maxforce)
                            end
                        end

                        return loss, xyzfull
                    end

                    #optimization parameters
                    lb = repeat([problem["LowerBound"]], ne)
                    ub = repeat([problem["UpperBound"]], ne)
                    abstol = problem["AbsTolerance"]
                    freq = problem["UpdateFrequency"]
                    maxiter = problem["MaxIterations"]
                    show = problem["ShowIterations"]

                    #trace
                    i = 0
                    iters = Vector{Vector{Float64}}()
                    losses = Vector{Float64}()

                    #callback function
                    function cb(q::Vector{Float64}, loss::Float64, xyz::Matrix{Float64})
                        if show && i % freq == 0
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
                    opp = Optimization.OptimizationProblem(opf, q,
                        p = SciMLBase.NullParameters(),
                        lb = lb,
                        ub = ub)

                    sol = Optimization.solve(opp, NLopt.LD_LBFGS(),
                        abstol = abstol,
                        maxiters = maxiter,
                        callback = cb)

                    println("SOLUTION FOUND")
                    # PARSING SOLUTION
                    xyz_final = solve_explicit(sol.u, Cn, Cf, Pn, xyzf)
                    xyz_full_final = fullXYZ(xyz_final, xyzf, N, F)

                    # cb(sol.u, sol.minimum, xyz_full_final)

                    msgout = Dict("Finished" => true,
                        "Iter" => i,
                        "Loss" => sol.minimum,
                        "Q" => sol.u,
                        "X" => xyz_full_final[:, 1],
                        "Y" => xyz_full_final[:, 2],
                        "Z" => xyz_full_final[:, 3],
                        "Losstrace" => losses)

                    WebSockets.send(ws, json(msgout))
                else
                    println("INVALID INPUT")
                end
            catch
                println("INVALID INPUT")
            end

            println("DONE")
        end
    end
end