function FDMsolve!(;host = "127.0.0.1", port = 2000)
    #start server
    println("###############################################")
    println("SERVER OPENED--FIRST RUN MIGHT TAKE A WHILE :) ")
    println("###############################################")

    ## initialize variable
    msgout = Dict
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
                Pn = hcat(px, py, pz)

                # solving
                xyz_new = solve_explicit(q, Cn, Cf, Pn, xyzf)
                xyzfull = fullXYZ(xyznew, xyzf, N, F)
                
            end
        end
    end
end