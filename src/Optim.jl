function FDMsolve!(;host = "127.0.0.1", port = 2000)
    #start server
    println("###############################################")
    println("###############SERVER OPENED###################")
    println("###############################################")

    counter = 0

    ## PERSISTENT LOOP
    server = WebSockets.listen!(host, port) do ws

        # FOR EACH MESSAGE SENT FROM CLIENT
        for msg in ws

            # ACKNOWLEDGE
            println("MSG RECEIVED")

            # CLOSE INSTRUCTION
            if msg == "CLOSE"
                WebSockets.send(ws, "CONNECTION ENDED")
                close(server)
                println("SERVER ENDED BY CLIENT")
                return
            end

            # FIRST MESSAGE
            if msg == "init" || msg == "Hello World"
                println("CONNECTION INITIALIZED")
                continue
            end

            # ANALYSIS
            try
                # DESERIALIZE MESSAGE
                problem = JSON.parse(msg)

                # MAIN ALGORITHM
                # IF PROBLEM VALID
                if haskey(problem, "Valid") && problem["Valid"]

                    println("READING DATA")

                    # CONVERT MESSAGE
                    receiver = Receiver(problem)

                    # SOLVE
                    println("OPTIMIZING")
                    if counter == 0
                        println("First run will take a while! :--)")
                        counter += 1
                    end
                    
                    # OPTIMIZATION
                    FDMoptim!(receiver, ws)
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