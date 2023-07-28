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
