# Functions to import and export LCP problems

function importlcp(s :: IOStream)
    # Line 1: number of variables
    # Line 2: number of nonzero matrix elements
    # Next numnz lines: coordinates and values of nonzero matrix elements (M)
    # Next n lines: constant vector (q)
    # Next n lines: starting vector (z0) - not yet implemented
    # Next n lines: solution vector (sol) - not yet implemented
    const n :: Integer = integer(strip(readline(s)))
    const numnz :: Integer = integer(strip(readline(s)))
    I :: Array{Int32,1} = zeros(Int32, numnz)
    J :: Array{Int32,1} = zeros(Int32, numnz)
    V :: Array{Float64,1} = zeros(Float64, numnz)
    q :: Array{Float64,1} = zeros(Float64, n)
    z0 :: Array{Float64,1} = zeros(Float64, n)
    sol :: Array{Float64,1} = zeros(Float64, n)
    for i in 1:numnz
        cell = split(strip(readline(s)))
        I[i] = integer(cell[1])
        J[i] = integer(cell[2])
        V[i] = float(cell[3])
    end
    for i in 1:n
        q[i] = float(strip(readline(s)))
    end
    for i in 1:n
        z0[i] = float(strip(readline(s)))
    end
    for i in 1:n
        sol[i] = float(strip(readline(s)))
    end
    return sparse(I, J, V, n, n), q, z0, sol
end

     
