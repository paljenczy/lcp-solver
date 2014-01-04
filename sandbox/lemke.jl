# I'm trying to code a simple lemke algorithm from the CMEM book (Ch.8.)
# This is a learning exercise, so don't expect nice code :)


function lemke(M::AbstractArray{Float64,2}, q::Array{Float64,1}; maxiter::Integer=1000)
    const n = length(q)
    d = sparse(1.0 * (q.<0))
    if sum(d) == 0.0
        println("Trivial solution exists.")
        return zeros((n, 1))
    end
    A::AbstractArray{Float64,2} = [copy(M) spzeros(n) d]
    b::Array{Float64,1} = copy(q)
    basic::Array{Integer,1} = [n+1:2n]
    enter::Integer = 2n+1
    pivot::Integer = indmin(b)
    leave::Integer = basic[pivot]
    @time pivot!(A, b, basic, pivot, enter)
    for iter in 1:maxiter
        if leave==2n+1
            println("I've found a solution!")
            solution::Array{Float64,1} = zeros(n)
            for i in 1:n
                if basic[i] <= n
                    solution[basic[i]] = b[i]
                end
            end
            return solution
        end
        enter = cplm(leave, n)
        if sum(A[:, enter] .< 0)==0
            println("Algorithm ended up on a secondary ray")
            return
        end
        pivot = 0
        v = Inf
        for i in 1:n
            if A[i, enter] < 0 && -b[i] / A[i, enter] < v
                pivot = i
                v = -b[i] / A[i, enter]
            end
        end
        leave = basic[pivot]
        @time pivot!(A, b, basic, pivot, enter)
    end

    println("Maximum number of iterations reached!")
end



function cplm(i :: Integer, n :: Integer)
    if i <= n
        i+n
    else
        i-n
    end
end



function pivot!(A :: AbstractArray{Float64,2}, b :: Array{Float64,1}, basic :: Array{Integer,1}, pivot :: Integer, enter :: Integer)
    const n :: Int32 = length(b)
    const d :: Float64 = -A[pivot, enter]
    const v :: Float64 = b[pivot] / d
    A[pivot, basic[pivot]] = -1.0
    A[pivot, enter] = 0.0
    I, J, V = findnz(A)
    const allpivcols :: Array{Int32,1} = [1:2n+1][bitunpack(slice(A, pivot, :) .!= 0.0)]
    const allpivrows :: Array{Int32,1} = [1:n][bitunpack(slice(A, :, enter) .!= 0.0)]
    const pivotrow::Array{Bool,1} = I.==pivot
    const entercol::Array{Bool,1} = J.==enter
    const numpivrows :: Integer = length(allpivrows)
    const numpivcols :: Integer = length(allpivcols)
    V[pivotrow] /= d    
    I2 :: Array{Int32,1} = zeros(numpivrows * numpivcols)
    J2 :: Array{Int32,1} = zeros(numpivrows * numpivcols)
    V2 :: Array{Float64,1} = zeros(numpivrows * numpivcols)
    for j in 1:numpivcols, i in 1:numpivrows
        ix = (j - 1) * numpivrows + i
        I2[ix] = allpivrows[i]
        J2[ix] = allpivcols[j]
        V2[ix] = A[i, enter] * A[pivot, j]
    end
    for i in allpivrows
        b[i] += A[i, enter] * v
    end
    V[entercol] = 0.0
    b[pivot] = v
    basic[pivot] = enter
    A = sparse([I; I2], [J; J2], [V; V2])
end

