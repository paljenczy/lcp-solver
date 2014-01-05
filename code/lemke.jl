# I'm trying to code a simple lemke algorithm from the CMEM book (Ch.8.)
# This is a learning exercise, so don't expect perfect code

function lemke(M :: AbstractArray{Float64,2}, q :: Array{Float64,1}; maxiter :: Integer = 1000, zerotol :: Float64 = 1e-5, verbose :: Bool = false)
    # z, err = lemke(M, q; kwargs)
    # Error codes:
    # 0: Found a solution
    # 1: Maximum number of iterations exceeded (zero vector returned)
    # 2: Algorithm ended up on a secondary ray (zero vector returned)
    const n = length(q)
    d = sparse(1.0 * (q.<0))
    if sum(d) == 0.0
        if verbose; println("Trivial solution exists."); end
        return zeros((n, 1)), 0
    end
    A :: AbstractArray{Float64,2} = [copy(M) spzeros(n) d]
    b :: Array{Float64,1} = copy(q)
    basic :: Array{Integer,1} = [n+1:2n]
    enter :: Integer = 2n+1
    pivrow :: Integer = indmin(b)
    leave :: Integer = basic[pivrow]
    pivot!(A, b, basic, pivrow, enter)
    for iter in 1:maxiter
        if leave==2n+1
            if verbose; println("I've found a solution!"); end
            solution :: Array{Float64,1} = zeros(n)
            for i in 1:n
                if basic[i] <= n
                    solution[basic[i]] = b[i]
                end
            end
            return solution, 0
        end
        enter = cplm(leave, n)
        if sum(A[:, enter] .< 0)==0
            if verbose; println("Algorithm ended up on a secondary ray"); end
            return zeros((n, 1)), 2
        end
        pivrow = 0
        v = Inf
        for i in 1:n
            if A[i, enter] < 0 && -b[i] / A[i, enter] < v
                pivrow = i
                v = -b[i] / A[i, enter]
            end
        end
        leave = basic[pivrow]
        pivot!(A, b, basic, pivrow, enter)
    end

    if verbose; println("Maximum number of iterations reached!"); end
    return zeros((n, 1)), 1
end


function pivot!(A :: AbstractArray{Float64,2}, b :: Array{Float64,1}, basic :: Array{Integer,1}, pivrow :: Integer, enter :: Integer; zerotol :: Float64 = 1e-5)
    const n :: Int32 = length(b)
    const d :: Float64 = -A[pivrow, enter]
    const v :: Float64 = b[pivrow] / d
    A[pivrow, basic[pivrow]] = -1.0
    A[pivrow, enter] = 0.0
    A[pivrow, :] /= d
    const pivrowixs :: Array{Int32,1} = [1:n][bitunpack(slice(A, :, enter) .!= 0.0)]
    const pivcolixs :: Array{Int32,1} = [1:2n+1][bitunpack(slice(A, pivrow, :) .!= 0.0)]
    
    Ap :: Array{Float64,2} = full(A[pivrowixs, pivcolixs] + A[pivrowixs, enter] * A[pivrow, pivcolixs])
    setindex!(Ap, 0.0, -zerotol .< Ap .< zerotol)
    setindex!(A, sparse(Ap), pivrowixs, pivcolixs)

    setindex!(b, b[pivrowixs] + A[pivrowixs, enter] * v, pivrowixs)

    setindex!(A, 0.0, pivrowixs, [enter])

    b[pivrow] = v

    basic[pivrow] = enter
end


function cplm(i :: Integer, n :: Integer)
    if i <= n
        i+n
    else
        i-n
    end
end

