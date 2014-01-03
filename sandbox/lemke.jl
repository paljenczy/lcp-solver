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
    pivot!(A, b, basic, pivot, enter)
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
        pivot!(A, b, basic, pivot, enter)
    end

    println("Maximum number of iterations reached!")
end



function cplm(i::Integer, n::Integer)
    if i <= n
        i+n
    else
        i-n
    end
end



function pivot!(A::AbstractArray{Float64,2}, b::Array{Float64,1}, basic::Array{Integer,1}, pivot::Integer, enter::Integer)
    const d = -A[pivot, enter]
    const v = b[pivot] / d
    A[pivot, basic[pivot]] = -1.0
    const Ip, Jp, Vp = findnz(A[pivot, :])
    for j in Jp
        A[pivot, j] /= d
    end
    A[pivot, enter] = 0.0
    const Ie, Je, Ve = findnz(A[:, enter])
    for j in Jp
        for i in Ie
            A[i, j] += A[i, enter] * A[pivot, j]
        end
    end
    for i in Ie
        b[i] += A[i, enter] * v
        A[i, enter] = 0.0
    end
    b[pivot] = v
    basic[pivot] = enter
end

