# I'm trying to code a simple lemke algorithm from the CMEM book (Ch.8.)
# This is a learning exercise, so don't expect nice code :)


function lemke(M::Array{Float64,2}, q::Array{Float64,1}; maxiter::Int64=1000)
    const n = length(q)
    d = 1.0 * (q.<0)
    if sum(d) == 0.0
        println("Trivial solution exists.")
        return zeros((n, 1))
    end
    A::Array{Float64,2} = [copy(M) spzeros(n) d]
    b::Array{Float64,1} = copy(q)
    basic::Array{Int64,1} = [n+1:2n]
    enter::Int64 = 2n+1
    pivot::Int64 = indmin(b)
    leave::Int64 = basic[pivot]
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



function cplm(i::Int64, n::Int64)
    if i <= n
        i+n
    else
        i-n
    end
end



function pivot!(A::Array{Float64,2}, b::Array{Float64,1}, basic::Array{Int64,1}, pivot::Int64, enter::Int64)
    const d = -A[pivot, enter]
    const v = b[pivot] / d
    A[pivot, basic[pivot]] = -1.0
    A[pivot, :] /= d
    A[pivot, enter] = 0.0
    for i in 1:length(basic)
        if A[i, enter] != 0.0
            A[i, :] += A[i, enter] * A[pivot, :]
            b[i] += A[i, enter] * v
            A[i, enter] = 0.0
        end
    end
    b[pivot] = v
    basic[pivot] = enter
end

