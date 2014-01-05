# Test code for the function "lemke(M, q, kwargs)" found in "lemke.jl"
# It assumes that all test problems have optimal solutions

function lemketest(;precision :: Integer = 8)
    const testdir :: String = ""
    const codedir :: String = "../"
    const testproblems :: Array{String,1} = ["lemketest1.prb", "lemketest2.prb"]
    const numtests :: Integer = length(testproblems)
    
    include(codedir * "problemio.jl")
    include(codedir * "lemke.jl")
    
    for i in 1:numtests
        print("Test " * string(i) * " / " * string(numtests) * " ... ")
        M, q, z0, sol = open(importlcp, testdir * testproblems[i])
        z, err = lemke(M, q)
        if round(z, precision) == sol && err == 0
            println("passed.")
        else
            println("FAILED!")
        end
    end
end
