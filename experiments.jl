#author: Norbert Jaśniewicz
include("blockmat.jl")
include("blocksys.jl")

using .matrixgen
using .blocksys

function luTest(mtx::SparseMatrix, bVector::Vector{Float64})::Vector{Float64}
  luDecomposition!(mtx)
  return computeFromLU!(mtx, bVector)
end

function luPivotTest(mtx::SparseMatrix, bVector::Vector{Float64})::Vector{Float64}
  swp = luDecompositionPivot!(mtx)
  return computeFromLUPivot!(mtx, swp, bVector)
end

gaussTest(mtx::SparseMatrix, bVector::Vector{Float64})::Vector{Float64} =
  gaussElimination!(mtx, bVector)

gaussPivotTest(mtx::SparseMatrix, bVector::Vector{Float64})::Vector{Float64} =
  gaussEliminationPivot!(mtx, bVector)

function generateStats(computationMethods::Dict{String, Function}, sizeRange)
  blockSize = 5
  testRuns = 10
  condCoeff::Float64 = 10.0
  filename::String = "plotsTemp.txt"

  resultTime = Dict([(fun, []) for fun in values(computationMethods)])
  resultAlloc = Dict([(fun, []) for fun in values(computationMethods)])

  for n in sizeRange
    blockmat(n, blockSize, condCoeff, filename)
    mtx = readMatrix(filename)
    bVector = rand(Float64, mtx.size)
    
    for fun in values(computationMethods)
      sumTime = 0.0
      allocAmount = 0
      
      for i in 1:testRuns
        println("n: $(n) fun: $(fun), run: $(i)")
        mtxCopy = deepcopy(mtx)
        bCopy = deepcopy(bVector)
        sumTime += @elapsed(fun(mtxCopy, bCopy))
        mtxCopy = deepcopy(mtx)
        bCopy = deepcopy(bVector)
        allocAmount += @allocated(fun(mtxCopy, bCopy))
      end
      sumTime /= testRuns
      allocAmount /= testRuns
      push!(resultTime[fun], sumTime) 
      push!(resultAlloc[fun], allocAmount)
    end

    
  end

  return resultTime, resultAlloc
end

function main()
  funs = Dict("gauss" => gaussTest, "gauss_pivot" => gaussPivotTest,
              "lu" => luTest, "lu_pivot" => luPivotTest)
  rng = 10_000:10_000:100_000
  res1, res2 = generateStats(funs, rng)
  
  funNames = keys(funs)

  open("experimentsTime.csv", "w") do file
    write(file, "n,")
    write(file, join(funNames, ","))
    write(file, "\n")
    index = 1
    for n in rng
      write(file, "$(n),")
      temp = [res1[funs[f]][index] for f in funNames]
      write(file, join(temp, ","))
      write(file, "\n")
      index += 1
    end
  end

  open("experimentsMemory.csv", "w") do file
    write(file, "n,")
    write(file, join(funNames, ","))
    write(file, "\n")
    index = 1
    for n in rng
      write(file, "$(n),")
      temp = [res2[funs[f]][index] for f in funNames]
      write(file, join(temp, ","))
      write(file, "\n")
      index += 1
    end
  end
end

main()