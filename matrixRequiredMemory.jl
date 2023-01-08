#author: Norbert Jaśniewicz
include("blockmat.jl")
include("blocksys.jl")

using .matrixgen
using .blocksys

function generateStats(sizeRange)
  blockSize = 5
  testRuns = 10
  condCoeff::Float64 = 10.0
  filename::String = "plotsTemp.txt"

  resultAlloc = []

  for n in sizeRange
    temp = 0.0
    println(n)
    for _ in 1:testRuns
      blockmat(n, blockSize, condCoeff, filename)
      temp += @allocated(mtx = readMatrix(filename))
    end
    temp /= testRuns
    push!(resultAlloc, temp)
  end

  return resultAlloc
end

function main()
  rng = 10_000:10_000:100_000
  res1 = generateStats(rng)
  
  open("experimentsLoad.csv", "w") do file
    write(file, "n,allocated\n")
    index = 1
    for n in rng
      write(file, "$(n),$(res1[index])\n")
      index += 1
    end
  end
end

main()
