#author: Norbert Jaśniewicz
include("blocksys.jl")
using .blocksys
using SparseArrays

"""
  read matrix A=ARGS[1] then compute b = A * x.
  Save it to ARGS[2] where x = (1,1,...,1) transposed
"""
function main()
  m = readMatrix(ARGS[1])
  result = m.vals * ones(Float64, m.size)
  writeVector(ARGS[2], result)
end

main()