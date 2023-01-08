#author: Norbert Jaśniewicz

include("blocksys.jl")
include("blockmat.jl")

import LinearAlgebra

using .blocksys
using .matrixgen

"""
Perform tests on blocksys functions with random matrix and vector. The vector is parametrised by:
  ARGS[1] = size
  ARGS[2] = l
  ARGS[3] = cond
"""
function main()
  tempFileName::String = "randomTestTempFile.txt"
  blockmat(parse(Int, ARGS[1]), parse(Int, ARGS[2]), parse(Float64, ARGS[3]), tempFileName)
  m = readMatrix(tempFileName)
  b = rand(Float64, m.size)

  function relError(x::Vector{Float64})
    temp = m.vals * x
    return LinearAlgebra.norm(b - temp) / LinearAlgebra.norm(b)
  end

  xPrim = gaussElimination!(deepcopy(m), deepcopy(b))
  println("Gauss elimination relError: $(relError(xPrim))")
  xPrim = gaussEliminationPivot!(deepcopy(m), deepcopy(b))
  println("Gauss elimination with pivot relError: $(relError(xPrim))")
  temp = deepcopy(m)
  luDecomposition!(temp)
  xPrim = computeFromLU!(temp, deepcopy(b))
  println("LU + mult relError: $(relError(xPrim))")
  temp = deepcopy(m)
  swp = luDecompositionPivot!(temp)
  xPrim = computeFromLUPivot!(temp, swp, deepcopy(b))
  println("LU with pivot + mult relError: $(relError(xPrim))")
end




main()