#author: Norbert Jaśniewicz
include("blocksys.jl")
import LinearAlgebra

using .blocksys
using SparseArrays

"""
  read matrix A=ARGS[1], compute b = A * x where x = (1,1,...,1)^T.
  Compute x' with given method ARGS[2]:
  1-Gauss
  2-Gauss with pivot
  3-LU
  4-LU with pivot
  Save result with relative error to file ARGS[3]
"""
function main()
  m = readMatrix(ARGS[1])
  x = ones(Float64, m.size)
  b = m.vals * x
  flag::String = ARGS[2]
  xPrim = []
  if flag == "1"
    gaussElimination!(m, b)
    xPrim = b
  elseif flag == "2"
    xPrim = gaussEliminationPivot!(m, b)
  elseif flag == "3"
    luDecomposition!(m)
    computeFromLU!(m, b)
    xPrim = b
  elseif flag == "4"
    swp = luDecompositionPivot!(m)
    xPrim = computeFromLUPivot!(m, swp, b)
  else
    error("Wrong flag. Pass x in {1,2,3,4} as 2nd argument")
  end

  relError::Float64 = LinearAlgebra.norm(x - xPrim) / LinearAlgebra.norm(x)
  
  open(ARGS[3], "w") do file
    write(file, "$(relError)\n")
    for val in xPrim
      write(file, "$(val)\n")
    end
  end  
end

main()