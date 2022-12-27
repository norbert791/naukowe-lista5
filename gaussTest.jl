#author: Norbert Jaśniewicz
include("blocksys.jl")
using .blocksys

"""
  read matrix A=ARGS[1] and vector b=ARGS[2] then compute x : b = A * x. Save x to file ARGS[3]_normal.txt and ARGS[3]_pivot.txt
"""
function main()
  m1 = readMatrix(ARGS[1])
  b1 = readVector(ARGS[2])
  m2 = readMatrix(ARGS[1])
  b2 = readVector(ARGS[2])
  gaussElimination!(m1, b1)
  b2 = gaussEliminationPivot!(m2, b2)
  writeVector("$(ARGS[3])_normal.txt", b1)
  writeVector("$(ARGS[3])_pivot.txt", b2)
end

main()