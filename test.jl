include("blocksys.jl")
using .blocksys

function main() 
  
  t = SparseMatrix(UInt64(12), UInt64(3))
  #=
  setCell!(t, UInt64(6), UInt64(6), 3.0)
  println(getCell(t, UInt64(6), UInt64(6)))
  setCell!(t, UInt64(6), UInt64(9), 1.0)
  println(getCell(t, UInt64(6), UInt64(9)))
  printMatrix(t)
  =#
  t = readMatrix("testMatrix.txt")
  b = readVector("testVector.txt")
  printMatrix(t)
  println(b)
  gaussElimination!(t, b)
  println("---------------------")
  printMatrix(t)
  println("---------------------")
  println(b) #expected resutl ~[-0.11664806432945962, 0.45355931057479104, -0.10887220790385349, 0.13142512919740418, -0.026441062380215718, 0.054427259042966754]
  t = readMatrix("testMatrix.txt")
  b = readVector("testVector.txt")
  w = gaussEliminationMajor!(t, b)
  printMatrix(t)
  println(w)
end

main()