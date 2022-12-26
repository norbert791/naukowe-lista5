include("blocksys.jl")
using .blocksys

function test1() 
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

function test2()
 t = readMatrix("examples/dane500_000/A.txt")
 b = readVector("examples/dane500_000/b.txt")
 #printMatrix(t)
 gaussElimination!(t, b)
 #println(b)
 #t = readMatrix("examples/dane500_000/A.txt")
 #b = readVector("examples/dane500_000/b.txt")
 #w = gaussEliminationMajor!(t, b)
 #println(w)
end

function test3()
  t = readMatrix("testMatrix.txt")
  b = readVector("testVector.txt")
  luDecomposition!(t)
  computeFromLU!(t, b)
  @show b
  #@show t.vals
end

function main()
  #@time test2()
  test3()
end

main()