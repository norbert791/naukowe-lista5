include("blocksys.jl")
using .blocksys

function main() 
  t = SparseMatrix(UInt64(12), UInt64(3))
  setCell!(t, UInt64(6), UInt64(6), 3.0)
  printMatrix(t)
  println(getCell(t, UInt64(6), UInt64(6)))
  setCell!(t, UInt64(6), UInt64(9), 1.0)
  println(getCell(t, UInt64(6), UInt64(9)))
end

main()