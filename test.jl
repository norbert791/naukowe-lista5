include("blocksys.jl")
using .blocksys

function main() 
  t = SparseMatrix(UInt64(12), UInt64(4))
  printMatrix(t)
end

main()