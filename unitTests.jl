include("blocksys.jl")
using Test
using .blocksys

function main() 
  t = SparseMatrix(UInt64(12), UInt64(3))
  setCell!(t, UInt64(6), UInt64(6), 3.0)
  @assert getCell(t, UInt64(6), UInt64(6)) == 3.0 "getCell != 3.0"
  setCell!(t, UInt64(6), UInt64(9), 1.0)
  @assert getCell(t, UInt64(6), UInt64(9)) == 1.0 "getCell != 1.0"

  @test_throws BoundsError getCell(t, UInt64(0), UInt64(9))
  @test_throws BoundsError getCell(t, UInt64(13), UInt64(9))
  @test_throws BoundsError getCell(t, UInt64(6), UInt64(0))
  @test_throws BoundsError getCell(t, UInt64(6), UInt64(13))

  @test_throws BoundsError getCell(t, UInt64(0), UInt64(9))
  @test_throws BoundsError getCell(t, UInt64(13), UInt64(9))
  @test_throws BoundsError getCell(t, UInt64(6), UInt64(0))
  @test_throws BoundsError getCell(t, UInt64(6), UInt64(13))

  println("All tests passed")
end

main()