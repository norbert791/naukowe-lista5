include("blocksys.jl")
using Test
using .blocksys

function accessTests()
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
end

function gaussEliminationTests()
  t = readMatrix("testMatrix.txt")
  b = readVector("testVector.txt")
  gaussElimination!(t, b)
  #computed with https://calculator-online.net/gaussian-elimination-calculator/
  expected::Vector{Float64} = [-0.11664806432945962, 0.45355931057479104, -0.10887220790385349, 0.13142512919740418, -0.026441062380215718, 0.054427259042966754]
  success = true
  for i in 1:(length(b))
    if abs(expected[i] - b[i]) > 0.001
      success = false
      break
    end
  end

  @assert success "Gauss Elimination failed"
end

function main() 
  accessTests()
  gaussEliminationTests()

  println("All tests passed")
end

main()