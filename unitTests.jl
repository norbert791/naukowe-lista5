#author Norbert Jaśniewicz
include("blocksys.jl")
using Test
using .blocksys

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


function gaussEliminationMajorTests()
  t = readMatrix("testMatrix.txt")
  b = readVector("testVector.txt")
  b = gaussEliminationMajor!(t, b)
  #computed with https://calculator-online.net/gaussian-elimination-calculator/
  expected::Vector{Float64} = [-0.11664806432945962, 0.45355931057479104, -0.10887220790385349, 0.13142512919740418, -0.026441062380215718, 0.054427259042966754]
  success = true
  
  for i in 1:(length(b))
    if abs(expected[i] - b[i]) > 0.001
      success = false
      break
    end
  end

  @assert success "Gaus Elimination with partial major failed"
end

function luDecompositionTest()
  t = readMatrix("testMatrix.txt")
  b = readVector("testVector.txt")
  luDecomposition!(t)
  computeFromLU!(t, b)
  
  println(t.vals)
  
  expected::Vector{Float64} = [-0.11664806432945962, 0.45355931057479104, -0.10887220790385349, 0.13142512919740418, -0.026441062380215718, 0.054427259042966754]
  success = true
  
  for i in 1:(length(b))
    if abs(expected[i] - b[i]) > 0.001
      success = false
      break
    end
  end

  @assert success "LU decomposition failed"
end

function main() 
  gaussEliminationTests()
  gaussEliminationMajorTests()
  luDecompositionTest()

  println("All tests passed")
end

main()