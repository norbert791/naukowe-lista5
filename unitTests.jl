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


function gaussEliminationPivotTests()
  t = readMatrix("testMatrix.txt")
  b = readVector("testVector.txt")
  b = gaussEliminationPivot!(t, b)
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

function luDecompositionPivotTest()
  t = readMatrix("testMatrix.txt")
  b = readVector("testVector.txt")
  swp = luDecompositionPivot!(t)
  result = computeFromLUPivot!(t, swp, b)
    
  expected::Vector{Float64} = [-0.11664806432945962, 0.45355931057479104, -0.10887220790385349, 0.13142512919740418, -0.026441062380215718, 0.054427259042966754]
  success = true
  
  for i in 1:(length(b))
    if abs(expected[i] - result[i]) > 0.001
      success = false
      break
    end
  end

  @assert success "LU decomposition with pivot failed"
end

function main() 
  gaussEliminationTests()
  gaussEliminationPivotTests()
  luDecompositionTest()
  luDecompositionPivotTest() 

  println("All tests passed")
end

main()