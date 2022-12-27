#author: Norbert Jaśniewicz
include("blocksys.jl")
using .blocksys

"""
  read matrix A=ARGS[1] and compute luDecomposition then read bVector file names and output file names
  from stdio and save x: x = A * b
"""
function main()
  m1 = readMatrix(ARGS[1])
  luDecomposition!(m1)
  m2 = readMatrix(ARGS[1])
  swp = luDecompositionPivot!(m2)
  while true
    println("input: <bVector filename> <outputName>")
    inp = readline()
    if inp == "quit" || inp == "q"
      break
    end
    temp = split(inp)
    input_file::String = temp[1]
    output_file::String = temp[2]
    bVector = readVector(input_file)
    x = copy(bVector)
    computeFromLU!(m1, x)
    bVector = computeFromLUPivot!(m2, swp, bVector)
    writeVector("$(output_file)_normal.txt", x)
    writeVector("$(output_file)_pivot.txt", bVector)
  end
end

main()