#author: Norbert Jaśniewicz
module blocksys
  import SparseArrays

  export SparseMatrix, printMatrix, gaussElimination!, gaussEliminationPivot!,
         readMatrix, readVector, luDecomposition!, computeFromLU!, luDecompositionPivot!,
         computeFromLUPivot!, writeVector

  mutable struct SparseMatrix
    size::UInt64
    subMatrixLength::UInt64
    vals::SparseArrays.SparseMatrixCSC{Float64, UInt64}
  end

  function printMatrix(mtx::SparseMatrix)
    @show mtx.vals
  end

  function gaussEliminationPriv!(mtx::SparseMatrix, bVector::Vector{Float64})
    @boundscheck if mtx.size != length(bVector)
      throw(DomainError("bVector and mtx have different length"))
    end

    #iterate over diagonal
    for rowIndex::UInt64 in 1:(mtx.size - 1)
      elem = mtx.vals[rowIndex,rowIndex]
      @boundscheck if iszero(elem)
        throw(DomainError("diagonal has 0 elem"))
      end

      blockEnd = min(rowIndex + mtx.subMatrixLength, mtx.size)
      #for every row within current block + 1 row
      for lowerRowIndex in (rowIndex + 1):(blockEnd)
        @inbounds temp = mtx.vals[lowerRowIndex,rowIndex]
        @inbounds multiplier = temp / elem
        @inbounds mtx.vals[lowerRowIndex, rowIndex] = 0.0

        lastIndex = min(rowIndex + mtx.subMatrixLength, mtx.size)
        
        #for every element in that row
        for index in (rowIndex + 1):lastIndex
          @inbounds currVal = mtx.vals[lowerRowIndex, index]
          @inbounds upperVal = mtx.vals[rowIndex, index]
          @inbounds mtx.vals[lowerRowIndex, index] = currVal - multiplier * upperVal
        end
        @inbounds bVector[lowerRowIndex] -= bVector[rowIndex] * multiplier
      end
    end
  end

  function gaussElimination!(mtx::SparseMatrix, bVector::Vector{Float64})::Vector{Float64}
    gaussEliminationPriv!(mtx::SparseMatrix, bVector::Vector{Float64})

    for i::UInt64 in length(bVector):-1:1
      for j::UInt64 in (i+1):(min(mtx.size, i + 2 * mtx.subMatrixLength - 1))
        @inbounds bVector[i] -= bVector[j] * mtx.vals[i, j]
      end
      @inbounds bVector[i] /= mtx.vals[i, i]
    end

    return bVector
  end

  function gaussEliminationPivotPriv!(mtx::SparseMatrix, bVector::Vector{Float64})::Vector{UInt64}
    @boundscheck if mtx.size != length(bVector)
      throw(DomainError("bVector and mtx have different length"))
    end
    rowPermutation::Vector{UInt64} = map(identity, 1:mtx.size)
    #iterate over diagonal
    for rowIndex in 1:(mtx.size - 1)
      
      #find major element row
      blockEnd = min(rowIndex + mtx.subMatrixLength, mtx.size)
      majorRowIndex::UInt64 = rowIndex
      
      for i in (rowIndex + 1):blockEnd
        firstElem = mtx.vals[rowPermutation[i], rowIndex]
        if abs(firstElem) > abs(mtx.vals[rowPermutation[majorRowIndex], rowIndex])
          majorRowIndex = i
        end
      end

      #update permutation vector
      rowPermutation[rowIndex], rowPermutation[majorRowIndex] = rowPermutation[majorRowIndex], rowPermutation[rowIndex]
      #Now majorRowIndex is the row number in actual representation
      majorRowIndex = rowPermutation[rowIndex]

      elem = mtx.vals[majorRowIndex, rowIndex]

      @boundscheck if iszero(elem)
        throw(DomainError("0 sub-column detected"))
      end

      #for every row within current block + 1 row
      for nonMajorRowIndex in rowIndex:blockEnd
        actualNonMajor = rowPermutation[nonMajorRowIndex]
        if actualNonMajor == majorRowIndex
          continue
        end

        #perform multiplication
        temp = mtx.vals[actualNonMajor, rowIndex]
        multiplier = temp / elem
        mtx.vals[actualNonMajor, rowIndex] =  0.0

        lastElem = min(rowIndex + mtx.subMatrixLength * 3 - 1, mtx.size)
        
        #for every element in that row
        for index in (rowIndex + 1):lastElem
          @inbounds currVal = mtx.vals[actualNonMajor, index]
          @inbounds upperVal = mtx.vals[majorRowIndex, index]
          @inbounds mtx.vals[ actualNonMajor, index] = currVal - multiplier * upperVal
        end

        @inbounds bVector[actualNonMajor] -= bVector[majorRowIndex] * multiplier
      end
    end

    return rowPermutation
  end

  function gaussEliminationPivot!(mtx::SparseMatrix, bVector::Vector{Float64})::Vector{Float64}
    swapVector = gaussEliminationPivotPriv!(mtx, bVector)

    for i::UInt64 in length(bVector):-1:1
      actualIndex = swapVector[i]
      for j::UInt64 in (i+1):(min(mtx.size, i + 3 * mtx.subMatrixLength - 1))
        bVector[actualIndex] -= bVector[swapVector[j]] * mtx.vals[actualIndex, j]
      end
      bVector[actualIndex] /= mtx.vals[actualIndex, i]
    end

    result = [bVector[swapVector[i]] for i in 1:(mtx.size)]
    
    return result
  end

  function luDecomposition!(mtx::SparseMatrix)
    #iterate over diagonal
    for rowIndex::UInt64 in 1:(mtx.size - 1)
      elem = mtx.vals[rowIndex,rowIndex]
      @boundscheck if iszero(elem)
        throw(DomainError("diagonal has 0 elem"))
      end

      blockEnd = min(rowIndex + mtx.subMatrixLength, mtx.size)
      #for every row within current block + 1 row
      for lowerRowIndex in (rowIndex + 1):(blockEnd)
        @inbounds temp = mtx.vals[lowerRowIndex,rowIndex]
        @inbounds multiplier = temp / elem
        @inbounds mtx.vals[lowerRowIndex, rowIndex] = multiplier

        lastIndex = min(rowIndex + mtx.subMatrixLength, mtx.size)
        
        #for every element in that row
        for index in (rowIndex + 1):lastIndex
          @inbounds currVal = mtx.vals[lowerRowIndex, index]
          @inbounds upperVal = mtx.vals[rowIndex, index]
          @inbounds mtx.vals[lowerRowIndex, index] = currVal - multiplier * upperVal
        end
      end
    end
  end

  function computeFromLU!(luMtx::SparseMatrix, bVector::Vector{Float64})
    for rowIndex::UInt64 in 2:(luMtx.size)
      # safe version of: rowIndex - luMtx.subMatrixLength > 0 
      columnBegin = rowIndex > luMtx.subMatrixLength ? rowIndex - luMtx.subMatrixLength : 1
      for columnIndex::UInt64 in (rowIndex - 1):-1:columnBegin
        bVector[rowIndex] -= luMtx.vals[rowIndex, columnIndex] * bVector[columnIndex]
      end
    end

    for i::UInt64 in (luMtx.size):-1:1
      for j::UInt64 in (i+1):(min(luMtx.size, i + 2 * luMtx.subMatrixLength - 1))
        @inbounds bVector[i] -= bVector[j] * luMtx.vals[i, j]
      end
      @inbounds bVector[i] /= luMtx.vals[i, i]
    end

    return bVector
  end

  function luDecompositionPivot!(mtx::SparseMatrix)::Vector{UInt64}
    rowPermutation::Vector{UInt64} = map(identity, 1:mtx.size)
    #iterate over diagonal
    for rowIndex in 1:(mtx.size - 1)
      
      #find major element row
      blockEnd = min(rowIndex + mtx.subMatrixLength, mtx.size)
      majorRowIndex::UInt64 = rowIndex
      
      for i in (rowIndex + 1):blockEnd
        firstElem = mtx.vals[rowPermutation[i], rowIndex]
        if abs(firstElem) > abs(mtx.vals[rowPermutation[majorRowIndex], rowIndex])
          majorRowIndex = i
        end
      end

      #update permutation vector
      rowPermutation[rowIndex], rowPermutation[majorRowIndex] = rowPermutation[majorRowIndex], rowPermutation[rowIndex]
      #Now majorRowIndex is the row number in actual representation
      majorRowIndex = rowPermutation[rowIndex]

      elem = mtx.vals[majorRowIndex, rowIndex]

      @boundscheck if iszero(elem)
        throw(DomainError("0 sub-column detected"))
      end

      #for every row within current block + 1 row
      for nonMajorRowIndex in rowIndex:blockEnd
        actualNonMajor = rowPermutation[nonMajorRowIndex]
        if actualNonMajor == majorRowIndex
          continue
        end
                
        #perform multiplication
        temp = mtx.vals[actualNonMajor, rowIndex]
        multiplier = temp / elem
        mtx.vals[actualNonMajor, rowIndex] = multiplier

        lastElem = min(rowIndex + mtx.subMatrixLength * 3 - 1, mtx.size)
        
        #for every element in that row
        for index in (rowIndex + 1):lastElem
          @inbounds currVal = mtx.vals[actualNonMajor, index]
          @inbounds upperVal = mtx.vals[majorRowIndex, index]
          @inbounds mtx.vals[ actualNonMajor, index] = currVal - multiplier * upperVal
        end
      end
    end

    return rowPermutation
  end

  function computeFromLUPivot!(luMtx::SparseMatrix, swapVector::Vector{UInt64}, bVector::Vector{Float64})::Vector{Float64}

    for rowIndex::UInt64 in 2:(luMtx.size)
      # safe version of: rowIndex - 3 * luMtx.subMatrixLength > 0 
      columnBegin = rowIndex > 3 * luMtx.subMatrixLength ? rowIndex - luMtx.subMatrixLength : 1
      
      for columnIndex::UInt64 in (rowIndex - 1):-1:columnBegin
        bVector[swapVector[rowIndex]] -= luMtx.vals[swapVector[rowIndex], columnIndex] * bVector[swapVector[columnIndex]]
      end
    end

    for i::UInt64 in (luMtx.size):-1:1
      actualIndex = swapVector[i]

      for j::UInt64 in (i+1):(min(luMtx.size, i + 3 * luMtx.subMatrixLength - 1))
        bVector[actualIndex] -= bVector[swapVector[j]] * luMtx.vals[actualIndex, j]
      end
      bVector[actualIndex] /= luMtx.vals[actualIndex, i]
    end

    result = [bVector[swapVector[i]] for i in 1:(luMtx.size)]
    
    return result
    
  end

  function readMatrix(filename::String)::SparseMatrix
    open(filename, "r") do file
      t = readline(file)
      fileStr = readlines(file)
      t = split(t)
      size::UInt64 = parse(UInt64, t[1])      
      blockSize::UInt64 = parse(UInt64, t[2])      
      rows, cols, vals = zeros(UInt64, length(fileStr)), zeros(UInt64, length(fileStr)), zeros(Float64, length(fileStr))

      index = 1
      for line in fileStr
        t = split(line)
        rows[index], cols[index], vals[index] = parse(UInt64, t[1]), parse(UInt64, t[2]), parse(Float64, t[3])
        index += 1
      end

      return SparseMatrix(size, blockSize, SparseArrays.sparse(rows, cols, vals))
    end
  end

  function readVector(filename::String)::Vector{Float64}
    open(filename, "r") do file
      t = readline(file)
      length = parse(UInt64, t)
      t = readlines(file)
      result::Vector{Float64} = zeros(length)
      index = 1

      for l in t
        result[index] = parse(Float64, l)
        index += 1
      end
      return result
    end
  end

  function writeVector(filename::String, xVector::Vector{Float64})
    open(filename, "w") do file
      write(file, "$(length(xVector))\n")
      for num in xVector
        write(file, "$num\n")
      end
    end
  end

end
