#author: Norbert Jaśniewicz
module blocksys

  export SparseMatrix, printMatrix, setCell!, getCell, gaussElimination!, gausseliminationMajor!, readMatrix, readVector

  mutable struct MatrixRow
    #actual first index
    firstIndex::UInt64
    #actual last index
    lastIndex::UInt64
    #values
    values::Vector{Float64}

    MatrixRow(firstIndex, lastIndex, ALength) = new(firstIndex, lastIndex, zeros(4 * ALength))
  end

  mutable struct SparseMatrix
    size::UInt64
    subMatrixLength::UInt64
    rows::Vector{MatrixRow}

    function SparseMatrix(size::UInt64, subMatrixLength::UInt64)
      if (size % subMatrixLength != 0)
        throw(DomainError((size, subMatrixLength), "size mod subMatrixLength must be equal 0"))
      end
      if (size == subMatrixLength)
        throw(DomainError((size, subMatrixLength), "size == subMatrixLength"))
      end
      
      rows::Vector{MatrixRow} = []
      
      #Insert A_0C_0 rows
      for i in 1:subMatrixLength
        push!(rows, MatrixRow(1, subMatrixLength + i, subMatrixLength))
      end

      #Number of remaining blocks
      numBlocks = div(size, subMatrixLength) - 1

      for i in 1:(numBlocks-1) #first block already inserted, last one ill be inserted manually
        offsetA = i * subMatrixLength + 1
        push!(rows, MatrixRow(offsetA - subMatrixLength, offsetA + subMatrixLength, subMatrixLength))
        
        for j in 1:(subMatrixLength - 1)
          push!(rows, MatrixRow(offsetA - 1, offsetA + subMatrixLength + j, subMatrixLength))
        end
      end

      #BA block
      offsetA = numBlocks * subMatrixLength + 1
      push!(rows, MatrixRow(offsetA - subMatrixLength, offsetA + subMatrixLength - 1, subMatrixLength))

      for _ in 1:(subMatrixLength - 1)
        push!(rows, MatrixRow(offsetA - 1, offsetA + subMatrixLength - 1, subMatrixLength))
      end

      return new(size, subMatrixLength, rows)
    end
    
  end

  function printMatrix(mtx::SparseMatrix)

    for row in mtx.rows
      prefix = "* " ^ (row.firstIndex == 0 ? 0 : row.firstIndex - 1)
      length = (row.lastIndex - row.firstIndex)
      rowView = row.values[1:(length + 1)]
      vals = map(x -> string(x), rowView)
      vals = join(vals, " ")
      suffix = "* " ^ (mtx.size - row.lastIndex)
      result = prefix * vals * " " * suffix * "\n"
      println(result)
    end
  end

  function setCell!(mtx::SparseMatrix, row::UInt64, column::UInt64, value::Float64)
    @boundscheck if row > mtx.size || row == 0
      throw(BoundsError(row, "Index of out bound"))
    end

    @boundscheck if column == 0 || column > mtx.size
      throw(BoundsError(row, "Index of out bound"))
    end
    
    leftOffset::UInt64 = mtx.rows[row].firstIndex
    rightOffset::UInt64 = leftOffset + 3 * mtx.subMatrixLength - 1

    if leftOffset <= column && column <= rightOffset
      mtx.rows[row].values[column - leftOffset + 1] = value
    else
      throw(BoundsError(row, "Attempt to assign value to 0-field"))
    end
  end

  function getCell(mtx::SparseMatrix, row::UInt64, column::UInt64)::Float64
    @boundscheck if row > mtx.size || row == 0
      throw(BoundsError(row, "Index of out bound"))
    end
    
    
    @boundscheck if column == 0 || column > mtx.size
      throw(BoundsError(row, "Index of out bound"))
    end

    leftOffset::UInt64 = mtx.rows[row].firstIndex
    rightOffset::UInt64 = leftOffset + 4 * mtx.subMatrixLength - 1
    
    if leftOffset <= column && column <= rightOffset
      return mtx.rows[row].values[column - leftOffset + 1]
    else
      return 0.0
    end
  end

  @inline function setCellNoCheck!(mtx::SparseMatrix, row::UInt64, column::UInt64, value::Float64)
    leftOffset::UInt64 = mtx.rows[row].firstIndex
    @inbounds (mtx.rows[row].values[column - leftOffset + 1] = value)
  end

  @inline function getCellNoCheck(mtx::SparseMatrix, row::UInt64, column::UInt64)::Float64
    leftOffset::UInt64 = mtx.rows[row].firstIndex
    return @inbounds (mtx.rows[row].values[column - leftOffset + 1])
  end

  function gaussEliminationPriv!(mtx::SparseMatrix, bVector::Vector{Float64})
    @boundscheck if mtx.size != length(bVector)
      throw(DomainError("bVector and mtx have different length"))
    end

    #iterate over diagonal
    for rowIndex in 1:(mtx.size - 1)
      elem = getCellNoCheck(mtx, rowIndex, rowIndex)
      @boundscheck if iszero(elem)
        throw(DomainError("diagonal has 0 elem"))
      end

      blockEnd = min(rowIndex + mtx.subMatrixLength, mtx.size)
      #for every row within current block + 1 row
      for lowerRowIndex in (rowIndex + 1):(blockEnd)
        temp = getCellNoCheck(mtx, lowerRowIndex, rowIndex)
        multiplier = temp / elem
        setCellNoCheck!(mtx, lowerRowIndex, rowIndex, 0.0)

        lastElem = mtx.rows[rowIndex].lastIndex
        
        #for every element in that row
        for index in (rowIndex + 1):lastElem
          currVal = getCellNoCheck(mtx, lowerRowIndex, index)
          upperVal = getCellNoCheck(mtx, rowIndex, index)
          setCellNoCheck!(mtx, lowerRowIndex, index, currVal - multiplier * upperVal)
        end
        @inbounds bVector[lowerRowIndex] -= bVector[rowIndex] * multiplier
      end
    end

  end

  function gaussElimination!(mtx::SparseMatrix, bVector::Vector{Float64})::Vector{Float64}
    gaussEliminationPriv!(mtx::SparseMatrix, bVector::Vector{Float64})
    for i::UInt64 in length(bVector):-1:1
      for j::UInt64 in (i+1):(min(mtx.size, i + 2 * mtx.subMatrixLength - 1))
        bVector[i] -= bVector[j] * getCellNoCheck(mtx, i, j)
      end
      bVector[i] /= getCellNoCheck(mtx, i, i)
    end

    return bVector
  end

  function gaussEliminationMajorPriv!(mtx::SparseMatrix, bVector::Vector{Float64})::Vector{UInt64}
    @boundscheck if mtx.size != length(bVector)
      throw(DomainError("bVector and mtx have different length"))
    end
    rowPermutation::Vector{UInt64} = map(identity, 1:mtx.length)
    #iterate over diagonal
    for rowIndex in 1:(mtx.size - 1)
      
      #find major element row
      blockEnd = min(rowIndex + mtx.subMatrixLength, mtx.size)
      majorRowIndex::UInt64 = rowIndex

      for i in (rowIndex + 1):blockEnd
        firstElem = getCell(mtx, i, rowPermutation[rowIndex])
        if abs(firstElem) > abs(mtx.rows[rowPermutation[majorRowIndex]].values[rowIndex])
          majorRowIndex = i
        end
      end

      
      majorRowView = mtx.rows[rowPermutation[majorRowIndex]]
      elem = getCellNoCheck(mtx, rowPermutation[majorRowIndex], rowIndex)
      #update permutation vector
      rowPermutation[rowIndex], rowPermutation[rowPermutation[majorRowIndex]] = majorRowIndex, rowIndex

      @boundscheck if iszero(elem)
        throw(DomainError("0 sub-column detected"))
      end

      #for every row within current block + 1 row
      for nonMajorRowIndex in rowIndex:blockEnd
        if nonMajorRowIndex == majorRowIndex
          continue
        end
        #align non major row to major row
        if (majorRowIndex == blockEnd)
          mtx.rows[nonMajorRowIndex].firstIndex = majorRowView.firstIndex
          mtx.rows[nonMajorRowIndex].lastIndex = majorRowView.lastIndex
        end

        temp = getCellNoCheck(mtx, nonMajorRowIndex, rowIndex)
        multiplier = temp / elem
        setCellNoCheck!(mtx, nonMajorRowIndex, rowIndex, 0.0)

        lastElem = majorRowView.lastIndex
        
        #for every element in that row
        for index in (rowIndex + 1):lastElem
          currVal = getCellNoCheck(mtx, nonMajorRowIndex, index)
          upperVal = getCellNoCheck(mtx, majorRowIndex, index)
          setCellNoCheck!(mtx, nonMajorRowIndex, index, currVal - multiplier * upperVal)
        end
        @inbounds bVector[nonMajorRowIndex] -= bVector[rowIndex] * multiplier
      end
    end

    return rowPermutation
  end

  function gaussEliminationMajor!(mtx::SparseMatrix, bVector::Vector{Float64})::Vector{Float64}
    swapVector = gaussEliminationMajorPriv!(mtx, bVector)
    for i::UInt64 in length(bVector):-1:1
      actualIndex = swapVector[i]
      for j::UInt64 in (i+1):(min(mtx.size, i + 2 * mtx.subMatrixLength - 1))
        bVector[actualIndex] -= bVector[j] * getCellNoCheck(mtx, actualIndex, j)
      end
      bVector[actualIndex] /= getCellNoCheck(mtx, actualIndex, i)
    end

    return bVector
  end

  function shiftLeft!(row::MatrixRow, shiftLength::UInt64)
    row.firstIndex += shiftLength
    
    for i in 1:(length(row.values) - shiftLength)
      row.values[i] = row.values[i + shiftLength]
    end
  end

  function readMatrix(filename::String)::SparseMatrix
    open(filename, "r") do file
      t = readline(file)
      fileStr = readlines(file)
      t = split(t)
      size::UInt64 = parse(UInt64, t[1])      
      blockSize::UInt64 = parse(UInt64, t[2])      
      result = SparseMatrix(size, blockSize)
      
      for line in fileStr
        t = split(line)
        row, col, val = parse(UInt64, t[1]), parse(UInt64, t[2]), parse(Float64, t[3])
        setCell!(result, row, col, val)
      end
      return result
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
end
