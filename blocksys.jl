#author: Norbert Jaśniewicz
module blocksys

  export SparseMatrix, printMatrix, setCell!, getCell, gaussElimination!, readMatrix, readVector

  @enum MatrixRowType begin
    NO_B
    SHORT_B
    LONG_B
  end

  mutable struct MatrixRow
    #index of first element in A row
    offsetA::UInt64
    #index of first element in B row, 0 means no B row
    offsetB::UInt64
    #index of first element in C row, 0 means no C row
    offsetC::UInt64
    values::Vector{Float64}

    function MatrixRow(offsetA::UInt64, length::UInt64, bType::MatrixRowType, cRowPresent::Bool)
      bLength::UInt64 = 0
      cLength::UInt64 = cRowPresent * length
      if bType == NO_B
        bLength = 0
      elseif bType == SHORT_B
        bLength = 1
      else
        bLength = length
      end
      
      values::Vector{Float64} = zeros(length + bLength + cLength)
      offsetB::UInt64 = bLength == 0 ? 0 : offsetA - bLength
      offsetC::UInt64 = (offsetA + length) * cRowPresent

      return new(offsetA, offsetB, offsetC, values)
    end
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
      for _ in 1:subMatrixLength
        push!(rows, MatrixRow(UInt64(1), subMatrixLength, NO_B, true))
      end

      #Number of remaining subMatrix rows to insert
      blockNum :: UInt64 = (div(size, subMatrixLength)) - 1

      #Insert all B_iA_iC_i block rows except the last one
      for index in 1:(blockNum - 1)
        offsetA::UInt64 = 1 + index * subMatrixLength
        push!(rows, MatrixRow(offsetA, subMatrixLength, LONG_B, true))
        
        for _ in 2:subMatrixLength
          push!(rows, MatrixRow(offsetA, subMatrixLength, SHORT_B, true))
        end
      end

      #Insert the last block row: B_vA_v
      push!(rows, MatrixRow(1 + blockNum * subMatrixLength, subMatrixLength, LONG_B, false))

      for _ in 2:subMatrixLength
        push!(rows, MatrixRow(1 + blockNum * subMatrixLength, subMatrixLength, SHORT_B, false))
      end

      return new(size, subMatrixLength, rows)
    end
    
  end

  function printMatrix(mtx::SparseMatrix)
    #for row in mtx.rows
    #  println(row)
    #end

    for row in mtx.rows
      prefix = "* " ^ (row.offsetB == 0 ? 0 : row.offsetB - 1)
      vals = map(x -> string(x), row.values)
      vals = join(vals, " ")
      suffix = row.offsetC == 0 ? "" : "* " ^ (mtx.size + 1 - row.offsetC - mtx.subMatrixLength)
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

    offsetA::UInt64 = mtx.rows[row].offsetA
    offsetB::UInt64 = mtx.rows[row].offsetB
    offsetC::UInt64 = mtx.rows[row].offsetC
    leftOffset::UInt64 = offsetB == 0 ? offsetA : offsetB
    rightOffset::UInt64 = (offsetC == 0 ? offsetA : offsetC) +
                           mtx.subMatrixLength - 1
    
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

    offsetA::UInt64 = mtx.rows[row].offsetA
    offsetB::UInt64 = mtx.rows[row].offsetB
    offsetC::UInt64 = mtx.rows[row].offsetC
    leftOffset::UInt64 = offsetB == 0 ? offsetA : offsetB
    rightOffset::UInt64 = (offsetC == 0 ? offsetA : offsetC) +
                           mtx.subMatrixLength - 1
    
    if leftOffset <= column && column <= rightOffset
      return mtx.rows[row].values[column - leftOffset + 1]
    else
      return 0.0
    end
  end

  @inline function setCellNoCheck!(mtx::SparseMatrix, row::UInt64, column::UInt64, value::Float64)
    offsetA::UInt64 = mtx.rows[row].offsetA
    offsetB::UInt64 = mtx.rows[row].offsetB
    leftOffset::UInt64 = offsetB == 0 ? offsetA : offsetB
    (mtx.rows[row].values[column - leftOffset + 1] = value)
  end

  @inline function getCellNoCheck(mtx::SparseMatrix, row::UInt64, column::UInt64)::Float64
    offsetA::UInt64 = mtx.rows[row].offsetA
    offsetB::UInt64 = mtx.rows[row].offsetB
    leftOffset::UInt64 = offsetB == 0 ? offsetA : offsetB
    (return mtx.rows[row].values[column - leftOffset + 1])
  end

  function gaussElimination!(mtx::SparseMatrix, bVector::Vector{Float64})
    @boundscheck if mtx.size != length(bVector)
      throw(DomainError("bVector and mtx have different length"))
    end

    #iterate over diagonal
    for rowIndex in 1:(mtx.size - 1)
      elem = getCellNoCheck(mtx, rowIndex, rowIndex)
      #printMatrix(mtx)
      #println("-------------------a-----------------------")
      @boundscheck if iszero(elem)
        #printMatrix(mtx)
        throw(DomainError("diagonal has 0 elem"))
      end

      #for every row within current block
      for lRowIndex in (rowIndex + 1):(rowIndex + mtx.subMatrixLength)
        println((lRowIndex, rowIndex))
        temp = getCellNoCheck(mtx, lRowIndex, rowIndex)
        multiplier = temp / elem
        setCellNoCheck!(mtx, lRowIndex, rowIndex, 0.0)

        lastElem = mtx.subMatrixLength + (mtx.rows[rowIndex].offsetC == 0) * 
          mtx.subMatrixLength - 1
        
        #for every element in that row
        for index in (rowIndex + 1):(rowIndex + lastElem)
          println((lRowIndex, index))
          currVal = getCellNoCheck(mtx, lRowIndex, index)
          println((lRowIndex, index))
          upperVal = getCellNoCheck(mtx, rowIndex, index)
          setCellNoCheck!(mtx, lRowIndex, index, currVal - multiplier * upperVal)
        end
      end
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
