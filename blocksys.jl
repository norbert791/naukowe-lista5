#author: Norbert Jaśniewicz
module blocksys

  export SparseMatrix, printMatrix, setCell!, getCell

  @enum MatrixRowType begin
    NO_B
    SHORT_B
    LONG_B
  end

  mutable struct MatrixRow
    #index of first element in A row
    offsetA :: UInt64
    #index of first element in B row, 0 means no B row
    offsetB :: UInt64
    #index of C, 0 means no C val
    offsetC :: UInt64
    abValues :: Vector{Float64}
    cVal :: Float64
    function MatrixRow(offsetA::UInt64, relOffsetC::UInt64, length::UInt64, type::MatrixRowType)
      bLength::UInt64 = 0
      if type == NO_B
        bLength = 0
      elseif type == SHORT_B
        bLength = 1
      else
        bLength = length
      end
      abValues::Vector{Float64} = zeros(length + bLength)
      offsetB::UInt64 = bLength == 0 ? 0 : offsetA - bLength

      offsetC::UInt64 = relOffsetC == 0 ? 0 : offsetA + length + (relOffsetC - 1)
      return new(offsetA, offsetB, offsetC, abValues, UInt64(0))
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
      for i in 1:subMatrixLength
        push!(rows, MatrixRow(UInt64(1), UInt64(i), subMatrixLength, NO_B))
      end

      #Number of remaining subMatrix rows to insert
      blockNum :: UInt64 = (div(size, subMatrixLength)) - 1

      #Insert all B_iA_iC_i block rows except the last one
      for index in 1:(blockNum - 1)
        offsetA::UInt64 = 1 + index * subMatrixLength
        push!(rows, MatrixRow(offsetA, UInt64(1), subMatrixLength, LONG_B))
        
        for i in 2:subMatrixLength
          push!(rows, MatrixRow(offsetA, UInt64(i), subMatrixLength, SHORT_B))
        end
      end

      #Insert the last block row: B_vA_v
      push!(rows, MatrixRow(1 + blockNum * subMatrixLength, UInt64(0), subMatrixLength, LONG_B))

      for _ in 2:subMatrixLength
        push!(rows, MatrixRow(1 + blockNum * subMatrixLength, UInt64(0), subMatrixLength, SHORT_B))
      end

      return new(size, subMatrixLength, rows)
    end
    
  end

  function printMatrix(mtx::SparseMatrix)
    #for row in mtx.rows
    #  printRow(row)
    #end
    for row in mtx.rows
      prefix = "* " ^ (row.offsetB == 0 ? 0 : row.offsetB - 1)
      vals = map(x -> string(x), row.abValues)
      vals = join(vals, " ")
      interfix = " " * "* " ^ (row.offsetC == 0 ? 0 : row.offsetC - row.offsetA - mtx.subMatrixLength)
      interfix = interfix == "" ? " " : interfix
      interfix *= row.offsetC == 0 ? "" : string(row.cVal) * " "
      suffix = row.offsetC == 0 ? "" : "* " ^ (mtx.size - row.offsetC)
      result = prefix * vals * interfix * suffix * "\n"
      println(result)
    end
  end

  function setCell!(mtx::SparseMatrix, row::UInt64, column::UInt64, value::Float64)
    if column == 0 || column > mtx.size
      throw(BoundsError(column, "Index of out bound"))
    end
    offsetA::UInt64 = mtx.rows[row].offsetA
    offsetB::UInt64 = mtx.rows[row].offsetB
    offsetC::UInt64 = mtx.rows[row].offsetC
    if column <= offsetA + mtx.subMatrixLength - 1 && column >= offsetB
      mtx.rows[row].abValues[column - offsetB + 1] = value
    elseif column == offsetC
      mtx.rows[row].cVal = value
    else
      throw(BoundsError(column, "Attempt to assign value to 0-field"))
    end
  end

  function getCell(mtx::SparseMatrix, row::UInt64, column::UInt64)::Float64
    if column == 0 || column > mtx.size
      throw(BoundsError(column, "Index of out bound"))
    end
    offsetA::UInt64 = mtx.rows[row].offsetA
    offsetB::UInt64 = mtx.rows[row].offsetB
    offsetC::UInt64 = mtx.rows[row].offsetC
    if column <= offsetA + mtx.subMatrixLength - 1 && column >= offsetB
      return mtx.rows[row].abValues[column - offsetB + 1]
    elseif column == offsetC
      return mtx.rows[row].cVal
    else
      return 0.0
    end
  end
end
