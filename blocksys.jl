#author: Norbert Jaśniewicz
module blocksys

  export SparseMatrix, printMatrix

  mutable struct MatrixRow
    #index of the first element of submatrix A row
    offsetA::UInt64
    #submatrix A row (dense, l-lenght)
    rowA::Vector{Float64}
    #submatrix B row (dense l or 1 length)
    rowB::Vector{Float64}
    #submatrix C row (1 length)
    rowC::Float64
    #relative position from end of A row
    offsetC::UInt64
    MatrixRow(offsetA::UInt64, offsetC::UInt64, rowALength::UInt64, rowBLength::UInt64) =
    rowBLength < 2 || rowBLength == rowALength ?
      new(offsetA, zeros(rowALength), zeros(rowBLength): [0.0], 0.0, offsetC) :
      throw(DomainError("Row b length is not in {0, 1, rowALength}"))
  end

  mutable struct SparseMatrix
    size::UInt64
    subMatrixLength::UInt64
    rows::Vector{MatrixRow}

    function SparseMatrix(size::UInt64, subMatrixLength::UInt64)
      if size % subMatrixLength != 0
        throw(DomainError((size, subMatrixLength), "size mod subMatrixLength != 0 is not allowed"))
      elseif size == subMatrixLength
        throw(DomainError((size, subMatrixLength), "subMatrixLength == size is not allowed"))
      end
      subSize::UInt64 = div(size, subMatrixLength)
      index::UInt64 = 0
      rows::Vector{MatrixRow} = []

      push!(rows, MatrixRow(1 + index * subMatrixLength, UInt64(1), subMatrixLength, UInt64(0)))
      for i in 2:subMatrixLength
        push!(rows, MatrixRow(1 + index * subMatrixLength, i, subMatrixLength, UInt64(1)))
      end
      index += 1

      for _ in 2:subSize
        push!(rows, MatrixRow(1 + index * subMatrixLength, UInt64(1), subMatrixLength, subMatrixLength))
        for i in 2:subMatrixLength
          push!(rows, MatrixRow(1 + index * subMatrixLength, i, subMatrixLength, UInt64(1)))
        end
        index += 1
      end

      return new(size, subMatrixLength, rows)
    end
  end

  function printMatrix(mtx::SparseMatrix)    
    for row in mtx.rows
      #prefix::String = "* " ^ (row.offsetA - length(row.rowB))
      println(row.offsetA)
      println(length(row.rowB))
      #println(prefix)
      #bRow = map(x -> parse(UInt64, x), row.rowB)
      #bRow = join(bRow, " ")
      #aRow = map(x -> parse(UInt64, x), row.rowA)
      #aRow = join(aRow, " ")
      #interfix = "* " ^ (row.offsetA - 1)
      #cRow = parse(UInt64, row.rowC)
      #postfix = "#" ^ (mtx.size - (row.offsetA + row.offsetC))
      #println(prefix + bRow + aRow + interfix + cRow + postfix)
    end
  end

end