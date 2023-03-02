"""
Currently a `HailBlockMatrix` is just a python wrapper, where data reading is 
achieved within python by the Hail software, then the result is passed into
Julia via PyCall.jl. The `bm` field behaves as a `BlockMatrix` from hail. 
In the future, a native Julia reader that directly reads
the compreseed format is desired. 

The Hail Block Matrix format is described here
https://hail.is/docs/0.2/_modules/hail/linalg/blockmatrix.html#BlockMatrix

A more specific file description is provided here
https://discuss.hail.is/t/blockmatrix-specification/3118
"""
struct HailBlockMatrix <: AbstractMatrix{Float64}
    bm_files::String
    bm::PyObject
end

# do not allow getindex for a single dimension
Base.getindex(x::HailBlockMatrix, i::Int) = 
    error("Can only access HailBlockMatrix using 2 indices (i, j)")
Base.getindex(x::HailBlockMatrix, ::AbstractRange) = 
    error("Can only access HailBlockMatrix using 2 dimension index, e.g. bm[1:10, 1:10]")

function Base.getindex(x::HailBlockMatrix, i::Int, j::Int)
    i ≤ 0 || j ≤ 0 && error("attempt to access matrix at index [$i, $j]")
    i > j ? 0.0 : x.bm[i, j]
end
function Base.getindex(x::HailBlockMatrix, idx1::UnitRange, idx2::UnitRange)
    _get_block(x, first(idx1), last(idx1), first(idx2), last(idx2))
end
function _get_block(x::HailBlockMatrix, 
    x_start::Int, x_end::Int, y_start::Int, y_end::Int)
    1 ≤ x_start ≤ size(x, 1) && 1 ≤ y_start ≤ size(x, 1) && 
        1 ≤ x_end ≤ size(x, 1) && 1 ≤ y_end ≤ size(x, 1) || 
        error("x_start out of bounds")
    # convert function uses python indexing, i.e. 0 based
    return py"convert"(x.bm, x_start - 1, x_end, y_start - 1, y_end)
end

Base.size(x::HailBlockMatrix, k::Int) = 
    k == 1 ? x.bm.n_rows : k == 2 ? x.bm.n_cols : k > 2 ? 1 : error("Dimension k out of range")
Base.size(x::HailBlockMatrix) = (size(x, 1), size(x, 2))

function hail_block_matrix(bm_files::String)
    isdir(bm_files) || error("Directory $bm_files does not exist")
    return HailBlockMatrix(bm_files, hail.BlockMatrix.read(bm_files))
end
