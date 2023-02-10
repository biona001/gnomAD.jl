"""
The Hail Block Matrix format: https://hail.is/docs/0.2/_modules/hail/linalg/blockmatrix.html#BlockMatrix

By default, a gnomAD LD matrix is a banded matrix with band width 1 and blocks of 
size 4096 by 4096. The blocks look like this:
"""
struct HailBlockMatrix <: AbstractMatrix{Float64}
    bm_files::String
    bm::PyObject
end

Base.getindex(x::HailBlockMatrix, i::Int) = 
    error("Can only access HailBlockMatrix using 2 indices (i, j)")
function Base.getindex(x::HailBlockMatrix, i::Int, j::Int)
    i == 0 || j == 0 && error("attempt to access matrix at index [$i, $j]")
    i > j ? 0.0 : x.bm[i, j]
end
Base.size(x::HailBlockMatrix, k::Int) = 
    k == 1 ? x.bm.n_rows : k == 2 ? x.bm.n_cols : k > 2 ? 1 : error("Dimension k out of range")
Base.size(x::HailBlockMatrix) = (size(x, 1), size(x, 2))

function hail_block_matrix(bm_files::String)
    isdir(bm_files) || error("Directory $bm_files does not exist")
    return HailBlockMatrix(bm_files, hail.BlockMatrix.read(bm_files))
end
