module gnomAD

using AWSS3
using ProgressMeter
using PyCall
using LinearAlgebra

export download_LD_matrices,
    get_all_filenames,
    hail_block_matrix,
    get_block

include("download.jl")
include("hailBM.jl")

# This to allow precompilation
# Unlike the @pyimport macro, this does not define a Julia module and members cannot be accessed with s.name.
# @see https://github.com/JuliaPy/PyCall.jl/issues/328
const hail = PyNULL()
const numpy = PyNULL()

function __init__()
    try
        pyimport("hail")
    catch
        error("The gnomAD module is correctly installed, but your python installation is missing the 'hail' module.")
    end
    try
        pyimport("numpy")
    catch
        error("The gnomAD module is correctly installed, but your python installation is missing the 'numpy' module.")
    end
    copy!(hail, pyimport("hail.linalg"))
    copy!(numpy, pyimport("numpy"))
    
    # function for reading a block of BlockMatrix
    py"""
    def convert(bm, range_start, range_end):
        return bm[range_start:range_end, range_start:range_end].to_numpy()
    """
end


end # module
