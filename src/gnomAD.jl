module gnomAD

using AWSS3
using ProgressMeter
using PyCall

export download_gnomad_LD_matrices,
    download_ukb_LD_matrices,
    get_gnomad_filenames, 
    get_ukb_filenames, 
    hail_block_matrix

include("download.jl")
include("hailBM.jl")

# This to allow precompilation
# Unlike the @pyimport macro, this does not define a Julia module and members cannot be accessed with s.name.
# @see https://github.com/JuliaPy/PyCall.jl/issues/328
const hail = PyNULL()

function __init__()
    try
        pyimport("hail")
    catch
        error("The gnomAD module is correctly installed, but your python installation is missing the 'hail' module.")
    end
    copy!(hail, pyimport("hail.linalg"))
    
    # function for reading a block of BlockMatrix
    py"""
    def convert(bm, x_start, x_end, y_start, y_end):
        return bm[x_start:x_end, y_start:y_end].to_numpy()
    """
end


end # module
