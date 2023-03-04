module EasyLD

using AWSS3
using ProgressMeter
using PyCall
using CSV
using DataFrames

export download_gnomad_LD_matrices,
    download_ukb_LD_matrices,
    download_gnomad_variant_index_tables,
    download_ukb_variant_index_tables,
    get_gnomad_filenames, 
    get_ukb_filenames, 
    hail_block_matrix,
    read_variant_index_tables

include("download.jl")
include("hailBM.jl")

# This to allow precompilation
# Unlike the @pyimport macro, this does not define a Julia module and members cannot be accessed with s.name.
# @see https://github.com/JuliaPy/PyCall.jl/issues/328
const hail_linalg = PyNULL()
const hail = PyNULL()

function __init__()
    try
        pyimport("hail")
    catch
        error("The EasyLD module is correctly installed, but your python installation is missing the 'hail' module.")
    end
    copy!(hail_linalg, pyimport("hail.linalg"))
    copy!(hail, pyimport("hail"))

    # function for reading a block of BlockMatrix
    py"""
    def convert(bm, x_start, x_end, y_start, y_end):
        return bm[x_start:x_end, y_start:y_end].to_numpy()
    """
end

end # module
