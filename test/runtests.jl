using EasyLD
using Test

@testset "gnomAD download and read" begin
    data = joinpath(
        dirname(pathof(EasyLD)), 
        "../", 
        "data/gnomad.genomes.r2.1.1.nfe.common.adj.ld.bm"
    )

    # download 1 file
    file1 = joinpath(data, "parts", 
        "part-00000-14-0-0-7e282c91-9060-db76-1d4a-a02877b62910")
    if !isfile(file1)
        mkpath(joinpath(data, "parts"))
        population = "nfe"
        download_gnomad_LD_matrices(population, data, start_from=1, num_files=1)
    end
    @test isfile(file1)
    @test isfile(joinpath(data, "_SUCCESS"))
    @test isfile(joinpath(data, "metadata.json"))

    # basic abstract matrix interface
    bm = hail_block_matrix(data);
    @test size(bm) == (14207204, 14207204)
    @test bm[1, 1] ≈ 1
    @test bm[1, 2] ≈ 0.032636100824383195
    @time Sigma = bm[1:1000, 1:1000]
    @test Sigma[1, 1] ≈ bm[1, 1]
    @test Sigma[100, 500] ≈ bm[100, 500]
end

@testset "pan-ukbb download and read" begin
    data = joinpath(
        dirname(pathof(EasyLD)), 
        "../", 
        "data/UKBB.EUR.ldadj.bm"
    )

    # download 1 file
    file1 = joinpath(data, "parts", 
        "part-000000-44-0-0-22f828c8-17c3-7c3c-1fa1-1fc113144aca")
    if !isfile(file1)
        mkpath(joinpath(data, "parts"))
        population = "EUR"
        download_ukb_LD_matrices(population, data, start_from=1, num_files=1)
    end
    @test isfile(file1)
    @test isfile(joinpath(data, "_SUCCESS"))
    @test isfile(joinpath(data, "metadata.json"))

    # basic abstract matrix interface
    bm = hail_block_matrix(data);
    @test size(bm) == (23960350, 23960350)
    @test bm[1, 1] ≈ 0.999959841549239
    @test bm[1, 2] ≈ 0.0009855904101964482
    @time Sigma = bm[1:1000, 1:1000]
    @test Sigma[1, 1] ≈ bm[1, 1]
    @test Sigma[100, 500] ≈ bm[100, 500]
end
