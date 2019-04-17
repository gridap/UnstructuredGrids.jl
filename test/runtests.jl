module UnstructuredGridsTests

using Test

@testset "UnstructuredGrids.jl" begin

  @testset "Kernels" begin include("KernelsTests.jl") end

end

end #module UnstructuredGridsTests
