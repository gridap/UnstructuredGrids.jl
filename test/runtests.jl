module UnstructuredGridsTests

using Test

@testset "UnstructuredGrids.jl" begin

  @testset "Kernels" begin include("KernelsTests.jl") end
  @testset "Factories" begin include("FactoriesTests.jl") end

end

end #module UnstructuredGridsTests
