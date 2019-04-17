module FactoriesTests

using Test
using UnstructuredGrids
using UnstructuredGrids.Factories

grid = generate(domain=(0,1,-1,0),partition=(2,2))

@test isa(grid,UGrid)

p = [0.0 1.0 2.0 0.0 1.0 2.0 0.0 1.0 2.0;
     -1.0 -1.0 -1.0 -0.5 -0.5 -0.5 0.0 0.0 0.0] 

c = ([1, 2, 4, 5, 2, 3, 5, 6, 4, 5, 7, 8, 5, 6, 8, 9], [1, 5, 9, 13, 17]) 

@test gridpoints(grid) ≈ p

@test gridcells(grid) == c

end # module FactoriesTests
