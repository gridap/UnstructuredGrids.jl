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

@test celltypes(grid) == ones(Int,length(c[2])-1)

grid = generate(domain=(0,1,-1,0,2,4),partition=(1,2,3))

p = [0.0 1.0 0.0 1.0 0.0 1.0 0.0 1.0 0.0 1.0 0.0 1.0 0.0 1.0 0.0 1.0 0.0 1.0 0.0 1.0 0.0 1.0 0.0 1.0;
     -1.0 -1.0 -0.5 -0.5 0.0 0.0 -1.0 -1.0 -0.5 -0.5 0.0 0.0 -1.0 -1.0 -0.5 -0.5 0.0 0.0 -1.0 -1.0 -0.5 -0.5 0.0 0.0;
     2.0 2.0 2.0 2.0 2.0 2.0 2.5 2.5 2.5 2.5 2.5 2.5 3.0 3.0 3.0 3.0 3.0 3.0 3.5 3.5 3.5 3.5 3.5 3.5]

c = ([1, 2, 3, 4, 7, 8, 9, 10, 3, 4, 5, 6, 9, 10, 11, 12, 7, 8, 9, 10, 13, 14, 15, 16, 9, 10, 11, 12,
      15, 16, 17, 18, 13, 14, 15, 16, 19, 20, 21, 22, 15, 16, 17, 18, 21, 22, 23, 24],
     [1, 9, 17, 25, 33, 41, 49])

@test gridpoints(grid) ≈ p

@test gridcells(grid) == c

@test celltypes(grid) == ones(Int,length(c[2])-1)

end # module FactoriesTests
