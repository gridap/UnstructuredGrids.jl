module VTKTests

using Test
using UnstructuredGrids
using UnstructuredGrids.Factories

d = mktempdir()
f = joinpath(d,"grid")

grid = generate(domain=(0,1,-1,0),partition=(2,2))

writevtk(grid,f)

grid = generate(domain=(0,1,-1,0,2,4),partition=(2,3,4))

writevtk(grid,f)

rm(d,recursive=true)

end # module VTKTests
