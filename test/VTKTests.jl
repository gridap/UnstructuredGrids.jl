module VTKTests

using Test
using UnstructuredGrids
using UnstructuredGrids.Factories

grid = generate(domain=(0,1,-1,0),partition=(2,2))

d = mktempdir()
f = joinpath(d,"grid")

writevtk(grid,f)

rm(d,recursive=true)

end # module VTKTests
