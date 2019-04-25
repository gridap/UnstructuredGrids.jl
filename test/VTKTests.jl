module VTKTests

using Test
using UnstructuredGrids.VTK
using UnstructuredGrids.Core
using UnstructuredGrids.Factories
using UnstructuredGrids.Factories: HEXAHEDRON

d = mktempdir()
f = joinpath(d,"grid")

grid = Grid(domain=(0,1,-1,0),partition=(2,2))

writevtk(grid,f)

grid = Grid(domain=(0,1,-1,0,2,4),partition=(2,3,4))

writevtk(grid,f)

fgrid = Grid(grid,dim=2)

writevtk(fgrid,f)

fgrid = Grid(grid,dim=1)

writevtk(fgrid,f)

writevtk(HEXAHEDRON,"hex")

rm(d,recursive=true)

end # module VTKTests
