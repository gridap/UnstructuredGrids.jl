module UnstructuredGrids

include("Helpers.jl")
include("Kernels.jl")
include("Core.jl")
include("Factories.jl")
include("VTK.jl")

export UGrid, RefCell, VERTEX
export gridpoints, gridcells, celltypes, refcells
export writevtk

using UnstructuredGrids.Core: UGrid, RefCell, VERTEX
using UnstructuredGrids.Core: gridpoints, gridcells, celltypes, refcells
using UnstructuredGrids.VTK: writevtk

end # module UnstructuredGrids
