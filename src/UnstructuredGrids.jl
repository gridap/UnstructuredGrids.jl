module UnstructuredGrids

include("Helpers.jl")
include("Kernels.jl")
include("Core.jl")
include("CoreNew.jl")
include("Factories.jl")
include("VTK.jl")

export UGrid, RefCell, VERTEX, GridGraph
export gridpoints, gridcells, celltypes, refcells
export connections
export writevtk

using UnstructuredGrids.Core: UGrid, RefCell, VERTEX, GridGraph
using UnstructuredGrids.Core: gridpoints, gridcells, celltypes, refcells
using UnstructuredGrids.Core: connections
using UnstructuredGrids.VTK: writevtk

end # module UnstructuredGrids
