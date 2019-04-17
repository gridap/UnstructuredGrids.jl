module UnstructuredGrids

include("Helpers.jl")
include("Kernels.jl")
include("Core.jl")
include("Factories.jl")

export UGrid, gridpoints, gridcells, celltypes
using UnstructuredGrids.Core: UGrid, gridpoints, gridcells, celltypes

end # module UnstructuredGrids
