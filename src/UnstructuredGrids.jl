module UnstructuredGrids

include("Helpers.jl")
include("Kernels.jl")
include("Core.jl")
include("Factories.jl")
include("VTK.jl")

export Connections
export RefCell
export UGrid
export VERTEX
export connections
export coordinates
export writevtk
export list, ptrs
export generate_dual_connections
export generate_cell_to_faces
export find_cell_to_faces

using UnstructuredGrids.Core: Connections
using UnstructuredGrids.Core: UGrid
using UnstructuredGrids.Core: RefCell
using UnstructuredGrids.Core: list, ptrs
using UnstructuredGrids.Core: connections
using UnstructuredGrids.Core: coordinates
using UnstructuredGrids.Core: VERTEX
using UnstructuredGrids.VTK: writevtk
using UnstructuredGrids.Kernels: generate_dual_connections
using UnstructuredGrids.Kernels: generate_cell_to_faces
using UnstructuredGrids.Kernels: find_cell_to_faces

end # module UnstructuredGrids
