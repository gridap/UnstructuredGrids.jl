module UnstructuredGrids

include("Helpers.jl")
include("Kernels.jl")
include("Core.jl")
include("Factories.jl")
include("VTK.jl")

export Connections
export Grid
export connections
export coordinates
export writevtk
export list, ptrs
export generate_face_to_cells
export generate_cell_to_faces

using UnstructuredGrids.Core: Connections
using UnstructuredGrids.Core: Grid
using UnstructuredGrids.Core: list, ptrs
using UnstructuredGrids.Core: connections
using UnstructuredGrids.Core: coordinates
using UnstructuredGrids.VTK: writevtk
using UnstructuredGrids.Kernels: generate_face_to_cells
using UnstructuredGrids.Kernels: generate_cell_to_faces

end # module UnstructuredGrids
