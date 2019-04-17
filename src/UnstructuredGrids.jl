module UnstructuredGrids

export UnstructuredGrid
export Connections

include("Kernels.jl")

struct Connections
  data::Vector{Int}
  ptrs::Vector{Int}
end

struct UnstructuredGrid
  point_to_coords::Array{Float64,2}
  cell_to_points::Connections
  cell_to_ctype::Vector{Int}
  ctype_to_grid::Vector{UnstructuredGrid}
end

include("ReferenceGrids.jl")

end # module UnstructuredGrids
