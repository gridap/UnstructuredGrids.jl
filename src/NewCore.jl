module NewCore

using UnstructuredGrids.Helpers
using UnstructuredGrids.Kernels

import Base: ndims
import Base: show

export Connections
export RefCell
export Grid

struct Connections{L<:AbstractVector{<:Integer},P<:AbstractVector{<:Integer}}
  list::L
  ptrs::P
end

get_data(c::Connections) = c.list

get_ptrs(c::Connections) = c.ptrs

function Connections(c::AbstractVector{<:AbstractVector{<:Integer}})
  list, ptrs = generate_data_and_ptrs(c)
  Connections(list,ptrs)
end

function show(io::IO,c::Connections)
  clist = list(c)
  cptrs = ptrs(c)
  ncells = length(cptrs)-1
  for cell in 1:ncells
    a = cptrs[cell]
    b = cptrs[cell+1]-1
    println(io,"$cell -> $(clist[a:b])")
  end
end

struct RefCell
  ndims::Int
  dim_to_face_to_vertices::Connections{Vector{Int},Vector{Int}}
  dim_to_face_to_ftype::Vector{Int}
  dim_to_ftype_to_refface::Vector{Vector{RefCell}}
  vertex_to_coords::Array{Float64,2}
  vtkid::Int
  vtknodes::Vector{Int}
end

ndims(r::RefCell)::Integer = r.ndims

get_face_to_vertices(r::RefCell, dim::Integer) =r.dim_to_face_to_vertices[dim+1]

get_face_to_ftype(r::RefCell, dim::Integer) = r.dim_to_face_to_ftype[dim+1]

get_ftype_to_refface(r::RefCell, dim::Integer) =r.dim_to_ftype_to_refface[dim+1]

get_vertex_to_coords(r::RefCell) = r.vertex_to_coords

get_vtkid(r::RefCell) = r.vtkid

get_vtknodes(r::RefCell) = r.vtknodes

function RefCell(;
  ndims::Int,
  dim_to_face_to_vertices::Vector{Vector{Vector{Int}}},
  dim_to_face_to_ftype::Vector{Vector{Int}} = fill(Int[],ndims),
  dim_to_ftype_to_refface::Vector{Vector{RefCell}} = fill(RefCell[],ndims),
  vertex_to_coords::Array{Float64,2} = zeros(ndims,0),
  vtkid::Int = UNSET,
  vtknodes::Vector{Int} = Int[])

  @assert ndims == length(dim_to_face_to_vertices)
  @assert ndims == length(dim_to_face_to_ftype)
  @assert ndims == length(dim_to_ftype_to_refface)
  @assert ndims == size(vertex_to_coords,1)

  RefCell(
    ndims,
    dim_to_face_to_vertices,
    dim_to_face_to_ftype,
    dim_to_ftype_to_refface,
    vertex_to_coords,
    vtkid,
    vtknodes)

end

const VERTEX = RefCell(
  ndims = 0,
  dim_to_face_to_vertices = fill([Int[]],0),
  vtkid = 1,
  vtknodes = [1] )

struct Grid{
  C<:Connections,
  T<:AbstractVector{<:Integer},
  X<:AbstractArray{<:Number,2}}

  cell_to_vertices::C
  cell_to_ctype::T
  ctype_to_refcell::Vector{RefCell}
  vertex_to_coords::X

end

get_cell_to_vertices(g::Grid) = g.cell_to_vertices

get_cell_to_ctype(g::Grid) = g.cell_to_ctype   

get_ctype_to_refcell(g::Grid) = g.ctype_to_refcell

get_vertex_to_coords(g::Grid) = g.vertex_to_coords

function get_ctype_to_lface_to_lvertices(g::Grid,dim::Integer)
  get_ctype_to_lface_to_lvertices(get_ctype_to_refcell(g),dim)
end

function get_ctype_to_lface_to_lvertices(
  ctype_to_refcell::Vector{RefCell}, dim::Integer)
  [ get_face_to_vertices(refcell,dim) for refcell in ctype_to_refcell ]
end

function generate_face_to_cells(cell_to_faces::Connections)
  cell_to_faces_data = get_data(cell_to_faces)
  cell_to_faces_ptrs = get_ptrs(cell_to_faces)
  face_to_cells_data, face_to_cells_ptrs = generate_face_to_cells(
  cell_to_faces_data, cell_to_faces_ptrs)
  Connections(face_to_cells_data, face_to_cells_ptrs)
end

function generate_cell_to_faces(
  cell_to_vertices::Connections,
  vertex_to_cells::Connections,
  cell_to_ctype::AbstractVector{<:Integer},
  ctype_to_refcell::Vector{RefCell},
  dim::Integer)

  ctype_to_lface_to_lvertices = get_ctype_to_lface_to_lvertices(
    ctype_to_refcell, dim)

  generate_cell_to_faces(
    cell_to_vertices,
    vertex_to_cells,
    cell_to_ctype,
    ctype_to_lface_to_lvertices)

end

function generate_cell_to_faces(
  cell_to_vertices::Connections,
  vertex_to_cells::Connections,
  cell_to_ctype::AbstractVector{<:Integer},
  ctype_to_lface_to_lvertices::Vector{<:Connections})

  cell_to_vertices_data = get_data(cell_to_vertices)
  cell_to_vertices_ptrs = get_ptrs(cell_to_vertices)

  vertex_to_cells_data = get_data(vertex_to_cells)
  vertex_to_cells_ptrs = get_ptrs(vertex_to_cells)

  ctype_to_lface_to_lvertices_data, ctype_to_lface_to_lvertices_ptrs = (
    _split_vector_of_indextoindices( ctype_to_lface_to_lvertices) )

  cell_to_faces_data, cell_to_faces_ptrs = generate_cell_to_faces(
    cell_to_vertices_data,
    cell_to_vertices_ptrs,
    ctype_to_lface_to_lvertices_data,
    ctype_to_lface_to_lvertices_ptrs,
    cell_to_ctype,
    vertex_to_cells_data,
    vertex_to_cells_ptrs)

  Connections(cell_to_faces_data, cell_to_faces_ptrs)

end

function generate_cell_to_faces_from_faces(
  cell_to_vertices::Connections,
  vertex_to_faces::Connections,
  cell_to_ctype::AbstractVector{<:Integer},
  ctype_to_refcell::Vector{RefCell},
  dim::Integer)

  ctype_to_lface_to_lvertices = get_ctype_to_lface_to_lvertices(
    ctype_to_refcell, dim)

  generate_cell_to_faces(
    cell_to_vertices,
    vertex_to_cells,
    cell_to_ctype,
    ctype_to_lface_to_lvertices)

end

function generate_cell_to_faces_from_faces(
  cell_to_vertices::Connections,
  vertex_to_faces::Connections,
  cell_to_ctype::AbstractVector{<:Integer},
  ctype_to_lface_to_lvertices::Vector{<:Connections})

  cell_to_vertices_data = get_data(cell_to_vertices)
  cell_to_vertices_ptrs = get_ptrs(cell_to_vertices)

  vertex_to_faces_data = get_data(vertex_to_faces)
  vertex_to_faces_ptrs = get_ptrs(vertex_to_faces)

  ctype_to_lface_to_lvertices_data, ctype_to_lface_to_lvertices_ptrs = (
    _split_vector_of_indextoindices( ctype_to_lface_to_lvertices) )

  cell_to_faces_data, cell_to_faces_ptrs = generate_cell_to_faces_from_faces(
    cell_to_vertices_data,
    cell_to_vertices_ptrs,
    ctype_to_lface_to_lvertices_data,
    ctype_to_lface_to_lvertices_ptrs,
    cell_to_ctype,
    vertex_to_faces_data,
    vertex_to_faces_ptrs)

  Connections(cell_to_faces_data, cell_to_faces_ptrs)

end

# Helpers

function _split_vector_of_indextoindices(
  ctype_to_lface_to_lvertices::Vector{Connections})

  ctype_to_lface_to_lvertices_data = [
    get_data(lface_to_lvertices)
    for lface_to_lvertices in ctype_to_lface_to_lvertices ]

  ctype_to_lface_to_lvertices_ptrs = [
    get_ptrs(lface_to_lvertices)
    for lface_to_lvertices in ctype_to_lface_to_lvertices ]

  (ctype_to_lface_to_lvertices_data, ctype_to_lface_to_lvertices_ptrs)
end

end # module NewCore
