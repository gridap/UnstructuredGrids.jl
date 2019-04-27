module Core

using UnstructuredGrids.Helpers
using UnstructuredGrids.Kernels

import Base: ndims
import Base: show

export Connections
export RefCell
export Grid
export VERTEX
export list
export ptrs
export coordinates
export connections
export reffaces
export facetypes
export vtkid
export vtknodes
export celltypes
export refcells


struct Connections{L<:AbstractVector{<:Integer},P<:AbstractVector{<:Integer}}
  list::L
  ptrs::P
end

list(c::Connections) = c.list

ptrs(c::Connections) = c.ptrs

function Connections(c::AbstractVector{<:AbstractVector{<:Integer}})
  list, ptrs = generate_data_and_ptrs(c)
  Connections(list,ptrs)
end

function split_vector_of_connections(
  cs::Vector{<:Connections})
  cdata = [ list(c) for c in cs ]
  cptrs = [ ptrs(c) for c in cs ]
  (cdata, cptrs)
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
  faces::Vector{Connections{Vector{Int},Vector{Int}}}
  facetypes::Vector{Vector{Int}}
  reffaces::Vector{Vector{RefCell}}
  coordinates::Array{Float64,2}
  vtkid::Int
  vtknodes::Vector{Int}
end

ndims(r::RefCell)::Integer = r.ndims

connections(r::RefCell, dim::Integer) =r.faces[dim+1]

facetypes(r::RefCell, dim::Integer) = r.facetypes[dim+1]

reffaces(r::RefCell, dim::Integer) =r.reffaces[dim+1]

coordinates(r::RefCell) = r.coordinates

vtkid(r::RefCell) = r.vtkid

vtknodes(r::RefCell) = r.vtknodes

function RefCell(;
  ndims::Int,
  faces::Vector{Vector{Vector{Int}}},
  facetypes::Vector{Vector{Int}} = fill(Int[],ndims),
  reffaces::Vector{Vector{RefCell}} = fill(RefCell[],ndims),
  coordinates::Array{Float64,2} = zeros(ndims,0),
  vtkid::Int = UNSET,
  vtknodes::Vector{Int} = Int[])

  @assert ndims == length(faces)
  @assert ndims == length(facetypes)
  @assert ndims == length(reffaces)
  @assert ndims == size(coordinates,1)

  _faces = [ Connections(f) for f in faces ]

  RefCell( ndims, _faces, facetypes, reffaces, coordinates, vtkid, vtknodes)

end

const VERTEX = RefCell(
  ndims = 0, faces = fill([Int[]],0), vtkid = 1, vtknodes = [1] )

struct Grid{
  C<:Connections,
  T<:AbstractVector{<:Integer},
  X<:AbstractArray{<:Number,2}}

  cells::C
  celltypes::T
  refcells::Vector{RefCell}
  coordinates::X

end

connections(g::Grid) = g.cells

celltypes(g::Grid) = g.celltypes

refcells(g::Grid) = g.refcells

coordinates(g::Grid) = g.coordinates

function refconnections(g::Grid,dim::Integer)
  refconnections(refcells(g),dim)
end

function refconnections(refcells::Vector{RefCell}, dim::Integer)
  [ connections(refcell,dim) for refcell in refcells ]
end

function Grid( cellsdata, cellsptrs, celltypes, refcells, coordinates)
  cells = Connections(cellsdata, cellsptrs)
  Grid( cells, celltypes, refcells, coordinates)
end

function Grid(r::RefCell;dim::Integer)
  cells = connections(r,dim)
  celltypes = facetypes(r,dim)
  refcells = reffaces(r,dim)
  coords = coordinates(r)
  Grid( cells, celltypes, refcells, coords)
end

function Grid(
  grid::Grid,
  vertex_to_cells=_generate_face_to_cells(connections(grid));
  dim::Integer)
  cell_to_vertices = connections(grid)
  cell_to_ctype = celltypes(grid)
  ctype_to_refcell = refcells(grid)
  cell_to_faces = _generate_cell_to_faces(
      cell_to_vertices, vertex_to_cells, cell_to_ctype, ctype_to_refcell,dim)
  ftype_to_refface, ctype_to_lface_to_ftype = _prepare_ftypes(dim,ctype_to_refcell)
  face_to_ftype = _generate_face_to_ftype(
    cell_to_faces, cell_to_ctype, ctype_to_lface_to_ftype)
  face_to_vertices = _generate_face_to_vertices(
    dim, cell_to_vertices, cell_to_faces, cell_to_ctype, ctype_to_refcell)
  point_to_coords = coordinates(grid)
  Grid(face_to_vertices, face_to_ftype, ftype_to_refface, point_to_coords)
end

function _generate_face_to_cells(cell_to_faces::Connections)
  cell_to_faces_data = list(cell_to_faces)
  cell_to_faces_ptrs = ptrs(cell_to_faces)
  face_to_cells_data, face_to_cells_ptrs = generate_face_to_cells(
  cell_to_faces_data, cell_to_faces_ptrs)
  Connections(face_to_cells_data, face_to_cells_ptrs)
end

function _generate_cell_to_faces(
  cell_to_vertices::Connections,
  vertex_to_cells::Connections,
  cell_to_ctype::AbstractVector{<:Integer},
  ctype_to_refcell::Vector{RefCell},
  dim::Integer)

  ctype_to_lface_to_lvertices = refconnections(ctype_to_refcell, dim)

  _generate_cell_to_faces(
    cell_to_vertices,
    vertex_to_cells,
    cell_to_ctype,
    ctype_to_lface_to_lvertices)

end

function _generate_cell_to_faces(
  cell_to_vertices::Connections,
  vertex_to_cells::Connections,
  cell_to_ctype::AbstractVector{<:Integer},
  ctype_to_lface_to_lvertices::Vector{<:Connections})

  cell_to_vertices_data = list(cell_to_vertices)
  cell_to_vertices_ptrs = ptrs(cell_to_vertices)

  vertex_to_cells_data = list(vertex_to_cells)
  vertex_to_cells_ptrs = ptrs(vertex_to_cells)

  ctype_to_lface_to_lvertices_data, ctype_to_lface_to_lvertices_ptrs = (
    split_vector_of_connections( ctype_to_lface_to_lvertices) )

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

function _generate_cell_to_faces_from_faces(
  cell_to_vertices::Connections,
  vertex_to_faces::Connections,
  cell_to_ctype::AbstractVector{<:Integer},
  ctype_to_refcell::Vector{RefCell},
  dim::Integer)

  ctype_to_lface_to_lvertices = refconnections(
    ctype_to_refcell, dim)

  _generate_cell_to_faces_from_faces(
    cell_to_vertices,
    vertex_to_faces,
    cell_to_ctype,
    ctype_to_lface_to_lvertices)

end

function _generate_cell_to_faces_from_faces(
  cell_to_vertices::Connections,
  vertex_to_faces::Connections,
  cell_to_ctype::AbstractVector{<:Integer},
  ctype_to_lface_to_lvertices::Vector{<:Connections})

  cell_to_vertices_data = list(cell_to_vertices)
  cell_to_vertices_ptrs = ptrs(cell_to_vertices)

  vertex_to_faces_data = list(vertex_to_faces)
  vertex_to_faces_ptrs = ptrs(vertex_to_faces)

  ctype_to_lface_to_lvertices_data, ctype_to_lface_to_lvertices_ptrs = (
    split_vector_of_connections( ctype_to_lface_to_lvertices) )

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

function _generate_face_to_vertices(
  dim::Integer,
  cell_to_vertices::Connections,
  cell_to_faces::Connections,
  cell_to_ctype::AbstractVector{<:Integer},
  ctype_to_refcell::Vector{RefCell})
  _generate_face_to_vertices(
    cell_to_vertices,
    cell_to_faces,
    cell_to_ctype,
    [ connections(refcell,dim) for refcell in ctype_to_refcell])
end

function _generate_face_to_vertices(
  cell_to_vertices::Connections,
  cell_to_faces::Connections,
  cell_to_ctype::AbstractVector{<:Integer},
  ctype_to_lface_to_lvertices::AbstractVector{<:Connections})
  l, p = generate_face_to_vertices(
    list(cell_to_vertices),
    ptrs(cell_to_vertices),
    list(cell_to_faces),
    ptrs(cell_to_faces),
    cell_to_ctype,
    [ list(c) for c in ctype_to_lface_to_lvertices ],
    [ ptrs(c) for c in ctype_to_lface_to_lvertices ])
  Connections(l,p)
end

function _generate_face_to_ftype(
  cell_to_faces::Connections,
  cell_to_ctype::AbstractVector{<:Integer},
  ctype_to_lface_to_ftype::AbstractVector{<:AbstractVector{<:Integer}})
  generate_face_to_ftype(
    list(cell_to_faces),
    ptrs(cell_to_faces),
    cell_to_ctype,
    ctype_to_lface_to_ftype)
end

function _prepare_ftypes(dim,ctype_to_refcell)

  i_to_refface = Vector{RefCell}(undef,0)
  nctypes = length(ctype_to_refcell)
  ctype_to_lftype_to_i = Vector{Vector{Int}}(undef,nctypes)

  i = 1
  for (ctype,refcell) in enumerate(ctype_to_refcell)
    lftype_to_refface = reffaces(refcell,dim)
    lftype_to_i = Vector{Int}(undef,length(lftype_to_refface))
    for (lftype,refface) in enumerate(lftype_to_refface)
      push!(i_to_refface,refface)
      lftype_to_i[lftype] = i
      i +=1
    end
    ctype_to_lftype_to_i[ctype] = lftype_to_i
  end

  ftype_to_refface = unique(i_to_refface)
  i_to_ftype = indexin(i_to_refface,ftype_to_refface)

  ctype_to_lftype_to_ftype = copy(ctype_to_lftype_to_i)
  for ctype in 1:length(ctype_to_lftype_to_i)
    for lftype in 1:length(ctype_to_lftype_to_i[ctype])
      i = ctype_to_lftype_to_i[ctype][lftype]
      ftype = i_to_ftype[i]
      ctype_to_lftype_to_ftype[ctype][lftype] = ftype
    end
  end

  ctype_to_lface_to_ftype = Vector{Vector{Int}}(undef,nctypes)
  for (ctype,refcell) in enumerate(ctype_to_refcell)
    lface_to_lftype = facetypes(refcell,dim)
    lftype_to_ftype = ctype_to_lftype_to_ftype[ctype]
    lface_to_ftype = lftype_to_ftype[lface_to_lftype]
    ctype_to_lface_to_ftype[ctype] = lface_to_ftype
  end

  (ftype_to_refface, ctype_to_lface_to_ftype)

end

end # module Core
