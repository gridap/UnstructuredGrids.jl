module CoreNew

using UnstructuredGrids.Helpers
using UnstructuredGrids.Kernels

import Base: ndims
import Base: ==
import Base: show

export AbstractRefCell
export AbstractGrid
export AbstractConnections
export AbstractGridGraph
export RefCell
export Grid
export Connections
export GridGraph
export VERTEX
export coordinates
export connections
export celltypes
export list
export ptrs
export vtkid
export vtknodes
export refcells

# Interfaces

abstract type AbstractConnections end

list(::AbstractConnections)::AbstractVector{<:Integer} = @abstractmethod

ptrs(::AbstractConnections)::AbstractVector{<:Integer} = @abstractmethod

function show(io::IO,c::AbstractConnections)
  clist = list(c)
  cptrs = ptrs(c)
  ncells = length(cptrs)-1
  for cell in 1:ncells
    a = cptrs[cell]
    b = cptrs[cell+1]-1
    println(io,"$cell -> $(clist[a:b])")
  end
end

abstract type AbstractCellData end

connections(::AbstractCellData)::AbstractConnections = @abstractmethod

celltypes(::AbstractCellData)::AbstractVector{<:Integer} = @abstractmethod

abstract type AbstractPointData end

coordinates(::AbstractPointData)::AbstractArray{<:Number,2} = @abstractmethod

ndims(p::AbstractPointData) = size(coordinates(p),1)

abstract type AbstractVtkData end

vtkid(::AbstractVtkData)::Integer = @abstractmethod

vtknodes(::AbstractVtkData)::AbstractVector{<:Integer} = @abstractmethod

abstract type AbstractRefCell end

celldata(::AbstractRefCell,dim::Integer)::AbstractCellData = @abstractmethod

reffaces(::AbstractRefCell,dim::Integer)::AbstractVector{<:AbstractRefCell} = @abstractmethod

ndims(::AbstractRefCell)::Integer = @abstractmethod

pointdata(::AbstractRefCell)::AbstractPointData = @abstractmethod

vtkdata(::AbstractRefCell)::AbstractVtkData = @abstractmethod

coordinates(refcell::AbstractRefCell) = coordinates(pointdata(refcell))

connections(refcell::AbstractRefCell,dim::Integer) = connections(celldata(refcell,dim))

celltypes(refcell::AbstractRefCell,dim::Integer) = celltypes(celldata(refcell,dim))

vtkid(r::AbstractRefCell) = vtkid(vtkdata(r))

vtknodes(r::AbstractRefCell) = vtknodes(vtkdata(r))

abstract type AbstractGrid end

celldata(::AbstractGrid)::AbstractCellData = @abstractmethod

refcells(::AbstractGrid)::AbstractVector{<:AbstractRefCell} = @abstractmethod

pointdata(::AbstractGrid)::AbstractPointData = @abstractmethod

coordinates(grid::AbstractGrid) = coordinates(pointdata(grid))

connections(grid::AbstractGrid) = connections(celldata(grid))

celltypes(grid::AbstractGrid) = celltypes(celldata(grid))

function ndims(grid::AbstractGrid)
  rcs = refcells(grid)
  rc = rcs[1]
  d = ndims(rc)
  for rc2 in rcs
    @assert d == ndims(rc2)
  end
  d
end

abstract type AbstractGridGraph end

ndims(::AbstractGridGraph)::Integer = @abstractmethod

connections(::AbstractGridGraph,dim::Integer)::AbstractConnections = @abstractmethod

dualconnections(::AbstractGridGraph,dim::Integer)::AbstractConnections = @abstractmethod

function connections(graph::AbstractGridGraph;from::Integer,to::Integer)
  if ndims(graph) == from 
    return connections(graph,to)
  elseif ndims(graph) == to
    return dualconnections(graph,from)
  else
    @unreachable
  end
end

# Concrete implementations

struct Connections{
  L<:AbstractVector{<:Integer},
  P<:AbstractVector{<:Integer}} <: AbstractConnections
  list::L
  ptrs::P
end

list(c::Connections) = c.list

ptrs(c::Connections) = c.ptrs

function Connections(c::AbstractVector{<:AbstractVector{<:Integer}})
  list, ptrs = generate_data_and_ptrs(c)
  Connections(list,ptrs)
end

struct CellData{
  A<:AbstractConnections,
  B<:AbstractVector{<:Integer}} <: AbstractCellData
  connections::A
  celltypes::B
end

connections(d::CellData) = d.connections

celltypes(d::CellData) = d.celltypes

function CellData(
  connections::AbstractVector{<:AbstractVector{<:Integer}},
  celltypes::AbstractVector{<:Integer})
  c = Connections(connections)
  CellData(c,celltypes)
end

struct PointData{P<:AbstractArray{<:Number,2}} <: AbstractPointData
  coords::P
end

coordinates(p::PointData) = p.coords

struct VtkData{I<:Integer,V<:AbstractVector{<:Integer}} <: AbstractVtkData
  vtkid::I
  vtknodes::V
end

vtkid(v::VtkData) = v.vtkid

vtknodes(v::VtkData) = v.vtknodes

struct RefCell{
  C<:AbstractVector{<:AbstractCellData},
  R<:AbstractVector{<:AbstractRefCell},
  I<:Integer,
  P<:AbstractPointData,
  V<:AbstractVtkData} <: AbstractRefCell
  cdata::C
  rfaces::R
  ndims::I
  pdata::P
  vtkdata::V
end

celldata(r::RefCell,dim::Integer) = r.cdata[dim]

reffaces(r::RefCell,dim::Integer) = r.rfaces[dim]

ndims(r::RefCell)::Integer = r.ndims

pointdata(r::RefCell) = r.pdata

vtkdata(r::RefCell) = r.vtkdata

function RefCell(;
  ndims::Integer,
  faces::AbstractVector{<:AbstractVector{<:AbstractVector{<:Integer}}},
  reffaces::AbstractVector{<:AbstractRefCell} = Vector{RefCell}(undef,0),
  facetypes::AbstractVector{<:AbstractVector{<:Integer}} = fill(Int[],ndims),
  points::AbstractArray{<:Number,2} = zeros(ndims,0),
  vtkid::Integer = UNSET,
  vtknodes::AbstractVector{<:Integer} = Int[])
  @assert ndims == length(faces)
  cdata = Vector{CellData}(undef,ndims)
  for dim in 1:ndims
    _faces = faces[dim]
    _facetypes = facetypes[dim]
    cdata[dim] = CellData(_faces,_facetypes)
  end
  rfaces = reffaces
  ndims
  pdata = PointData(points)
  vtkdata = VtkData(vtkid,vtknodes)
  RefCell(cdata,rfaces,ndims,pdata,vtkdata)
end

struct Grid{
  C<:AbstractCellData,
  R<:AbstractVector{<:AbstractRefCell},
  P<:AbstractPointData} <: AbstractGrid
  cdata::C
  rcells::R
  pdata::P
end

celldata(g::Grid) = g.cdata

refcells(g::Grid) = g.rcells

pointdata(g::Grid) = g.pdata

function Grid(
  points::AbstractArray{<:Number,2},
  celllist::AbstractVector{<:Integer},
  cellptrs::AbstractVector{<:Integer},
  celltypes::AbstractVector{<:Integer},
  refcells::AbstractVector{<:AbstractRefCell})
  c = Connections(celllist,cellptrs)
  cdata = CellData(c,celltypes)
  pdata = PointData(points)
  Grid(cdata,refcells,pdata)
end

struct GridGraph{
  I<:Integer,
  C<:AbstractVector{<:AbstractConnections},
  D<:AbstractVector{<:AbstractConnections}} <:AbstractGridGraph
  ndims::I
  conn::C
  dualconn::D
end

ndims(g::GridGraph) = g.ndims

connections(g::GridGraph,dim::Integer) = g.conn[dim+1]

dualconnections(g::GridGraph,dim::Integer) = g.dualconn[dim+1]

function GridGraph(grid::AbstractGrid)
  dim = ndims(grid)
  conn = Vector{Connections}(undef,dim)
  dualconn = Vector{Connections}(undef,dim)
  cell_to_vertices = connections(grid)
  vertex_to_cells = _generate_face_to_cells(cell_to_vertices)
  conn[0+1] = cell_to_vertices
  dualconn[0+1] = vertex_to_cells
  cell_to_ctype = celltypes(grid)
  ctype_to_refcell = refcells(grid)
  for d in 1:dim-1
    conn[d+1] = _generate_cell_to_faces(
      dim, cell_to_vertices, vertex_to_cells, cell_to_ctype, ctype_to_refcell)
    dualconn[d+1] = _generate_face_to_cells(conn[d+1])
  end
  GridGraph(dim,conn,dualconn)
end

# Definition of vertex

const VERTEX = RefCell(
  ndims = 0, faces = fill([Int[]],0), vtkid = 1, vtknodes = [1] )

# Helpers

function _generate_face_to_cells(c::AbstractConnections)
  _data, _ptrs = generate_face_to_cells(list(c),ptrs(c))
  Connections(_data,_ptrs)
end

function _generate_cell_to_faces(
  dim::Integer,
  cell_to_vertices::AbstractConnections,
  vertex_to_cells::AbstractConnections,
  cell_to_ctype::Vector{<:Integer},
  ctype_to_refcell::Vector{<:AbstractRefCell})

  _generate_cell_to_faces(
    cell_to_vertices,
    vertex_to_cells,
    cell_to_ctype,
    [ connections(refcell,dim) for refcell in ctype_to_refcell])

end

function _generate_cell_to_faces(
  cell_to_vertices::AbstractConnections,
  vertex_to_cells::AbstractConnections,
  cell_to_ctype::AbstractVector{<:Integer},
  ctype_to_lface_to_lvertices::AbstractVector{<:AbstractConnections})

  l, p = generate_cell_to_faces(
    list(cell_to_vertices),
    ptrs(cell_to_vertices),
    [ list(c) for c in ctype_to_lface_to_lvertices ],
    [ ptrs(c) for c in ctype_to_lface_to_lvertices ],
    cell_to_ctype,
    list(vertex_to_cells),
    ptrs(vertex_to_cells))

  Connections(l,p)

end


end # module Core

