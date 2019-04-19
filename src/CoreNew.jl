module CoreNew

using UnstructuredGrids.Helpers
using UnstructuredGrids.Kernels

import Base: ndims
import Base: ==
import Base: show

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

abstract type AbstractGrid end

celldata(::AbstractGrid)::AbstractCellData = @abstractmethod

refcells(::AbstractGrid)::AbstractVector{<:AbstractRefCell} = @abstractmethod

pointdata(::AbstractGrid)::AbstractPointData = @abstractmethod

coordinates(grid::AbstractGrid) = coordinates(pointdata(grid))

connections(grid::AbstractGrid) = connections(celldata(grid))

celltypes(grid::AbstractGrid) = celltypes(celldata(grid))

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

celldata(r::RefCell,dim::Integer) = r.cdata

reffaces(r::RefCell,dim::Integer) = r.rfaces

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

struct GridGraph{
  I<:Integer,
  C<:AbstractVector{<:AbstractConnections},
  D<:AbstractVector{<:AbstractConnections}} <:AbstractGridGraph
  ndims::I
  conn::C
  dualconn::D
end

ndims(g::GridGraph) = g.ndims

connections(g::GridGraph,dim::Integer) = g.conn[dim]

dualconnections(g::GridGraph,dim::Integer) = g.dualconn[dim]

# Definition of vertex

const VERTEX = RefCell(
  ndims = 0, faces = fill([Int[]],0), vtkid = 1, vtknodes = [1] )

end # module Core

