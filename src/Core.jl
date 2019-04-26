module Core

using UnstructuredGrids.Helpers
using UnstructuredGrids.Kernels

import Base: ndims
import Base: show

export RefCell
export Grid
export Connections
export Mesh
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
export append
export celldata
export reffaces
export pointdata

struct Connections{L<:AbstractVector{<:Integer},P<:AbstractVector{<:Integer}}
  list::L
  ptrs::P
end

list(c::Connections) = c.list

ptrs(c::Connections) = c.ptrs

function Connections(c::Vector{Vector{Int}})
  list, ptrs = generate_data_and_ptrs(c)
  Connections(list,ptrs)
end

function append(a::Connections{A,B},b::Connections{A,B}) where {A,B}
  p = append_ptrs(ptrs(a),ptrs(b))
  l = vcat(list(a),list(b))
  Connections(l,p)
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

struct CellData{
  L<:AbstractVector{<:Integer},
  P<:AbstractVector{<:Integer},
  C<:AbstractVector{<:Integer}}
  connections::Connections{L,P}
  celltypes::C
end

connections(d::CellData) = d.connections

celltypes(d::CellData) = d.celltypes

function CellData(
  connections::AbstractVector{<:AbstractVector{<:Integer}},
  celltypes::AbstractVector{<:Integer})
  c = Connections(connections)
  CellData(c,celltypes)
end

function append(a::CellData{A,B,C},b::CellData{A,B,C}) where {A,B,C}
  c = append(connections(a),connections(b))
  t = vcat(celltypes(a),celltypes(b))
  CellData(c,t)
end

struct PointData{P<:AbstractArray{<:Number,2}}
  coords::P
end

coordinates(p::PointData) = p.coords

ndims(p::PointData) = size(coordinates(p),1)

struct VtkData
  vtkid::Int
  vtknodes::Vector{Int}
end

vtkid(v::VtkData) = v.vtkid

vtknodes(v::VtkData) = v.vtknodes

struct RefCell
  cdata::Vector{CellData{Vector{Int},Vector{Int},Vector{Int}}}
  rfaces::Vector{RefCell}
  ndims::Int
  pdata::PointData{Array{Float64,2}}
  vtkdata::VtkData
end

celldata(r::RefCell,dim::Integer) = r.cdata[dim+1]

reffaces(r::RefCell) = r.rfaces

ndims(r::RefCell)::Integer = r.ndims

pointdata(r::RefCell) = r.pdata

vtkdata(r::RefCell) = r.vtkdata

coordinates(refcell::RefCell) = coordinates(pointdata(refcell))

connections(refcell::RefCell,dim::Integer) = connections(celldata(refcell,dim))

celltypes(refcell::RefCell,dim::Integer) = celltypes(celldata(refcell,dim))

vtkid(r::RefCell) = vtkid(vtkdata(r))

vtknodes(r::RefCell) = vtknodes(vtkdata(r))

function RefCell(;
  ndims::Int,
  faces::Vector{Vector{Vector{Int}}},
  reffaces::Vector{RefCell} = Vector{RefCell}(undef,0),
  facetypes::Vector{Vector{Int}} = fill(Int[],ndims),
  points::Array{Float64,2} = zeros(ndims,0),
  vtkid::Int = UNSET,
  vtknodes::Vector{Int} = Int[])
  @assert ndims == length(faces)
  cdata = Vector{CellData}(undef,ndims)
  for dim in 1:ndims
    _faces = faces[dim]
    _facetypes = facetypes[dim]
    cdata[dim] = CellData(_faces,_facetypes)
  end
  rfaces = reffaces
  pdata = PointData(points)
  vtkdata = VtkData(vtkid,vtknodes)
  RefCell(cdata,rfaces,ndims,pdata,vtkdata)
end

struct Grid{ C<:CellData, P<:PointData}
  cdata::C
  rcells::Vector{RefCell}
  pdata::P
end

celldata(g::Grid) = g.cdata

refcells(g::Grid) = g.rcells

pointdata(g::Grid) = g.pdata

coordinates(grid::Grid) = coordinates(pointdata(grid))

connections(grid::Grid) = connections(celldata(grid))

celltypes(grid::Grid) = celltypes(celldata(grid))

function ndims(grid::Grid)
  rcs = refcells(grid)
  rc = rcs[1]
  d = ndims(rc)
  for rc2 in rcs
    @assert d == ndims(rc2)
  end
  d
end

function Grid(
  points::AbstractArray{<:Number,2},
  celllist::AbstractVector{<:Integer},
  cellptrs::AbstractVector{<:Integer},
  celltypes::AbstractVector{<:Integer},
  refcells::Vector{RefCell})
  c = Connections(celllist,cellptrs)
  cdata = CellData(c,celltypes)
  pdata = PointData(points)
  Grid(cdata,refcells,pdata)
end

function Grid(
  points::AbstractArray{<:Number,2},
  cells::Connections,
  celltypes::AbstractVector{<:Integer},
  refcells::Vector{RefCell})
  cdata = CellData(cells,celltypes)
  pdata = PointData(points)
  Grid(cdata,refcells,pdata)
end

function Grid(
  grid::Grid,
  vertex_to_cells=_generate_face_to_cells(connections(grid));
  dim::Integer)
  cell_to_vertices = connections(grid)
  cell_to_ctype = celltypes(grid)
  ctype_to_refcell = refcells(grid)
  cell_to_faces = _generate_cell_to_faces(
      dim, cell_to_vertices, vertex_to_cells, cell_to_ctype, ctype_to_refcell)
  ftype_to_refface, ctype_to_lface_to_ftype = _prepare_ftypes(dim,ctype_to_refcell)
  face_to_ftype = _generate_face_to_ftype(
    cell_to_faces, cell_to_ctype, ctype_to_lface_to_ftype)
  face_to_vertices = _generate_face_to_vertices(
    dim, cell_to_vertices, cell_to_faces, cell_to_ctype, ctype_to_refcell)
  point_to_coords = coordinates(grid)
  Grid(point_to_coords, face_to_vertices, face_to_ftype, ftype_to_refface)
end

function Grid(r::RefCell;dim::Integer)
  cdata = celldata(r,dim)
  rcells = reffaces(r)
  pdata = pointdata(r)
  Grid(cdata,rcells,pdata)
end

struct GridGraph{ C<:Vector{<:Connections}, D<:Vector{<:Connections}}
  ndims::Int
  conn::C
  dualconn::D
end

ndims(g::GridGraph) = g.ndims

connections(g::GridGraph,dim::Integer) = g.conn[dim+1]

dualconnections(g::GridGraph,dim::Integer) = g.dualconn[dim+1]

function connections(graph::GridGraph;from::Integer,to::Integer)
  if ndims(graph) == from 
    return connections(graph,to)
  elseif ndims(graph) == to
    return dualconnections(graph,from)
  else
    @unreachable
  end
end

function GridGraph(grid::Grid)
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
      d, cell_to_vertices, vertex_to_cells, cell_to_ctype, ctype_to_refcell)
    dualconn[d+1] = _generate_face_to_cells(conn[d+1])
  end
  GridGraph(dim,conn,dualconn)
end

struct Mesh{C<:Vector{<:CellData},D<:Vector{<:Connections},P<:PointData}
  cdata::C
  dualc::D
  rcells::Vector{RefCell}
  ndims::Int
  pdata::P
end

dualconnections(r::Mesh,dim::Integer) = r.dualc[dim+1]

connections(r::Mesh,dim::Integer) = connections(celldata(r,dim))

celltypes(r::Mesh,dim::Integer) = celltypes(celldata(r,dim))

celldata(r::Mesh,dim::Integer) = r.cdata[dim+1]

refcells(r::Mesh) = r.rcells

ndims(r::Mesh)::Integer = r.ndims

pointdata(r::Mesh) = r.pdata

function Mesh(grid::Grid)
  dim = ndims(grid)
  cell_to_vertices = connections(grid)
  vertex_to_cells = _generate_face_to_cells(cell_to_vertices)
  cdata = CellData[]
  rcells = RefCell[]
  dualc = Connections[]
  for d in 0:(dim-1)
    fgrid = Grid(grid,vertex_to_cells,dim=d)
    c = connections(fgrid)
    t = celltypes(fgrid) .+ length(rcells)
    dc = _generate_face_to_cells(c)
    push!(cdata,CellData(c,t))
    push!(dualc,dc)
    rcells = vcat(rcells, refcells(fgrid))
  end
  c = connections(grid)
  t = celltypes(grid) .+ length(rcells)
  push!(cdata,CellData(c,t))
  push!(dualc,vertex_to_cells)
  rcells = vcat(rcells, refcells(grid))
  pdata = PointData(coordinates(grid))
  Mesh(cdata,dualc,rcells,dim,pdata)
end

function Grid(r::Mesh;dim::Integer)
  cdata = celldata(r,dim)
  rcells = refcells(r)
  pdata = pointdata(r)
  Grid(cdata,rcells,pdata)
end

#function GridGraph(m::Mesh;dim::Integer)
#  conn = Vector{Connections}(undef,dim)
#  dualconn = Vector{Connections}(undef,dim)
#  cell_to_vertices = connections(m,dim)
#  vertex_to_cells = dualconnections(m,dim)
#  conn[0+1] = cell_to_vertices
#  dualconn[0+1] = vertex_to_cells
#  cell_to_ctype = celltypes(m,dim)
#  ctype_to_refcell = refcells(m)
#  for d in 1:(dim-1)
#    vertex_to_faces = dualconnections(m,d)
#    conn[d+1] = _generate_cell_to_faces_from_faces(
#      dim, cell_to_vertices, vertex_to_faces, cell_to_ctype, ctype_to_refcell)
#    dualconn[d+1] = _generate_face_to_cells(conn[d+1])
#  end
#  GridGraph(dim,conn,dualconn)
#end

struct Model{M<:Mesh}
  mesh::M
  dim_to_cell_to_entity::Vector{Vector{Int}}
  entity_to_tags::Vector{Int}
  tags_to_names::Vector{String}
end

function Model(grid::Grid)
  mesh = Mesh(grid)
end

# Definition of vertex

const VERTEX = RefCell(
  ndims = 0, faces = fill([Int[]],0), vtkid = 1, vtknodes = [1] )

# Helpers

function _generate_face_to_cells(c::Connections)
  _data, _ptrs = generate_face_to_cells(list(c),ptrs(c))
  Connections(_data,_ptrs)
end

function _generate_cell_to_faces(
  dim::Integer,
  cell_to_vertices::Connections,
  vertex_to_cells::Connections,
  cell_to_ctype::Vector{<:Integer},
  ctype_to_refcell::Vector{RefCell})

  _generate_cell_to_faces(
    cell_to_vertices,
    vertex_to_cells,
    cell_to_ctype,
    [ connections(refcell,dim) for refcell in ctype_to_refcell])

end

function _generate_cell_to_faces(
  cell_to_vertices::Connections,
  vertex_to_cells::Connections,
  cell_to_ctype::AbstractVector{<:Integer},
  ctype_to_lface_to_lvertices::AbstractVector{<:Connections})

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

function _generate_cell_to_faces_from_faces(
  dim::Integer,
  cell_to_vertices::Connections,
  vertex_to_faces::Connections,
  cell_to_ctype::Vector{<:Integer},
  ctype_to_refcell::Vector{RefCell})

  _generate_cell_to_faces_from_faces(
    cell_to_vertices,
    vertex_to_faces,
    cell_to_ctype,
    [ connections(refcell,dim) for refcell in ctype_to_refcell])

end

function _generate_cell_to_faces_from_faces(
  cell_to_vertices::Connections,
  vertex_to_faces::Connections,
  cell_to_ctype::AbstractVector{<:Integer},
  ctype_to_lface_to_lvertices::AbstractVector{<:Connections})

  l, p = generate_cell_to_faces_from_faces(
    list(cell_to_vertices),
    ptrs(cell_to_vertices),
    [ list(c) for c in ctype_to_lface_to_lvertices ],
    [ ptrs(c) for c in ctype_to_lface_to_lvertices ],
    cell_to_ctype,
    list(vertex_to_faces),
    ptrs(vertex_to_faces))

  Connections(l,p)

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
    lftype_to_refface = reffaces(refcell)
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
    lface_to_lftype = celltypes(refcell,dim)
    lftype_to_ftype = ctype_to_lftype_to_ftype[ctype]
    lface_to_ftype = lftype_to_ftype[lface_to_lftype]
    ctype_to_lface_to_ftype[ctype] = lface_to_ftype
  end

  (ftype_to_refface, ctype_to_lface_to_ftype)

end

end # module Core

