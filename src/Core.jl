module Core

using UnstructuredGrids.Helpers
using UnstructuredGrids.Kernels

export UGrid
export RefCell
export Grid
export GridGraph
export Connections
export VERTEX
export gridpoints, gridcells, celltypes, refcells
export connections
export vtkid, vtknodes
import Base: ==
import Base: show

# Structs

struct Connections
  data::Vector{Int}
  ptrs::Vector{Int}
end

struct GridData
  points::Array{Float64,2}
  cells::Connections
  celltypes::Vector{Int}
end

struct RefCell
  data::GridData
  refcells::Vector{RefCell}
  vtkid::Int
  vtknodes::Vector{Int}
end

struct UGrid
  data::GridData
  refcells::Vector{RefCell}
end

const Grid = Union{UGrid,RefCell}

struct GridGraph
  maxdim::Int
  cell_to_ctype::Vector{Int}
  dim_to_refcells::Dict{Int,Vector{RefCell}}
  dim_to_cell_to_faces::Dict{Int,Connections}
  dim_to_face_to_cells::Dict{Int,Connections}
end

# Constructors

function RefCell(;
  points::Array{Float64,2},
  cells::Vector{Vector{Int}},
  celltypes::Vector{Int},
  refcells::Vector{RefCell},
  vtkid::Int,
  vtknodes::Vector{Int})
  data = GridData(points,cells,celltypes)
  RefCell( data, refcells, vtkid, vtknodes)
end

function UGrid(
  points::Array{Float64,2},
  cellsdata::Vector{Int},
  cellsptrs::Vector{Int},
  celltypes::Vector{Int},
  refcells::Vector{RefCell})
  cells = Connections(cellsdata,cellsptrs)
  UGrid(points,cells,celltypes,refcells)
end

function UGrid(
  points::Array{Float64,2},
  cells::Connections,
  celltypes::Vector{Int},
  refcells::Vector{RefCell})
  data = GridData(points,cells,celltypes)
  UGrid(data,refcells)
end

function UGrid(grid::Grid;dim::Int,graph::GridGraph=GridGraph(grid))
  point_to_coords = gridpoints(grid)
  cell_to_ctype = celltypes(grid)
  cell_to_vertices = connections(grid)
  maxdim = graph.maxdim
  cell_to_faces = connections(graph,from=maxdim,to=dim)
  ctype_to_refcell = graph.dim_to_refcells[dim]
  ftype_to_refface, ctype_to_lface_to_ftype = _prepare_ftypes(ctype_to_refcell)
  face_to_ftype = _generate_face_to_ftype(
    cell_to_faces, cell_to_ctype, ctype_to_lface_to_ftype)
  face_to_vertices = _generate_face_to_vertices(
    cell_to_vertices, cell_to_faces, cell_to_ctype, ctype_to_refcell)
  UGrid(point_to_coords, face_to_vertices, face_to_ftype, ftype_to_refface)
end

function Connections(c::Vector{Vector{Int}})
  data, ptrs = generate_data_and_ptrs(c)
  Connections(data,ptrs)
end

function GridData(points,cells::Vector{Vector{Int}},celltypes)
  co = Connections(cells)
  GridData(points,co,celltypes)
end

function GridGraph(grid::Grid)
  cell_to_vertices = connections(grid)
  vertex_to_cells = _generate_face_to_cells(cell_to_vertices)
  cell_to_ctype = celltypes(grid)
  dim_to_refcells = Dict{Int,Vector{RefCell}}()
  # TODO Don't like to extract this from points since
  # in some situations points will not be set
  # WARNING!! this does not work when the dimension of the
  # cells is different from the points
  d = size(gridpoints(grid),1)
  dim_to_refcells[d-1] = refcells(grid)
  _gridgraph(d,cell_to_vertices,vertex_to_cells,cell_to_ctype,dim_to_refcells)
end

# Definition of vertex

const VERTEX = RefCell(
  points = zeros(0,1),
  cells = [Int[]],
  celltypes = Int[],
  refcells = Vector{RefCell}(undef,0),
  vtkid = 1,
  vtknodes = [1])

# Behavior

function show(io::IO,c::Connections)
  ncells = length(c.ptrs)-1
  for cell in 1:ncells
    a = c.ptrs[cell]
    b = c.ptrs[cell+1]-1
    println(io,"$cell -> $(c.data[a:b])")
  end
end

function (==)(a::Connections,b::Connections)
  a.data == b.data && a.ptrs == b.ptrs
end

function (==)(a::GridData,b::GridData)
  !(a.points ≈ b.points) && return false
  a.cells != b.cells && return false
  a.celltypes != b.celltypes && return false
  return true
end

function (==)(a::RefCell,b::RefCell)
  a.data == b.data && a.refcells == b.refcells
end

function (==)(a::UGrid,b::UGrid)
  a.data == b.data && a.refcells == b.refcells
end

gridpoints(grid::Grid) = grid.data.points

gridcells(grid::Grid) = (grid.data.cells.data, grid.data.cells.ptrs)

connections(grid::Grid) = grid.data.cells

celltypes(grid::Grid) = grid.data.celltypes

refcells(grid::Grid) = grid.refcells

vtkid(grid::RefCell) = grid.vtkid

vtknodes(grid::RefCell) = grid.vtknodes

function connections(graph::GridGraph;from::Int,to::Int)
  if graph.maxdim == from 
    return graph.dim_to_cell_to_faces[to]
  elseif graph.maxdim == to
    return graph.dim_to_face_to_cells[from]
  else
    @unreachable
  end
end

# Helpers

function _gridgraph(
  maxdim,cell_to_vertices,vertex_to_cells,cell_to_ctype,dim_to_refcells)
  dim_to_cell_to_faces = Dict{Int,Connections}()
  dim_to_face_to_cells = Dict{Int,Connections}()
  dim_to_cell_to_faces[0] = cell_to_vertices
  dim_to_face_to_cells[0] = vertex_to_cells
  for (d,refcells) in dim_to_refcells
    dim_to_cell_to_faces[d] = _generate_cell_to_faces(
      cell_to_vertices, vertex_to_cells, cell_to_ctype, refcells)
    dim_to_face_to_cells[d] = _generate_face_to_cells(dim_to_cell_to_faces[d])
  end
  GridGraph(
    maxdim,
    cell_to_ctype,
    dim_to_refcells,
    dim_to_cell_to_faces,
    dim_to_face_to_cells)
end

function _generate_face_to_cells(c::Connections)
  data, ptrs = generate_face_to_cells(c.data,c.ptrs)
  Connections(data,ptrs)
end

function _generate_cell_to_faces(
  cell_to_vertices::Connections,
  vertex_to_cells::Connections,
  cell_to_ctype::Vector{Int},
  ctype_to_refcell::Vector{RefCell})

  _generate_cell_to_faces(
    cell_to_vertices,
    vertex_to_cells,
    cell_to_ctype,
    [ connections(refcell) for refcell in ctype_to_refcell])

end

function _generate_cell_to_faces(
  cell_to_vertices::Connections,
  vertex_to_cells::Connections,
  cell_to_ctype::Vector{Int},
  ctype_to_lface_to_lvertices::Vector{Connections})

  data, ptrs = generate_cell_to_faces(
    cell_to_vertices.data,
    cell_to_vertices.ptrs,
    [ c.data for c in ctype_to_lface_to_lvertices ],
    [ c.ptrs for c in ctype_to_lface_to_lvertices ],
    cell_to_ctype,
    vertex_to_cells.data,
    vertex_to_cells.ptrs)

  Connections(data,ptrs)

end

function _generate_face_to_vertices(
  cell_to_vertices::Connections,
  cell_to_faces::Connections,
  cell_to_ctype::Vector{Int},
  ctype_to_refcell::Vector{RefCell})
  _generate_face_to_vertices(
    cell_to_vertices,
    cell_to_faces,
    cell_to_ctype,
    [ connections(refcell) for refcell in ctype_to_refcell])
end

function _generate_face_to_vertices(
  cell_to_vertices::Connections,
  cell_to_faces::Connections,
  cell_to_ctype::Vector{Int},
  ctype_to_lface_to_lvertices::Vector{Connections})
  data, ptrs = generate_face_to_vertices(
    cell_to_vertices.data,
    cell_to_vertices.ptrs,
    cell_to_faces.data,
    cell_to_faces.ptrs,
    cell_to_ctype,
    [ c.data for c in ctype_to_lface_to_lvertices ],
    [ c.ptrs for c in ctype_to_lface_to_lvertices ])
  Connections(data,ptrs)
end

function _generate_face_to_ftype(
  cell_to_faces::Connections,
  cell_to_ctype::Vector{Int},
  ctype_to_lface_to_ftype::Vector{Vector{Int}})
  generate_face_to_ftype(
    cell_to_faces.data,
    cell_to_faces.ptrs,
    cell_to_ctype,
    ctype_to_lface_to_ftype)
end

function _prepare_ftypes(ctype_to_refcell)

  i_to_refface = Vector{RefCell}(undef,0)
  nctypes = length(ctype_to_refcell)
  ctype_to_lftype_to_i = Vector{Vector{Int}}(undef,nctypes)

  i = 1
  for (ctype,refcell) in enumerate(ctype_to_refcell)
    lftype_to_refface = refcells(refcell)
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
    lface_to_lftype = celltypes(refcell)
    lftype_to_ftype = ctype_to_lftype_to_ftype[ctype]
    lface_to_ftype = lftype_to_ftype[lface_to_lftype]
    ctype_to_lface_to_ftype[ctype] = lface_to_ftype
  end

  (ftype_to_refface, ctype_to_lface_to_ftype)

end

end #module Core
