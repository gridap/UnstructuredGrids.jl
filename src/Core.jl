module Core

using UnstructuredGrids.Kernels: generate_data_and_ptrs

export UGrid
export RefCell
export Grid
export gridpoints, gridcells, celltypes, refcells
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

# Convenience constructors

function RefCell(
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
  data = GridData(points,cells,celltypes)
  UGrid(data,refcells)
end

function Connections(c::Vector{Vector{Int}})
  data, ptrs = generate_data_and_ptrs(c)
  Connections(data,ptrs)
end

function GridData(points,cells::Vector{Vector{Int}},celltypes)
  co = Connections(cells)
  GridData(points,co,celltypes)
end

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

celltypes(grid::Grid) = grid.data.celltypes

refcells(grid::Grid) = grid.refcells

vtkid(grid::RefCell) = grid.vtkid

vtknodes(grid::RefCell) = grid.vtknodes

end #module Core
