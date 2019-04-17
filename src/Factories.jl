module Factories

using UnstructuredGrids.Helpers
using UnstructuredGrids.Kernels
using UnstructuredGrids.Core

# **DISCLAIMER**
# This library is not supposed to be a mesh generator.
# The following mesh generation routines are mainly for
# testing purposes.

export generate

function generate(;domain,partition)
  _cartesian_grid(domain,partition)
end

# Helpers

points = zeros(0,1)
cells = [Int[]]
celltypes = Int[]
refcells = Vector{RefCell}(undef,0)
vtkid = 1
vtknodes = [1]

const VERTEX = RefCell(points,cells,celltypes,refcells,vtkid,vtknodes)

points = Float64[ -1 1; ]
cells = [[1,],[2,]]
celltypes = [1,1]
refcells = [VERTEX]
vtkid = 3
vtknodes = [1,2]

const SEGMENT = RefCell(points,cells,celltypes,refcells,vtkid,vtknodes)

points = Float64[ 0 1 0; 0 0 1]
cells = [[1,2],[2,3],[3,1]]
celltypes = [1,1,1]
refcells = [SEGMENT]
vtkid = 5
vtknodes = [1,2,3]

const TRIANGLE = RefCell(points,cells,celltypes,refcells,vtkid,vtknodes)

points = Float64[ -1 1 -1 1; -1 -1 1 1]
cells = [[1,2],[3,4],[1,3],[2,4]]
celltypes = [1,1,1,1]
refcells = [SEGMENT]
vtkid = 9
vtknodes = [1,2,4,3]

const SQUARE = RefCell(points,cells,celltypes,refcells,vtkid,vtknodes)

_cartesian_grid(domain,partition) = @notimplemented

function _cartesian_grid(domain::NTuple{4,Int},partition::NTuple{2,Int})
  ncells = prod(partition)
  npoints = prod([ n+1 for n in partition])
  d = 2
  n = 4
  points = Array{Float64,2}(undef,(d,npoints))
  celltypes = ones(Int,ncells)
  cellptrs = fill(n,ncells+1)
  length_to_ptrs!(cellptrs)
  celldata = Vector{Int}(undef,n*ncells)
  refcells = [SQUARE]
  _cartesian_2d_fill_points!(points,domain,partition)
  _cartesian_2d_fill_cells!(celldata,partition)
  UGrid(points,celldata,cellptrs,celltypes,refcells)
end

function _cartesian_2d_fill_points!(points,domain,partition)
  ncx = partition[1]
  ncy = partition[2]
  x0 = domain[1]
  x1 = domain[2]
  y0 = domain[3]
  y1 = domain[4]
  dx = x1-x0/ncx
  dy = y1-y0/ncy
  k = 1
  for j in 1:ncy+1
    for i in 1:ncx+1
      points[1,k] = x0 + (i-1)*dx
      points[2,k] = y0 + (j-1)*dy
      k += 1
    end
  end
end

function _cartesian_2d_fill_cells!(celldata,partition)
  ncx = partition[1]
  ncy = partition[2]
  k = 1
  for j in 1:ncy
    for i in 1:ncx
      for b in 0:1
        for a in 0:1
          celldata[k] =  i+a + (j+b-1)*(ncx+1)
          k += 1
        end
      end
    end
  end
end

end # module Factories
