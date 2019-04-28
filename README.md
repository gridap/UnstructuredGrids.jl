# UnstructuredGrids

*Helper routines for topological operations on unstructured grids in julia*

[![Build Status](https://travis-ci.com/lssc-team/UnstructuredGrids.jl.svg?branch=master)](https://travis-ci.com/lssc-team/UnstructuredGrids.jl)
[![Codecov](https://codecov.io/gh/lssc-team/UnstructuredGrids.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/lssc-team/UnstructuredGrids.jl)

**UnstructuredGrids** provides a set of functions providing common topological operations associated with unstructured meshes/grids such as:

- Find the lower dimensial objects (e.g., edges and faces) on the boundary of each cell in the grid
- Find the vertices on low dimensional objects of the grid (e.g., the vertices on each face, the vertices on each edge)
- Find dual connections (e.g., cells arround a face, cells around a vertex, faces around an edge, etc.)
- Identify objects on the boundary of the grid
- Export unstructured grids into `.vtu` files (using the `WriteVTK` package).

## Quick Start
```julia
#Pkg.add("https://github.com/lssc-team/UnstructuredGrids.jl")
using UnstructuredGrids

# Generate a toy structured grid of hexahedra with 2x3x1 cells
grid = Grid(domain=(0,1,0,1,0,1),partition=(2,3,1))

# Get the connectivities of the cells of the grid (i.e., for each cell, the ids of its vertices)
cell_to_vertices = connections(grid)
```
