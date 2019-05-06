# UnstructuredGrids

*Helper routines for topological operations on unstructured grids in julia*

[![Build Status](https://travis-ci.com/lssc-team/UnstructuredGrids.jl.svg?branch=master)](https://travis-ci.com/lssc-team/UnstructuredGrids.jl)
[![Codecov](https://codecov.io/gh/lssc-team/UnstructuredGrids.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/lssc-team/UnstructuredGrids.jl)

If you ❤️ this project, give us a ⭐️!

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

# Generate a toy structured grid of the unit cube
# with 2x3x2 quadrilateral cells
grid = UGrid(domain=(0,1,0,1,0,1),partition=(2,3,2))

# Get the connectivity of the cells of the grid (i.e., for each cell, the ids of its vertices)
cell_to_vertices = connections(grid)

# The cell connectivity is represented by a `Connections` object
@assert isa(cell_to_vertices, Connections)

# `Connections` is a struct representing a vector of vectors
# in compressed form. It is fully described by a vector containing
# the flatted underlying data and another vector containing
# the start and end of each of the sub-vectors.
# That is:

cell_to_vertices_data = list(cell_to_vertices)
cell_to_vertices_ptrs = ptrs(cell_to_vertices)

cell = 2
ibeg = cell_to_vertices_ptrs[cell]
iend = cell_to_vertices_ptrs[cell+1]-1
vertices = cell_to_vertices_data[ibeg:iend] # Vertices of cell 2
@assert vertices == [2, 3, 5, 6, 14, 15, 17, 18]

# Get the coordinates of the mesh vertices
vertex_to_coords = coordinates(grid)
@assert isa(vertex_to_coords, Array{Float64,2})

# First axis for space dimensions, second one for number of vertices
@assert size(vertex_to_coords) == (3,36)

# Find cells around each vertex
vertex_to_cells = generate_dual_connections(cell_to_vertices)

# The returned data is a vector of vectors represented by a Connections object
@assert isa(vertex_to_cells, Connections)

# Generate a global numbering for the edges (1d objects) of the grid,
# and find which edges are on the boundary of each cell
n=1
cell_to_edges = generate_cell_to_faces(n,grid)

# The returned data again a Connections object
@assert isa(cell_to_edges, Connections)

# Idem for faces (2d objects)
n=2
cell_to_faces = generate_cell_to_faces(n,grid)
@assert isa(cell_to_faces, Connections)
```
