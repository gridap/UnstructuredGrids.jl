module CoreTests

using Test
using UnstructuredGrids
using UnstructuredGrids.Core
using UnstructuredGrids.Kernels
using UnstructuredGrids.Factories
using UnstructuredGrids.VTK

c = Connections([[1,2,6,3,],[1,4,4],[1]])
s = """
    1 -> [1, 2, 6, 3]
    2 -> [1, 4, 4]
    3 -> [1]
    """
@test s == string(c)

grid = UGrid(domain=(0,1,-1,0),partition=(2,2))

cell_to_ctype = celltypes(grid)

ctype_to_refcell = refcells(grid)

cell_to_vertices = connections(grid)

c = cell_to_vertices
l = [1, 2, 4, 5, 2, 3, 5, 6, 4, 5, 7, 8, 5, 6, 8, 9]
p = [1, 5, 9, 13, 17]
@test list(c) == l
@test ptrs(c) == p

vertex_to_cells = generate_dual_connections(cell_to_vertices)

c = vertex_to_cells
l = [1, 1, 2, 2, 1, 3, 1, 2, 3, 4, 2, 4, 3, 3, 4, 4]
p = [1, 2, 4, 5, 7, 11, 13, 14, 16, 17]
@test list(c) == l
@test ptrs(c) == p

cell_to_faces = generate_cell_to_faces(
  cell_to_vertices,vertex_to_cells,cell_to_ctype,ctype_to_refcell,1)

c = cell_to_faces
l = [1, 2, 3, 4, 5, 6, 4, 7, 2, 8, 9, 10, 6, 11, 10, 12]
p = [1, 5, 9, 13, 17]
@test list(c) == l
@test ptrs(c) == p

fgrid = UGrid(grid,1)

face_to_vertices = connections(fgrid)

c = face_to_vertices
l = [1, 2, 4, 5, 1, 4, 2, 5, 2, 3, 5, 6, 3, 6, 7, 8, 4, 7, 5, 8, 8, 9, 6, 9]
p = [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25]
@test list(c) == l
@test ptrs(c) == p

fgrid = UGrid(grid,1,vertex_to_cells)

face_to_vertices = connections(fgrid)

c = face_to_vertices
@test list(c) == l
@test ptrs(c) == p

fgrid = UGrid(grid,1,vertex_to_cells,cell_to_faces)

face_to_vertices = connections(fgrid)

c = face_to_vertices
@test list(c) == l
@test ptrs(c) == p

vertex_to_faces = generate_dual_connections(face_to_vertices)

cell_to_faces = generate_cell_to_faces_from_faces(
  cell_to_vertices, vertex_to_faces, cell_to_ctype, ctype_to_refcell, 1)

c = cell_to_faces
l = [1, 2, 3, 4, 5, 6, 4, 7, 2, 8, 9, 10, 6, 11, 10, 12]
p = [1, 5, 9, 13, 17]
@test list(c) == l
@test ptrs(c) == p

face_to_cells = generate_dual_connections(cell_to_faces)

_face_to_isboundary = generate_facet_to_isboundary(face_to_cells)

face_to_isboundary = Bool[
  true, false, true, false, true, false, true, true, true, false, true, true]

@test face_to_isboundary == _face_to_isboundary

_vertex_to_isboundary = generate_face_to_isboundary(
  face_to_isboundary, vertex_to_faces)

vertex_to_isboundary = Bool[
  true, true, true, true, false, true, true, true, true]

@test vertex_to_isboundary == _vertex_to_isboundary

end # module CoreTests
