
#
#  7 ---9--- 8 ---12-- 9
#
#  |         |         |
#  |         |         |
#  |         |         |
#  10   4    11   5    13
#  |         |         |
#  |         |         |
#  |         |         |
#
#  4 ---2--- 5 ---8--- 6
#
#  |         | \   3   |
#  |         |  \      |
#  |         |   \     |
#  3    1    4    6    7
#  |         |     \   |
#  |         |  2   \  |
#  |         |       \ |
#
#  1 ---1--- 2 ---5--- 3
#


                        #1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8
cell_to_vertices_data = [1,2,4,5,2,3,5,3,6,5,4,5,7,8,5,6,8,9]
cell_to_vertices_ptrs = [1,      5,    8,    11,     15,    19]

                       #1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8
vertex_to_cells_data = [1,1,2,2,3,1,4,1,2,3,4,5,3,5,4,4,5,5]
vertex_to_cells_ptrs = [1,2,  4,  6,  8,        13, 15,16, 18,19]

cell_to_faces_data = [1,2,3,4,5,4,6,7,6,8,2,9,10,11,8,12,11,13]
cell_to_faces_ptrs = [1,5,8,11,15,19]

quad_to_lface_to_lvertices_data = [1,2,3,4,1,3,2,4]
quad_to_lface_to_lvertices_ptrs = [1,3,5,7,9]

tri_to_lface_to_lvertices_data = [1,2,1,3,2,3]
tri_to_lface_to_lvertices_ptrs = [1,3,5,7]

quad = 1
tri = 2
cell_to_ctype = [quad, tri, tri, quad, quad]

ctype_to_lface_to_lvertices_data = [
  quad_to_lface_to_lvertices_data, tri_to_lface_to_lvertices_data]

ctype_to_lface_to_lvertices_ptrs = [
  quad_to_lface_to_lvertices_ptrs, tri_to_lface_to_lvertices_ptrs]

