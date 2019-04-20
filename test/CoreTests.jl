module CoreTests

using Test
using UnstructuredGrids
using UnstructuredGrids.Core
using UnstructuredGrids.Factories
using UnstructuredGrids.VTK

c = Connections([[1,2,6,3,],[1,4,4],[1]])
s = """
    1 -> [1, 2, 6, 3]
    2 -> [1, 4, 4]
    3 -> [1]
    """
@test s == string(c)

grid = Grid(domain=(0,1,-1,0),partition=(2,2))

graph = GridGraph(grid)

c = connections(graph,from=2,to=1)

l = [1, 2, 3, 4, 5, 6, 4, 7, 2, 8, 9, 10, 6, 11, 10, 12]
p = [1, 5, 9, 13, 17]
@test list(c) == l
@test ptrs(c) == p

c = connections(graph,from=2,to=0)

l = [1, 2, 4, 5, 2, 3, 5, 6, 4, 5, 7, 8, 5, 6, 8, 9]
p = [1, 5, 9, 13, 17]
@test list(c) == l
@test ptrs(c) == p

c = connections(graph,from=1,to=2)

l = [1, 1, 3, 1, 1, 2, 2, 2, 4, 2, 3, 3, 3, 4, 4, 4]
p = [1, 2, 4, 5, 7, 8, 10, 11, 12, 13, 15, 16, 17]
@test list(c) == l
@test ptrs(c) == p

c = connections(graph,from=0,to=2)

l = [1, 1, 2, 2, 1, 3, 1, 2, 3, 4, 2, 4, 3, 3, 4, 4]
p = [1, 2, 4, 5, 7, 11, 13, 14, 16, 17]
@test list(c) == l
@test ptrs(c) == p

grid = Grid(domain=(0,1,-1,0,2,3),partition=(2,4,3))

graph = GridGraph(grid)

c = connections(graph,from=3,to=2)
c = connections(graph,from=3,to=1)
c = connections(graph,from=3,to=0)
c = connections(graph,from=2,to=3)
c = connections(graph,from=1,to=3)
c = connections(graph,from=0,to=3)

grid = Grid(domain=(0,1,-1,0),partition=(2,2))

fgrid = Grid(grid,dim=1)

c = connections(fgrid)

l = [1, 2, 4, 5, 1, 4, 2, 5, 2, 3, 5, 6, 3, 6, 7, 8, 4, 7, 5, 8, 8, 9, 6, 9]
p = [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25]

@test list(c) == l
@test ptrs(c) == p

grid = Grid(domain=(0,1,-1,0,2,3),partition=(2,4,3))

fgrid = Grid(grid,dim=2)

fgrid = Grid(grid,dim=1)

HEX = UnstructuredGrids.Factories.HEXAHEDRON

grid = Grid(HEX,dim=2)

grid = Grid(HEX,dim=1)

grid = Grid(HEX,dim=0)

end # module CoreTests
