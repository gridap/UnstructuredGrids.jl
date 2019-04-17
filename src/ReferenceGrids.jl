module ReferenceGrids

using UnstructuredGrids

export VERTEX
export SEGMENT

x = zeros(0,1)
c = [Int[]]
t = Int[]
g = Int[]
const VERTEX = UnstructuredGrid(x,c,t,g)

x = [ -1 1; ]
c = [[1,],[2,]]
t = [1,1]
g = [VERTEX]
const SEGMENT = UnstructuredGrid(x,c,t,g)

x = [ 0 1 0; 0 0 1]
c = [[1,2],[2,3],[3,1]]
t = [1,1,1]
g = [SEGMENT]
const TRIANGLE = UnstructuredGrid(x,c,t,g)

x = [ -1 1 -1 1; -1 -1 1 1]
c = [[1,2],[3,4],[1,3],[2,4]]
t = [1,1,1,1]
g = [SEGMENT]
const SQUARE = UnstructuredGrid(x,c,t,g)

end # module ReferenceGrids
