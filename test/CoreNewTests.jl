module CoreNewTests

using UnstructuredGrids.CoreNew: RefCell
using UnstructuredGrids.CoreNew: VERTEX
using UnstructuredGrids.CoreNew: AbstractRefCell

const SEGMENT = RefCell(
  ndims = 1,
  faces = [ [[1],[2]] ],
  reffaces = [ VERTEX ],
  facetypes = [ [1,1] ],
  points = [-1 1;],
  vtkid = 3,
  vtknodes = [1,2] )

const SQUARE = RefCell(
  ndims = 2,
  faces = [ [[1],[2],[3],[4]], [[1,2],[3,4],[1,3],[2,4]] ],
  reffaces = [ VERTEX, SEGMENT ],
  facetypes = [ [1,1,1,1], [2,2,2,2] ],
  points = Float64[ -1 1 -1 1; -1 -1 1 1],
  vtkid = 9,
  vtknodes = [1,2,4,3])

const HEXAHEDRON = RefCell(
  ndims = 3,
  faces = [
    [[1],[2],[3],[4],[5],[6],[7],[8]],
    [[1,2],[3,4],[1,3],[2,4],[5,6],[7,8],[5,7],[6,8],[1,5],[2,6],[3,7],[4,8]],
    [[1,2,3,4],[5,6,7,8],[1,2,5,6],[3,4,7,8],[1,3,5,7],[2,4,6,8]]],
  reffaces = [ VERTEX, SEGMENT, SQUARE ],
  facetypes = [ fill(1,8), fill(2,12), fill(3,6) ],
  points = Float64[ -1 1 -1 1 -1 1 -1 1; -1 -1 1 1 -1 -1 1 1; -1 -1 -1 -1 1 1 1 1],
  vtkid = 12,
  vtknodes = [1,2,4,3,5,6,8,7])

end # module CoreNewTests
