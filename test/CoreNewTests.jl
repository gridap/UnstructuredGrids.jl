module CoreNewTests

using UnstructuredGrids.CoreNew: RefCell

faces_1d = [[1,2],[3,4],[1,3],[2,4]]

refcell = RefCell( ndims = 2, dim_to_faces = [faces_1d])

end # module CoreNewTests
