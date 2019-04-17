module ReferenceGridsTests

using Test
using UnstructuredGrids.ReferenceGrids

@test SEGMENT == SEGMENT

@show SEGMENT.point_to_coords
@show SEGMENT.cell_to_points

end # module ReferenceGridsTests
