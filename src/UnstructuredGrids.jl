module UnstructuredGrids

export face_to_cells
export rewind_ptrs!
export length_to_ptrs!

# Functionality given by this module

"""
Given the faces on the boundary of each cell,
find the cells around each face.
"""
function face_to_cells(
  cell_to_faces_data::AbstractVector{<:Integer},
  cell_to_faces_ptrs::AbstractVector{<:Integer},
  nfaces = maximum(cell_to_faces_data))
  _face_to_cells(cell_to_faces_data,cell_to_faces_ptrs,nfaces)
end

"""
Rewind the given vector.
"""
function rewind_ptrs!(ptrs::AbstractVector{<:Integer})
  @inbounds for i in (length(ptrs)-1):-1:1
    ptrs[i+1] = ptrs[i]
  end
  ptrs[1] = 1
end

"""
Given a vector of integers, mutate it from length state to pointer state.
"""
function length_to_ptrs!(ptrs::AbstractArray{<:Integer})
  ptrs[1] = 1
  @inbounds for i in 1:(length(ptrs)-1)
    ptrs[i+1] += ptrs[i]
  end
end

# Helpers

function _face_to_cells(cell_to_faces_data, cell_to_faces_ptrs, nfaces)

  face_to_cells_ptrs = zeros(Int,nfaces+1)

  _face_to_cells_count!(
    face_to_cells_ptrs, cell_to_faces_data, cell_to_faces_ptrs)

  length_to_ptrs!(face_to_cells_ptrs)

  ndata = face_to_cells_ptrs[end]-1

  face_to_cells_data = Vector{Int}(undef,ndata)

  _face_to_cells_fill!(
    face_to_cells_data, face_to_cells_ptrs,
    cell_to_faces_data, cell_to_faces_ptrs)

  rewind_ptrs!(face_to_cells_ptrs)

  (face_to_cells_data, face_to_cells_ptrs)

end

function _face_to_cells_count!(
  face_to_cells_ptrs, cell_to_faces_data, cell_to_faces_ptrs)

  ncells = length(cell_to_faces_ptrs) - 1
  for cell in 1:ncells
    a = cell_to_faces_ptrs[cell]
    b = cell_to_faces_ptrs[cell+1]-1
    for p in a:b
      face = cell_to_faces_data[p]
      face_to_cells_ptrs[face+1] += 1
    end
  end

end

function _face_to_cells_fill!(
  face_to_cells_data, face_to_cells_ptrs,
  cell_to_faces_data, cell_to_faces_ptrs)

  ncells = length(cell_to_faces_ptrs) - 1
  for cell in 1:ncells
    a = cell_to_faces_ptrs[cell]
    b = cell_to_faces_ptrs[cell+1]-1
    for p in a:b
      face = cell_to_faces_data[p]
      q = face_to_cells_ptrs[face]
      face_to_cells_data[q] = cell
      face_to_cells_ptrs[face] += 1
    end
  end

end

end # module
