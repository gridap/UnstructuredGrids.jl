module Kernels

export face_to_cells
export cell_to_faces
export rewind_ptrs!
export length_to_ptrs!
export generate_data_and_ptrs

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
Given the vertices on the boundary of each cell,
and, for each cell type, the local vertices on
the boundary of each local face, returns the faces
on the boundary of each cell
"""
function cell_to_faces(
  cell_to_vertices_data::AbstractVector{<:Integer},
  cell_to_vertices_ptrs::AbstractVector{<:Integer},
  ctype_to_lface_to_lvertices_data::AbstractVector{<:AbstractVector{<:Integer}},
  ctype_to_lface_to_lvertices_ptrs::AbstractVector{<:AbstractVector{<:Integer}},
  cell_to_ctype::AbstractVector{<:Integer},
  vertex_to_cells_data::AbstractVector{<:Integer},
  vertex_to_cells_ptrs::AbstractVector{<:Integer})
  _cell_to_faces(
    cell_to_vertices_data,
    cell_to_vertices_ptrs,
    ctype_to_lface_to_lvertices_data,
    ctype_to_lface_to_lvertices_ptrs,
    cell_to_ctype,
    vertex_to_cells_data,
    vertex_to_cells_ptrs)
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

"""
Given a vector of vectors compute the corresponding data and and ptrs
"""
function generate_data_and_ptrs(vv::Vector{Vector{Int}})
  ptrs = Vector{Int}(undef,length(vv)+1)
  _generate_data_and_ptrs_fill_ptrs!(ptrs,vv)
  length_to_ptrs!(ptrs)
  ndata = ptrs[end]-1
  data = Vector{Int}(undef,ndata)
  _generate_data_and_ptrs_fill_data!(data,vv)
  (data, ptrs)
end

# Helpers

const UNSET = 0

function _generate_data_and_ptrs_fill_ptrs!(ptrs,vv)
  for (i,v) in enumerate(vv)
    ptrs[i+1] = length(v)
  end
end

function _generate_data_and_ptrs_fill_data!(data,vv)
  k = 1
  for v in vv
    for vi in v
      data[k] = vi
      k += 1
    end
  end
end

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

function _cell_to_faces(
    cell_to_vertices_data,
    cell_to_vertices_ptrs,
    ctype_to_lface_to_lvertices_data,
    ctype_to_lface_to_lvertices_ptrs,
    cell_to_ctype,
    vertex_to_cells_data,
    vertex_to_cells_ptrs)

  ncells = length(cell_to_vertices_ptrs) - 1

  cell_to_faces_ptrs = zeros(Int,ncells+1)

  type_to_nlfaces = _type_to_nlfaces(ctype_to_lface_to_lvertices_ptrs)

  _cell_to_faces_count!(
    cell_to_faces_ptrs,
    type_to_nlfaces,
    cell_to_ctype)

  length_to_ptrs!(cell_to_faces_ptrs)

  ndata = cell_to_faces_ptrs[end]-1

  cell_to_faces_data = fill(UNSET,ndata)

  nvertices = max_nvertices_in_lface(ctype_to_lface_to_lvertices_ptrs)
  vertices = fill(UNSET,nvertices)
  vertices_scratch = fill(UNSET,nvertices)

  nvertices = max_cells_arround_vertex(vertex_to_cells_ptrs)
  cells_around = fill(UNSET,nvertices)
  cells_around_scratch = fill(UNSET,nvertices)
  lface_of_cells_around = fill(UNSET,nvertices)

  _cell_to_faces_fill!(
    cell_to_faces_data,
    cell_to_faces_ptrs,
    cell_to_vertices_data,
    cell_to_vertices_ptrs,
    ctype_to_lface_to_lvertices_data,
    ctype_to_lface_to_lvertices_ptrs,
    cell_to_ctype,
    vertex_to_cells_data,
    vertex_to_cells_ptrs,
    vertices,
    vertices_scratch,
    cells_around,
    cells_around_scratch,
    lface_of_cells_around)

  (cell_to_faces_data, cell_to_faces_ptrs)

end

function max_nvertices_in_lface(ctype_to_lface_to_lvertices_ptrs)
  n = 0
  for lface_to_lvertices_ptrs in ctype_to_lface_to_lvertices_ptrs
    nlfaces = length(lface_to_lvertices_ptrs)-1
    for lface in 1:nlfaces
      m = lface_to_lvertices_ptrs[lface+1] - lface_to_lvertices_ptrs[lface]
      n = max(n,m)
    end
  end
  n
end

function max_cells_arround_vertex(vertex_to_cells_ptrs)
  n = 0
  nvertices = length(vertex_to_cells_ptrs)-1
  for vertex in 1:nvertices
    m = vertex_to_cells_ptrs[vertex+1]-vertex_to_cells_ptrs[vertex]
    n = max(n,m)
  end
  n
end

function _type_to_nlfaces(ctype_to_lface_to_lvertices_ptrs)
  [ length(v)-1 for v in ctype_to_lface_to_lvertices_ptrs ]
end

function _cell_to_faces_count!(
    cell_to_faces_ptrs,
    type_to_nlfaces,
    cell_to_ctype)

  ncells = length(cell_to_ctype)
  for cell in 1:ncells
    ctype = cell_to_ctype[cell]
    nlfaces = type_to_nlfaces[ctype]
    cell_to_faces_ptrs[cell+1] = nlfaces
  end

end

function  _cell_to_faces_fill!(
    cell_to_faces_data,
    cell_to_faces_ptrs,
    cell_to_vertices_data,
    cell_to_vertices_ptrs,
    ctype_to_lface_to_lvertices_data,
    ctype_to_lface_to_lvertices_ptrs,
    cell_to_ctype,
    vertex_to_cells_data,
    vertex_to_cells_ptrs,
    vertices,
    vertices_scratch,
    cells_around,
    cells_around_scratch,
    lface_of_cells_around)

  face = 1

  ncells = length(cell_to_ctype)

  for cell in 1:ncells

    ctype = cell_to_ctype[cell]
    lface_to_lvertices_data = ctype_to_lface_to_lvertices_data[ctype]
    lface_to_lvertices_ptrs = ctype_to_lface_to_lvertices_ptrs[ctype]
    nlfaces = length(lface_to_lvertices_ptrs)-1
    a = cell_to_faces_ptrs[cell]-1

    for lface in 1:nlfaces

      if cell_to_faces_data[a+lface] != UNSET
        continue
      end

      _fill_vertices_in_lface!(
        vertices,
        lface,
        lface_to_lvertices_data,
        lface_to_lvertices_ptrs,
        cell,
        cell_to_vertices_data,
        cell_to_vertices_ptrs)

      _find_cells_around_vertices!(
        cells_around,
        cells_around_scratch,
        vertices,
        vertex_to_cells_data,
        vertex_to_cells_ptrs)

      _find_lface_of_cells_around_lface!(
        lface_of_cells_around,
        cells_around,
        lface,
        ctype_to_lface_to_lvertices_data,
        ctype_to_lface_to_lvertices_ptrs,
        cell,
        cell_to_vertices_data,
        cell_to_vertices_ptrs,
        cell_to_ctype,
        vertices,
        vertices_scratch)

      _fill_face_in_cells_arround!(
        cell_to_faces_data,
        cell_to_faces_ptrs,
        face,
        cells_around,
        lface_of_cells_around)

      face += 1
    end

  end

end

function _fill_face_in_cells_arround!(
  cell_to_faces_data,
  cell_to_faces_ptrs,
  face,
  cells_around,
  lface_of_cells_around)

  for icell_around in 1:length(cells_around)
    cell_around = cells_around[icell_around]
    if cell_around == UNSET
      continue
    end
    lface = lface_of_cells_around[icell_around]
    f = cell_to_faces_ptrs[cell_around]-1
    cell_to_faces_data[f+lface] = face
  end

end

function _find_cells_around_vertices!(
  cells_around,
  cells_around_scratch,
  vertices,
  vertex_to_cells_data,
  vertex_to_cells_ptrs)

  ncells_around = UNSET
  ncells_around_scratch = UNSET

  for ivertex in 1:length(vertices)
    vertex = vertices[ivertex]
    if vertex == UNSET
      continue
    end
    if ivertex == 1
      ncells_around = _fill_cells_around_scratch!(
        cells_around,
        vertex,
        vertex_to_cells_data,
        vertex_to_cells_ptrs)
    else
      ncells_around_scratch = _fill_cells_around_scratch!(
        cells_around_scratch,
        vertex,
        vertex_to_cells_data,
        vertex_to_cells_ptrs)
      _set_intersection!(
        cells_around,cells_around_scratch,
        ncells_around,ncells_around_scratch)
    end
  end

end

function _fill_cells_around_scratch!(
  cells_around_scratch,
  vertex,
  vertex_to_cells_data,
  vertex_to_cells_ptrs)

  cells_around_scratch .= UNSET
  d = vertex_to_cells_ptrs[vertex] - 1
  ncells_around = vertex_to_cells_ptrs[vertex+1] - (d + 1)
  for icell_around in 1:ncells_around
    cell_around = vertex_to_cells_data[d+icell_around]
    cells_around_scratch[icell_around] = cell_around
  end
  ncells_around

end

function _set_intersection!(
  cells_around,cells_around_scratch,
  ncells_around,ncells_around_scratch)
  for i in 1:ncells_around
    if cells_around[i] == UNSET
      continue
    end
    _find_eq!(i,cells_around,cells_around_scratch,ncells_around_scratch)
  end
end

function _find_eq!(i,cells_around,cells_around_scratch,ncells_around_scratch)
  for j in 1:ncells_around_scratch
    if cells_around[i] == cells_around_scratch[j]
      return
    end
  end
  cells_around[i] = UNSET
  return
end

function _find_lface_of_cells_around_lface!(
  lface_of_cells_around,
  cells_around,
  lface,
  ctype_to_lface_to_lvertices_data,
  ctype_to_lface_to_lvertices_ptrs,
  cell,
  cell_to_vertices_data,
  cell_to_vertices_ptrs,
  cell_to_ctype,
  vertices,
  vertices_scratch)

  lface_of_cells_around .= UNSET

  for icell_around in 1:length(cells_around)
    cell_around = cells_around[icell_around]
    if cell_around == UNSET
      continue
    end

    lface_of_cell_around = _find_lface_with_same_vertices(
      vertices,
      vertices_scratch,
      cell_around,
      cell_to_vertices_data,
      cell_to_vertices_ptrs,
      cell_to_ctype,
      ctype_to_lface_to_lvertices_data,
      ctype_to_lface_to_lvertices_ptrs)

    lface_of_cells_around[icell_around] = lface_of_cell_around

  end

end

function _find_lface_with_same_vertices(
  vertices,
  vertices_scratch,
  cell,
  cell_to_vertices_data,
  cell_to_vertices_ptrs,
  cell_to_ctype,
  ctype_to_lface_to_lvertices_data,
  ctype_to_lface_to_lvertices_ptrs)

  ctype = cell_to_ctype[cell]
  lface_to_lvertices_data = ctype_to_lface_to_lvertices_data[ctype]
  lface_to_lvertices_ptrs = ctype_to_lface_to_lvertices_ptrs[ctype]

  nlfaces = length(lface_to_lvertices_ptrs)-1
  c = cell_to_vertices_ptrs[cell]-1

  for lface in 1:nlfaces

    _fill_vertices_in_lface!(
      vertices_scratch,
      lface,
      lface_to_lvertices_data,
      lface_to_lvertices_ptrs,
      cell,
      cell_to_vertices_data,
      cell_to_vertices_ptrs)

    if _set_equal(vertices,vertices_scratch)
      return lface
    end

  end

  return UNSET

end

function _fill_vertices_in_lface!(
  vertices,
  lface,
  lface_to_lvertices_data,
  lface_to_lvertices_ptrs,
  cell,
  cell_to_vertices_data,
  cell_to_vertices_ptrs)

  vertices .= UNSET
  b = lface_to_lvertices_ptrs[lface] - 1
  c = cell_to_vertices_ptrs[cell] - 1
  nlfvertex = lface_to_lvertices_ptrs[lface+1] - (b + 1)
  for lfvertex in 1:nlfvertex
    lvertex = lface_to_lvertices_data[b+lfvertex]
    vertex = cell_to_vertices_data[c+lvertex]
    vertices[lfvertex] = vertex
  end

end

function _set_equal(vertices,vertices_scratch)

  b = _is_subset(vertices,vertices_scratch)
  if b == false; return false; end
  b = _is_subset(vertices_scratch,vertices)
  if b == false; return false; end
  return true

end

function _is_subset(vertices,vertices_scratch)
  for i in 1:length(vertices)
    v = vertices[i]
    if v == UNSET
      continue
    end
    b = _find_eq(v,vertices_scratch)
    if b == false; return false; end
  end
  return true
end

function _find_eq(v,vertices_scratch)
  for vs in vertices_scratch
    if v == vs
      return true
    end
  end
  return false
end


end # module Kernels