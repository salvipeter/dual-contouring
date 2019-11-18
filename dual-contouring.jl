module DualContouring

export isosurface, writeOBJ

ϵ = 1e-15

"""
    compute_cell(origin, dirs, vertices)

Computes an approximate position of the isosurface in a cell. `origin` is the bottom-left-front
corner of the cell, with `dirs` contains 3 axis-aligned vectors going to the opposite corner.
`vertices` gives the function values at the cell corners in this fashion:
```
      7-----8
     /|    /|     Y
    3-----4 |     ^  
    | 5---|-6     | Z
    |/    |/      |/
    1-----2       O----> X
```
It is assumed that the function values are never exactly zero.
When all values have the same sign, the function returns `nothing`.
"""
function compute_cell(origin, dirs, vertices)
    (all(v < 0 for v in vertices) || all(v > 0 for v in vertices)) && return nothing
    edges = [0, 4, 1, 5, 2, 6, 3, 7,
             0, 2, 1, 3, 4, 6, 5, 7,
             0, 1, 2, 3, 4, 5, 6, 7]
    mean = [0, 0, 0]
    count = 0
    for i in 1:12
        v1, v2 = edges[2i-1], edges[2i]
        a, b = vertices[v1+1], vertices[v2+1]
        a * b > 0 && continue
        denom = abs(b - a)
        x = denom < ϵ ? 0.5 : abs(a) / denom
        c = (i - 1) ÷ 4 + 1
        mean += dirs[1] * (v1 ÷ 4) + dirs[2] * (v1 % 4 ÷ 2) + dirs[3] * (v1 % 2) + dirs[c] * x
        count += 1
    end
    origin + mean / count
end

"""
    generate_cells(values, corner, delta, resolution)

Generates cells with dual points for an isosurface. `values` is a 3D matrix of dimension
`resolution[i]+1` in each direction, containing the function values. `corner` is the
bottom-left-front corner of the volume, with the `delta` vector going to the opposite corner.
`resolution` is a tuple of 3 values, describing the number of cells in each axis direction.

The return value is a `(points, cells)` pair, where `points` is a list of approximated pionts
of the isosurface, and `cells` is a 3D matrix of dimension `resolution[i]` in each direction,
containing an index to the `points` array, or 0 when empty.
"""
function generate_cells(values, corner, delta, resolution)
    cells = zeros(Int, resolution)
    points = []
    dirs = [[delta[1], 0, 0],
            [0, delta[2], 0],
            [0, 0, delta[3]]]
    for i in 1:resolution[1], j in 1:resolution[2], k in 1:resolution[3]
        vertices = []
        for di in 0:1, dj in 0:1, dk in 0:1
            push!(vertices, values[i+di,j+dj,k+dk])
        end
        origin = corner + [delta[1] * (i - 1), delta[2] * (j - 1), delta[3] * (k - 1)]
        surface_point = compute_cell(origin, dirs, vertices)
        if surface_point != nothing
            push!(points, surface_point)
            cells[i,j,k] = length(points)
        end
    end

    points, cells
end

"""
    generate_faces(values, cells, resolution)

Generate faces for the cells created by `generate_cells`, when possible. `values` and
`resolution` are defined the same way as in that function.

The return value is an array of quads, where each quad is an array containing 4 indices.
"""
function generate_faces(values, cells, resolution)
    quads = []
    for (r, value, cell) in [([1, 2, 3], (i, j, k) -> values[i,j,k], (i, j, k) -> cells[i,j,k]),
                             ([2, 3, 1], (i, j, k) -> values[j,k,i], (i, j, k) -> cells[j,k,i]),
                             ([3, 1, 2], (i, j, k) -> values[k,i,j], (i, j, k) -> cells[k,i,j])]
        for i in 1:resolution[r[1]], j in 2:resolution[r[2]], k in 2:resolution[r[3]]
            a, b, c, d = cell(i,j,k), cell(i,j-1,k), cell(i,j-1,k-1), cell(i,j,k-1)
            a * b * c * d == 0 && continue
            v1, v2 = value(i,j,k), value(i+1,j,k)
            v1 * v2 > 0 && continue
            if v1 < 0
                push!(quads, [a, b, c, d])
            else
                push!(quads, [d, c, b, a])
            end
        end
    end
    quads
end

"""
    isosurface(f, isolevel, bounding_box, resolution)

Generates the isosurface of `f` at the given `isolevel`. The `bounding_box` is a tuple
containing the bottom-left-front and top-right-back points. The `resolution` is given as a
tuple of 3 integers - the number of cells in each of the axis directions.

Sphere example:

    julia> using LinearAlgebra
    julia> isosurface(p -> norm(p) - 1, 0, ([-2, -2, -2], [2, 2, 2]), (20, 20, 20))
"""
function isosurface(f, isolevel, bounding_box, resolution)
    delta = (bounding_box[2] - bounding_box[1]) ./ resolution
    values = zeros(resolution .+ 1)
    for i in 0:resolution[1], j in 0:resolution[2], k in 0:resolution[3]
        v = f(bounding_box[1] + delta .* [i, j, k]) - isolevel
        values[i+1,j+1,k+1] = abs(v) < ϵ ? copysign(ϵ, v) : v
    end
    points, cells = generate_cells(values, bounding_box[1], delta, resolution)
    faces = generate_faces(values, cells, resolution)
    points, faces
end

"""
    writeOBJ(points, faces, filename)

Writes a quad mesh in OBJ format to the given file.
"""
function writeOBJ(points, faces, filename)
    open(filename, "w") do fp
        for p in points
            println(fp, "v $(p[1]) $(p[2]) $(p[3])")
        end
        for f in faces
            println(fp, "f $(f[1]) $(f[2]) $(f[3]) $(f[4])")
        end
    end
end

end # module
