# The model on a curvilinear (terrain-following) mesh of GeophysicalModelGenerator (curvilinear_grid): the grid levels
# follow the topography and the interfaces that span the whole domain, so that these are grid surfaces instead of
# staircases of a rectilinear grid.

"""
Names of the surfaces followed by the levels of the curvilinear mesh, from the top to the bottom (see
[`meshing_surfaces`](@ref)).
"""
const CURVILINEAR_SURFACES = (:topography, :fine_base, :base_cs, :conrad, :moho, :bottom)

"""
    meshing_surfaces(I, zbot; fine_depth = 1.0, fine_min = 0.25, depths = (fine_depth,),
                     interfaces = (:base_cs, :conrad, :moho)) -> NamedTuple of matrices (km)

The surfaces that the levels of the curvilinear mesh follow, from the top to the bottom, from the interfaces `I` of
[`model_interfaces`](@ref):
- `topography`: the free surface (swissALTIRegio);
- `fine_base` (and `below_2`, `below_3`, … for more `depths`): surfaces parallel to the topography, `depths` below it,
  bounding near-surface layers with their own vertical spacing. They are not geological interfaces: the first keeps the
  Quaternary fill of the valleys (at most ≈ 1 km thick) in the finely spaced layer. With `:base_cs` among the
  `interfaces`, surface `i` of `n` lies at most the fraction `i / (n + 1)` of the way down to the base of the
  consolidated sediments, but at least `i fine_min` below the topography (so that cells of these layers are never
  thinner than `fine_min / nz`);
- the geological `interfaces`, any of `base_cs` (base of the consolidated sediments, 3D-ALPS / swissBEDROCK), `conrad`
  (upper / lower crust boundary) and `moho` (European Moho, Spada et al. 2013);
- `bottom`: the flat base of the model at `zbot`.

These interfaces are continuous over the region and do not cross, so each is a grid surface. The base of the
Quaternary deposits (absent in most columns) and the Adriatic Moho (a wedge in the south-east that crosses the Conrad)
are not: they are material boundaries inside the layers, as on the rectilinear grid. Where the consolidated sediments
pinch out (crystalline massifs), the layers above `base_cs` keep a thickness of at least `fine_min` (plus
`min_thickness` per cell, see [`curvilinear_levels`](@ref)): the meshed `base_cs` then lies up to ≈ 0.3 km below the
interface, which the material rules still honour at every node.

Without `:base_cs`, the surfaces parallel to the topography must stay above the Conrad (or the next interface): an error
is thrown where they do not. This is the choice for wave simulations (SeismicWavePC), whose vertical spacing should
change only gradually between layers that are not separated by a material jump along the whole width.
"""
function meshing_surfaces(I, zbot; fine_depth = 1.0, fine_min = 0.25, depths = (fine_depth,),
                          interfaces = (:base_cs, :conrad, :moho))
    all(in((:base_cs, :conrad, :moho)), interfaces) || throw(ArgumentError("interfaces: any of :base_cs, :conrad, :moho"))
    issorted(depths) && all(>(0), depths) || throw(ArgumentError("depths must be positive and increasing"))
    n = length(depths)
    names = [i == 1 ? :fine_base : Symbol("below_", i) for i in 1:n]
    helper = map(1:n) do i
        if :base_cs in interfaces
            I.topo .- max.(min.(depths[i], (I.topo .- I.base_cs) .* (i / (n + 1))), i * fine_min)
        else
            I.topo .- depths[i]
        end
    end
    geo = (; base_cs = I.base_cs, conrad = I.conrad, moho = I.moho_eu)
    inter = [k for k in (:base_cs, :conrad, :moho) if k in interfaces]
    if !isempty(inter) && !(:base_cs in interfaces)
        below = geo[first(inter)]
        bad = count(helper[end] .<= below)
        bad == 0 || throw(ArgumentError("the surface $(depths[end]) km below the topography lies below the $(first(inter)) in $(bad) columns"))
    end
    return NamedTuple{(:topography, names..., inter..., :bottom)}((I.topo, helper..., (geo[k] for k in inter)...,
                                                                   fill(Float64(zbot), size(I.topo))))
end

# CartData surface of size (nx, ny, 1) for curvilinear_grid
function _surface(x, y, h)
    X, Y, _ = xyz_grid(x, y, 0.0)
    Z = reshape(Float64.(h), length(x), length(y), 1)
    return CartData(X, Y, Z, (z = Z,))
end

"""
    curvilinear_levels(x, y, S; dz = (0.025, 0.1, 0.5, 0.5, 1.0), nz = nothing, spacing = :uniform, min_thickness = 0.001)
        -> (; z, nz, levels)

The heights `z[i, j, k]` (km, increasing with `k`) of the curvilinear mesh with the columns `x`, `y` whose levels follow
the surfaces `S` of [`meshing_surfaces`](@ref), made by GeophysicalModelGenerator's `curvilinear_grid`. Layer `l`
(between surface `l` and `l + 1`) has `nz[l]` cells in every column: given, or such that the spacing is `dz[l]` where
the layer is thickest. `levels` gives the level index of every surface. `spacing = :sbp` places the points of each
layer at the nodes of a summation-by-parts grid (SeismicWavePC.jl). Where a layer pinches out, the surface below is
moved down such that each of its cells is at least `min_thickness` (km) thick. Every column is made on its own: the
heights of a part of the columns (with the same `nz`) equal that part of the whole mesh.
"""
function curvilinear_levels(x, y, S; dz = (0.025, 0.1, 0.5, 0.5, 1.0), nz = nothing, spacing = :uniform, min_thickness = 0.001)
    surfs = values(S)
    if nz === nothing
        length(dz) == length(surfs) - 1 || throw(ArgumentError("give one dz per layer ($(length(surfs) - 1))"))
        nz = layer_cells(S, dz; spacing)
    end
    length(nz) == length(surfs) - 1 || throw(ArgumentError("give one nz per layer ($(length(surfs) - 1))"))
    G = curvilinear_grid(x, y, Tuple(_surface(x, y, s) for s in surfs), nz; spacing, min_thickness)
    z = G.z.val
    levels = NamedTuple{keys(S)}(Tuple(size(z, 3) .- cumsum([0; nz])))
    return (; z, nz = collect(nz), levels)
end

"Cells per layer of the surfaces `S`, such that the vertical spacing is `dz[l]` where layer `l` is thickest."
function layer_cells(S, dz; spacing = :uniform)
    surfs = values(S)
    min_cells = spacing == :sbp ? 12 : 1
    return [max(min_cells, ceil(Int, maximum(surfs[l] .- surfs[l + 1]) / dz[l] - 1e-9)) for l in eachindex(dz)]
end

# bilinear interpolation of the raster A (cell centres xr, yr; clamped at the edges) at the points x, y
function _bilinear_raster(xr, yr, A, x, y)
    ix(v, g) = (q = clamp(searchsortedlast(g, v), 1, length(g) - 1); (q, clamp((v - g[q]) / (g[q + 1] - g[q]), 0, 1)))
    return [(((p, s), (q, t)) = (ix(a, xr), ix(b, yr));
             (1 - s) * (1 - t) * A[p, q] + s * (1 - t) * A[p + 1, q] + (1 - s) * t * A[p, q + 1] + s * t * A[p + 1, q + 1])
            for a in x, b in y]
end

"""
    valais_curvilinear_maps(r = valais_region(); dx = 0.25, dz = (0.025, 0.1, 0.5, 0.5, 1.0), fine_depth = 1.0,
                            fine_min = 0.25, depths = (fine_depth,), interfaces = (:base_cs, :conrad, :moho),
                            columns = :uniform, dx_data = dx, spacing = :uniform, min_thickness = 0.001,
                            interface_above = false, heights = true, morcles = true, glaciers = false,
                            lakes = false) -> NamedTuple

The 3D model on a curvilinear mesh, as maps: the column axes `x`, `y` (LV95, km), the heights `z` of the mesh (an array
`(nx, ny, nz)`, see [`curvilinear_levels`](@ref)), the layer `interfaces` ([`model_interfaces`](@ref)), the meshed
`surfaces` ([`meshing_surfaces`](@ref) with `fine_depth`, `fine_min`, `depths`, `interfaces`), the level of each
surface (`levels`), the cells per layer (`nz`) and the Morcles model. The default `dz` gives 25 m in the near-surface
layer (`fine_depth` km), 100 m down to the base of the consolidated sediments, 500 m in the crust and 1 km in the mantle.

- `columns = :uniform`: the columns are the cell centres of the rasters at spacing `dx`. `columns = :sbp`: the nodes of
  a summation-by-parts grid from edge to edge of the region (`sbp_nodes`, spacing `dx` inside), as SeismicWavePC's
  curvilinear solver needs; the rasters are read at `dx_data` and interpolated bilinearly onto the columns.
- `interface_above = true`: a node exactly on the base of the consolidated sediments, the Conrad or the Moho takes the
  material above it (SeismicWavePC's convention for the interfaces of a curvilinear mesh); by default the one below.
- `glaciers`, `lakes`: ice and lake water between the topography and the solid surface, as in
  [`valais_model_maps`](@ref). The top level stays the topography, so they lie inside the near-surface layer.
- `heights = false`: `z` is `nothing` (for meshes too large for one process): [`valais_model_block`](@ref) and
  [`curvilinear_heights`](@ref) make the heights of the columns they need.

[`valais_model_block`](@ref) fills any part of it with the same rules as the rectilinear grid; the mesh has no air
nodes: the top level is the topography.
"""
function valais_curvilinear_maps(r::ValaisRegion = valais_region(); dx = 0.25, dz = (0.025, 0.1, 0.5, 0.5, 1.0),
                                 fine_depth = 1.0, fine_min = 0.25, depths = (fine_depth,), interfaces = (:base_cs, :conrad, :moho),
                                 columns = :uniform, dx_data = dx, spacing = :uniform, min_thickness = 0.001,
                                 interface_above = false, heights = true, morcles = true, glaciers = false, lakes = false)
    columns in (:uniform, :sbp) || throw(ArgumentError("columns = :uniform or :sbp"))
    if columns === :uniform
        x, y, topo = fetch_topography(r; dx = dx)
        _, _, _, tmud = fetch_swissbedrock(r; dx = dx)
        B = ice_water_thickness(r, topo; dx, glaciers, lakes)
    else
        xr, yr, tr = fetch_topography(r; dx = dx_data)
        _, _, _, qr = fetch_swissbedrock(r; dx = dx_data)
        Br = ice_water_thickness(r, tr; dx = dx_data, glaciers, lakes)
        x = sbp_nodes(r.E[1], r.E[2], round(Int, (r.E[2] - r.E[1]) / dx) + 1)
        y = sbp_nodes(r.N[1], r.N[2], round(Int, (r.N[2] - r.N[1]) / dx) + 1)
        topo = _bilinear_raster(xr, yr, replace(tr, NaN => 0.0), x, y)
        tmud = _bilinear_raster(xr, yr, qr, x, y)
        B = map(a -> _bilinear_raster(xr, yr, a, x, y), Br)
    end
    I = model_interfaces(x, y, topo, tmud; B...)
    S = meshing_surfaces(I, r.z[1]; fine_depth, fine_min, depths, interfaces)
    nz = layer_cells(S, dz; spacing)
    L = heights ? curvilinear_levels(x, y, S; nz, spacing, min_thickness) : nothing
    levels = NamedTuple{keys(S)}(Tuple(sum(nz) + 1 .- cumsum([0; nz])))
    @info "Curvilinear mesh: $(length(x)) × $(length(y)) × $(sum(nz) + 1) = $(round(length(x) * length(y) * (sum(nz) + 1) / 1e6, digits = 1)) M nodes; cells per layer $(nz)"
    return (; x, y, z = heights ? L.z : nothing, interfaces = I, surfaces = S, levels, nz, spacing, min_thickness,
            interface_above, morcles = morcles ? read_morcles() : nothing)
end

"""
    curvilinear_heights(maps, ix = eachindex(maps.x), iy = eachindex(maps.y)) -> Array (length(ix), length(iy), nz)

The heights of the curvilinear mesh of [`valais_curvilinear_maps`](@ref) in the columns `maps.x[ix]`, `maps.y[iy]`,
made from the meshed surfaces of these columns: identical to that part of the heights of the whole mesh.
"""
function curvilinear_heights(maps, ix = eachindex(maps.x), iy = eachindex(maps.y))
    S = map(s -> s[ix, iy], maps.surfaces)
    return curvilinear_levels(maps.x[ix], maps.y[iy], S; nz = maps.nz, spacing = maps.spacing, min_thickness = maps.min_thickness).z
end

"""
    build_valais_curvilinear_model(r = valais_region(); dx = 0.25, dz = (0.025, 0.1, 0.5, 0.5, 1.0), fine_depth = 1.0,
                                   spacing = :uniform, min_thickness = 0.001, velocities = DEFAULT_VELOCITIES,
                                   morcles = true, temperature = true, glaciers = false, lakes = false,
                                   outdir = datadir("output"), write = true)
        -> (; model, surfaces, maps)

The Valais model of [`build_valais_model`](@ref) on the curvilinear mesh of [`valais_curvilinear_maps`](@ref): `model` is
a `CartData` whose `z` follows the topography, the base of the consolidated sediments, the Conrad and the European Moho,
with the fields `Phases`, `Vp`, `Vs`, `Rho`, `MorclesUnit` (and `Temperature`) at its nodes. `surfaces` is a `CartData`
with the height of each meshed surface. With `write`, both are written to ParaView (`valais_curvilinear_model.vts`,
`valais_curvilinear_surfaces.vts`).
"""
function build_valais_curvilinear_model(r::ValaisRegion = valais_region(); dx = 0.25, dz = (0.025, 0.1, 0.5, 0.5, 1.0),
                                        fine_depth = 1.0, fine_min = 0.25, spacing = :uniform, min_thickness = 0.001,
                                        velocities = DEFAULT_VELOCITIES, morcles = true, temperature = true,
                                        glaciers = false, lakes = false, outdir = datadir("output"), write = true)
    maps = valais_curvilinear_maps(r; dx, dz, fine_depth, fine_min, spacing, min_thickness, morcles, glaciers, lakes)
    fields = valais_model_block(maps; velocities, temperature)
    X, Y, _ = xyz_grid(maps.x, maps.y, 1:size(maps.z, 3))
    model = CartData(X, Y, maps.z, fields)
    Xs, Ys, _ = xyz_grid(maps.x, maps.y, 0.0)
    s3(a) = reshape(a, size(a)..., 1)
    names = keys(maps.levels)
    surfaces = CartData(Xs, Ys, s3(maps.z[:, :, end]),
                        NamedTuple{names}(Tuple(s3(maps.z[:, :, maps.levels[n]]) for n in names)))
    if write
        mkpath(outdir)
        write_paraview(model, joinpath(outdir, "valais_curvilinear_model"))
        write_paraview(surfaces, joinpath(outdir, "valais_curvilinear_surfaces"))
    end
    return (; model, surfaces, maps)
end

"""
    curvilinear_density_model(r = valais_region(); pad = 20.0, dx = 1.0, sub = 4, nsub_z = 4,
                              dz = (0.25, 0.5, 2.0, 2.5, 5.0), fine_depth = 1.0, fine_min = 0.25, min_thickness = 0.001,
                              glaciers = false, lakes = false)
        -> (; x, y, z, rho, levels, nz, topo, maps)

Density [kg/m³] of the model at the nodes of a coarse curvilinear mesh for gravity: columns every `dx` km over the region
`r` extended by `pad` km on every side, levels following the surfaces of [`meshing_surfaces`](@ref) with the spacings `dz`.
The density of a node is the mean over its dual cell: `sub × sub` columns of the model at `dx / sub` and `nsub_z` depths
in each, at the depths below the topography of the dual cell (so that every sample is below the surface of its own
column). In the padding, the model has no Quaternary deposits (no swissBEDROCK); topography, 3D-ALPS and Moho are used
as in the region. `maps` holds the interfaces at `dx / sub` over the padded region. `glaciers`, `lakes`: ice and lake water
over the whole padded region (see [`ice_water_thickness`](@ref)).
"""
function curvilinear_density_model(r::ValaisRegion = valais_region(); pad = 20.0, dx = 1.0, sub = 4, nsub_z = 4,
                                   dz = (0.25, 0.5, 2.0, 2.5, 5.0), fine_depth = 1.0, fine_min = 0.25, min_thickness = 0.001,
                                   glaciers = false, lakes = false)
    P = ValaisRegion((r.E[1] - pad, r.E[2] + pad), (r.N[1] - pad, r.N[2] + pad), r.z)
    dxf = dx / sub
    xf, yf, topof = fetch_topography(P; dx = dxf)
    topof = replace(topof, NaN => 0.0)
    # Quaternary thickness of the region, none in the padding
    xr, yr, _, tr = fetch_swissbedrock(r; dx = dxf)
    tmud = zeros(length(xf), length(yf))
    oi, oj = round(Int, (xr[1] - xf[1]) / dxf), round(Int, (yr[1] - yf[1]) / dxf)
    tmud[oi .+ eachindex(xr), oj .+ eachindex(yr)] .= replace(tr, NaN => 0.0)
    If = model_interfaces(xf, yf, topof, tmud; ice_water_thickness(P, topof; dx = dxf, glaciers, lakes)...)
    maps = (; x = xf, y = yf, z = nothing, interfaces = If, morcles = read_morcles())
    # coarse columns: block means of the fine surfaces
    nx, ny = length(xf) ÷ sub, length(yf) ÷ sub
    blk(i) = (i - 1) * sub + 1:i * sub
    bmean(A) = [mean(@view A[blk(i), blk(j)]) for i in 1:nx, j in 1:ny]
    xc = [mean(xf[blk(i)]) for i in 1:nx]; yc = [mean(yf[blk(j)]) for j in 1:ny]
    Ic = (; topo = bmean(If.topo), base_cs = bmean(If.base_cs), conrad = bmean(If.conrad), moho_eu = bmean(If.moho_eu))
    S = meshing_surfaces(Ic, r.z[1]; fine_depth, fine_min)
    L = curvilinear_levels(xc, yc, S; dz, min_thickness)
    z = L.z; nzl = size(z, 3)
    rho = zeros(nx, ny, nzl)
    Threads.@threads for j in 1:ny
        for i in 1:nx
            t = Ic.topo[i, j]
            for k in 1:nzl
                # depth range of the dual cell below the topography of the coarse column
                dtop = t - (k < nzl ? (z[i, j, k] + z[i, j, k + 1]) / 2 : z[i, j, k])
                dbot = t - (k > 1 ? (z[i, j, k] + z[i, j, k - 1]) / 2 : z[i, j, k])
                s = 0.0
                for jf in blk(j), i_f in blk(i), m in 1:nsub_z
                    d = max(dtop + (m - 0.5) / nsub_z * (dbot - dtop), 0.0)
                    s += _model_point(maps, i_f, jf, If.topo[i_f, jf] - d)[4]
                end
                rho[i, j, k] = s / (sub^2 * nsub_z)
            end
        end
    end
    @info "Curvilinear density mesh: $(nx) × $(ny) × $(nzl) = $(round(length(z) / 1e6, digits = 2)) M nodes; cells per layer $(L.nz)"
    return (; x = xc, y = yc, z, rho, levels = L.levels, nz = L.nz, topo = Ic.topo, maps)
end

"""
    curvilinear_bouguer_contrast(z, rho; rho_topo = 2670.0, dzbin = 0.5) -> (; contrast, zref, rhoref)

Density contrast of the nodes of a curvilinear mesh (heights `z`) whose gravity is the complete Bouguer anomaly:
`rho - rho_topo` above sea level (the masses that the Bouguer correction replaced by `rho_topo`), and below sea level
`rho` minus a laterally uniform reference `rhoref(zref)`, the mean density of the nodes in bins of `dzbin` km
(interpolated linearly to the height of every node). A laterally uniform reference only adds a smooth regional field.
"""
function curvilinear_bouguer_contrast(z, rho; rho_topo = RHO_BOUGUER, dzbin = 0.5)
    lo = floor(minimum(z) / dzbin) * dzbin
    zref = collect(lo + dzbin / 2:dzbin:0.0)
    s = zeros(length(zref)); n = zeros(length(zref))
    for c in eachindex(z)
        z[c] > 0 && continue
        b = clamp(floor(Int, (z[c] - lo) / dzbin) + 1, 1, length(zref))
        s[b] += rho[c]; n[b] += 1
    end
    rhoref = s ./ max.(n, 1)
    # empty bins (none in practice): nearest filled bin
    filled = findall(>(0), n)
    for b in findall(==(0), n)
        rhoref[b] = rhoref[filled[argmin(abs.(filled .- b))]]
    end
    function ref(zq)
        zq <= zref[1] && return rhoref[1]
        zq >= zref[end] && return rhoref[end]
        b = searchsortedlast(zref, zq); w = (zq - zref[b]) / (zref[b + 1] - zref[b])
        return (1 - w) * rhoref[b] + w * rhoref[b + 1]
    end
    contrast = [z[c] > 0 ? rho[c] - rho_topo : rho[c] - ref(z[c]) for c in CartesianIndices(z)]
    return (; contrast, zref, rhoref)
end

"""
    interpolate_curvilinear(xg, yg, Zg, F, x, y, Z) -> array of size(Z)

The field `F` at the nodes of a curvilinear mesh (columns `xg`, `yg`, heights `Zg`, top level = topography) at the nodes of
another one (columns `x`, `y`, heights `Z`, top level = topography), e.g. the density change of
[`curvilinear_density_model`](@ref) on the mesh of a wave simulation. Bilinear between the four nearest columns
(clamped at the edges), linear along each of them, at the depth below the topography of the target node: the
surface of both meshes, and the material under it, stay aligned. Below the bottom or above the top of a column `F` is
taken constant.
"""
function interpolate_curvilinear(xg, yg, Zg, F, x, y, Z)
    ix(v, g) = (q = clamp(searchsortedlast(g, v), 1, length(g) - 1); (q, clamp((v - g[q]) / (g[q + 1] - g[q]), 0, 1)))
    function colinterp(i, j, zq)
        c = view(Zg, i, j, :)
        zq <= c[1] && return F[i, j, 1]
        zq >= c[end] && return F[i, j, end]
        k = clamp(searchsortedlast(c, zq), 1, length(c) - 1)
        w = (zq - c[k]) / (c[k + 1] - c[k])
        return (1 - w) * F[i, j, k] + w * F[i, j, k + 1]
    end
    out = zeros(eltype(F), size(Z))
    Threads.@threads for j in eachindex(y)
        q, t = ix(y[j], yg)
        for i in eachindex(x)
            p, s = ix(x[i], xg)
            for k in axes(Z, 3)
                d = Z[i, j, end] - Z[i, j, k]                      # depth below the topography
                out[i, j, k] = (1 - s) * (1 - t) * colinterp(p, q, Zg[p, q, end] - d) + s * (1 - t) * colinterp(p + 1, q, Zg[p + 1, q, end] - d) +
                               (1 - s) * t * colinterp(p, q + 1, Zg[p, q + 1, end] - d) + s * t * colinterp(p + 1, q + 1, Zg[p + 1, q + 1, end] - d)
            end
        end
    end
    return out
end
