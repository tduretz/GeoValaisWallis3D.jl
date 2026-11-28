"""
Phase numbering of the 3D model.
"""
const PHASES = (
    air = 0,
    quaternary = 1,              # unconsolidated deposits (swissBEDROCK TMUD)
    consolidated_sediments = 2,  # Alpine sedimentary cover/nappes (3D-ALPS; Morcles sedimentary units)
    upper_crust = 3,
    lower_crust = 4,
    mantle = 5,                  # below the European Moho (Spada et al. 2013)
    adriatic_mantle = 6,         # Adriatic mantle wedge (Ivrea body): below Adriatic Moho, above European Moho
    ice = 7,                     # glaciers (Swiss Glacier Thickness R2020; `glaciers = true`)
    water = 8,                   # lakes (swissBATHY3D; `lakes = true`)
)
const PHASE_NAMES = Dict(v => k for (k, v) in pairs(PHASES))

"""
    model_interfaces(x, y, topo, tmud) -> NamedTuple of matrices (km / kg m⁻³)

Layer interfaces on the LV95 grid `x`, `y` from the combined data sets:
- `topo`    : surface (swissALTIRegio): the surface of the glaciers and lakes
- `ice`, `water` : ice thickness and water depth (km; keywords, default 0, see [`ice_water_thickness`](@ref))
- `solid`   : solid surface = topo - ice - water (the glacier bed, the lake floor, or the topography)
- `bed`     : bedrock = solid - Quaternary thickness (swissBEDROCK counts it from the glacier bed and the lake floor;
              0 thickness outside CH)
- `base_cs` : base of consolidated sediments = min(bedrock, 3D-ALPS top of upper crust)
- `moho_eu`, `moho_ad` : European / Adriatic Moho (Spada et al. 2013)
- `conrad`  : upper/lower crust boundary, placed at the 3D-ALPS fraction of the crystalline crust
             thickness between `base_cs` and the Spada Moho
- `rho_*`   : 3D-ALPS densities of consolidated sediments, upper/lower crust, mantle
"""
function model_interfaces(x, y, topo, tmud; ice = nothing, water = nothing)
    Ee, Ne, Ze = read_moho_spada(1)
    Ea, Na, Za = read_moho_spada(2)
    moho_eu = interpolate_scattered(Ee, Ne, Ze, x, y; radius = 7.5)
    moho_ad = interpolate_scattered(Ea, Na, Za, x, y; radius = 7.5, maxdist = 5.0)
    a = alps3d_on_grid(x, y)
    t = replace(topo, NaN => 0.0)
    q = clamp.(replace(tmud, NaN => 0.0), 0.0, Inf)
    h = ice === nothing ? zeros(size(t)) : Float64.(ice)
    w = water === nothing ? zeros(size(t)) : Float64.(water)
    solid = t .- h .- w
    bed = solid .- q
    base_cs = max.(min.(bed, a.top_uc), moho_eu .+ 2.0)
    f = clamp.((a.top_uc .- a.top_lc) ./ (a.top_uc .- a.top_lm), 0.05, 0.95)
    conrad = base_cs .- f .* (base_cs .- moho_eu)
    return (; topo = t, ice = h, water = w, solid, bed, quat_thickness = q, base_cs, conrad, moho_eu, moho_ad,
            rho_cs = a.rho_cs, rho_uc = a.rho_uc, rho_lc = a.rho_lc, rho_lm = a.rho_lm)
end

"""
    top_cell_height(z, t) -> Float64

Height of the top solid cell of a column with surface elevation `t` on the vertical levels `z` (increasing): the cell
of the highest level `z[k] ≤ t`, from the midpoint to the level below to the midpoint to the level above.
"""
function top_cell_height(z, t)
    n = length(z)
    k = clamp(searchsortedlast(z, t), 1, n)
    lo = k > 1 ? (z[k] - z[k-1]) / 2 : 0.0
    hi = k < n ? (z[k+1] - z[k]) / 2 : 0.0
    return k == 1 ? 2hi : k == n ? 2lo : lo + hi
end

"""
    drop_thin_deposits(topo, tmud, z, fraction) -> tmud

The Quaternary thickness `tmud` with the deposits dropped (thickness 0, the bedrock reaches the surface) in the columns
where they fill less than `fraction` of the top solid cell ([`top_cell_height`](@ref) on the levels `z`). A cell takes
the material at its level, so a deposit of a few metres would otherwise make a whole top cell of Quaternary when it
reaches below that level, and none when it does not. `NaN` (no data) stays `NaN`; `topo` is in km like `tmud` and `z`.
"""
function drop_thin_deposits(topo, tmud, z, fraction)
    out = copy(tmud)
    for c in eachindex(topo, tmud)
        q = tmud[c]
        (isnan(q) || q <= 0) && continue
        t = isnan(topo[c]) ? 0.0 : topo[c]
        q < fraction * top_cell_height(z, t) && (out[c] = zero(q))
    end
    return out
end

"""
    valais_model_maps(r = valais_region(); dx = 0.25, z = nothing, fine_band = :auto, dz_fine = 0.025,
                      dz_max_up = 0.25, dz_max_down = 2.0, growth = 1.08, morcles = true, min_fill = 0.0,
                      glaciers = false, lakes = false) -> NamedTuple

The 3D model as maps: the grid axes `x`, `y` (LV95, km; uniform spacing `dx`) and `z` (vertical levels, see
[`refined_levels`](@ref); `fine_band = :auto` derives the refined band from swissBEDROCK, or pass a tuple
`(zlow, zhigh)`; pass `z` to give the levels explicitly), the layer `interfaces` on `x`, `y`
([`model_interfaces`](@ref)) and the Morcles voxel model (`morcles`, or `nothing`).

`min_fill`: Quaternary deposits that fill less than this fraction of the column's top cell are dropped
([`drop_thin_deposits`](@ref)); e.g. `min_fill = 0.5` makes a top cell Quaternary only if the deposit fills at least
half of it. The default `0` keeps all deposits.

`glaciers = true`, `lakes = true`: the glaciers (ice, `PHASES.ice`) and lakes (water, `PHASES.water`) between the
topography and the solid surface, see [`ice_water_thickness`](@ref). Without them the model fills the glaciers with rock
and the lakes with Quaternary deposits.

This is all that is needed to fill any part of the 3D grid with [`valais_model_block`](@ref); it is small (2D), so every
process of a parallel run can hold it.
"""
function valais_model_maps(r::ValaisRegion = valais_region(); dx = 0.25, z = nothing, fine_band = :auto, dz_fine = 0.025,
                           dz_max_up = 0.25, dz_max_down = 2.0, growth = 1.08, morcles = true, min_fill = 0.0,
                           glaciers = false, lakes = false)
    x, y, topo = fetch_topography(r; dx = dx)
    _, _, _, tmud = fetch_swissbedrock(r; dx = dx)
    if z === nothing
        band = fine_band === :auto ? fine_band_from_bedrock(topo, tmud) : fine_band
        z = refined_levels(r.z[1], r.z[2]; fine = band, dz_fine, dz_max_up, dz_max_down, growth)
        @info "Vertical grid: $(length(z)) levels, dz = $(dz_fine) km in [$(round(band[1], digits=3)), $(round(band[2], digits=3))] km"
    end
    if min_fill > 0
        thin = drop_thin_deposits(topo, tmud, z, min_fill)
        @info "Quaternary deposits filling less than $(min_fill) of the top cell dropped in $(count(.!isequal.(thin, tmud))) of $(count(>(0), replace(tmud, NaN => 0.0))) columns with deposits"
        tmud = thin
    end
    B = ice_water_thickness(r, topo; dx, glaciers, lakes)
    interfaces = model_interfaces(x, y, topo, tmud; B...)
    return (; x, y, z, interfaces, morcles = morcles ? read_morcles() : nothing)
end

"""
    _model_point(maps, i, j, zk, V = DEFAULT_VELOCITIES) -> (phase, Vp, Vs, Rho, MorclesUnit)

Material of the model at the height `zk` (km) in the column `maps.x[i]`, `maps.y[j]` of [`valais_model_maps`](@ref) (or
[`valais_curvilinear_maps`](@ref)): the rules of [`build_valais_model`](@ref), shared by the rectilinear and the
curvilinear grids. Above the topography, the phase is `PHASES.air` and the other values are 0. Between the topography and
the solid surface (`interfaces.solid`), the phase is `PHASES.water` (lakes) or `PHASES.ice` (glaciers). The velocity laws
of the solid phases take the depth below the solid surface.
"""
# velocity law of a phase: from `V`, or the default for ice and water when `V` (user laws) has none
_law(V, name) = haskey(V, name) ? V[name] : DEFAULT_VELOCITIES[name]

function _model_point(maps, i, j, zk, V = DEFAULT_VELOCITIES)
    I, mm = maps.interfaces, maps.morcles
    t = I.topo[i, j]
    zk > t && return (Int8(PHASES.air), 0f0, 0f0, 0f0, Int8(0))
    sf = haskey(I, :solid) ? I.solid[i, j] : t
    if zk > sf                                 # lake water or glacier ice
        ph, law = I.water[i, j] > 0 ? (PHASES.water, _law(V, :water)) : (PHASES.ice, _law(V, :ice))
        vp, vs = velocity(law, t - zk)
        return (Int8(ph), Float32(vp), Float32(vs), Float32(law.rho), Int8(0))
    end
    d = sf - zk                                # depth below the solid surface
    mu = mm === nothing ? Int8(0) : morcles_unit_at(mm, maps.x[i], maps.y[j], zk)
    mad = I.moho_ad[i, j]; meu = I.moho_eu[i, j]
    # a node exactly on the base of the consolidated sediments, the Conrad or the Moho: the material below it, or with
    # `interface_above` (curvilinear meshes for SeismicWavePC) the one above
    ab = get(maps, :interface_above, false)
    above(a, b) = ab ? a >= b : a > b
    ph, law, rho = if zk > I.bed[i, j]
        PHASES.quaternary, V.quaternary, NaN
    elseif mu > 0
        mu == 1 ? (PHASES.upper_crust, V.upper_crust, I.rho_uc[i, j]) :
                  (PHASES.consolidated_sediments, V.consolidated_sediments, I.rho_cs[i, j])
    elseif above(zk, I.base_cs[i, j])
        PHASES.consolidated_sediments, V.consolidated_sediments, I.rho_cs[i, j]
    elseif !isnan(mad) && mad > meu && meu < zk <= mad
        PHASES.adriatic_mantle, V.adriatic_mantle, I.rho_lm[i, j]
    elseif above(zk, I.conrad[i, j])
        PHASES.upper_crust, V.upper_crust, I.rho_uc[i, j]
    elseif above(zk, meu)
        PHASES.lower_crust, V.lower_crust, I.rho_lc[i, j]
    else
        PHASES.mantle, V.mantle, I.rho_lm[i, j]
    end
    vp, vs = velocity(law, d)
    return (Int8(ph), Float32(vp), Float32(vs), Float32(isnan(law.rho) ? rho : law.rho), zk <= I.bed[i, j] ? mu : Int8(0))
end

"""
    valais_model_block(maps, ix = eachindex(maps.x), iy = eachindex(maps.y); velocities = DEFAULT_VELOCITIES,
                       temperature = false, air_value = NaN) -> NamedTuple

Fill the columns `maps.x[ix]`, `maps.y[iy]` (index ranges) of the 3D model of [`valais_model_maps`](@ref) over all
levels `maps.z`: arrays `Phases`, `Vp`, `Vs`, `Rho`, `MorclesUnit` (and `Temperature`) of size
`(length(ix), length(iy), length(maps.z))`, as in [`build_valais_model`](@ref). A block is identical to that part of
the whole model, so a model too large for one process can be built block by block, e.g. one per MPI rank.
"""
function valais_model_block(maps, ix = eachindex(maps.x), iy = eachindex(maps.y); velocities = DEFAULT_VELOCITIES,
                            temperature = false, air_value = NaN)
    # rectilinear levels (vector), the heights of a curvilinear mesh (global array), or none: then the heights of these
    # columns are made here (curvilinear_heights; local indices)
    loc = maps.z === nothing
    z = loc ? curvilinear_heights(maps, ix, iy) : maps.z
    nx, ny, nz = length(ix), length(iy), size(z, ndims(z))
    Phases = zeros(Int8, nx, ny, nz)
    Vp = zeros(Float32, nx, ny, nz); Vs = zeros(Float32, nx, ny, nz); Rho = zeros(Float32, nx, ny, nz)
    Morc = zeros(Int8, nx, ny, nz)
    Threads.@threads for jl in 1:ny
        j = iy[jl]
        for il in 1:nx, k in 1:nz
            i = ix[il]
            zk = z isa AbstractVector ? z[k] : loc ? z[il, jl, k] : z[i, j, k]
            ph, vp, vs, rho, mu = _model_point(maps, i, j, zk, velocities)
            ph == PHASES.air && continue                       # air (set to air_value below)
            Phases[il, jl, k] = ph
            Vp[il, jl, k] = vp; Vs[il, jl, k] = vs; Rho[il, jl, k] = rho
            Morc[il, jl, k] = mu
        end
    end
    air = Phases .== PHASES.air
    for A in (Vp, Vs, Rho)
        A[air] .= Float32(air_value)
    end
    fields = (Phases = Phases, Vp = Vp, Vs = Vs, Rho = Rho, MorclesUnit = Morc)
    if temperature
        T = Float32.(temperature_3dalps(maps.x[ix], maps.y[iy], z isa AbstractVector || loc ? z : z[ix, iy, :]))
        T[air] .= Float32(air_value)
        fields = merge(fields, (Temperature = T,))
    end
    return fields
end

"""
    build_valais_model(r = valais_region(); dx = 0.25, dx_surface = 0.1, z = nothing,
                       fine_band = :auto, dz_fine = 0.025, dz_max_up = 0.25, dz_max_down = 2.0,
                       growth = 1.08, velocities = DEFAULT_VELOCITIES, morcles = true, min_fill = 0.0,
                       temperature = true, air_value = NaN, glaciers = false, lakes = false,
                       outdir = datadir("output"), write = true)

Build the Valais 3D geological / seismic velocity model on a rectilinear Cartesian LV95 grid (km):
uniform horizontal spacing `dx`, vertical levels refined (`dz_fine`) in the elevation band of the
Quaternary valley fill (see [`refined_levels`](@ref); `fine_band = :auto` derives it from
swissBEDROCK, or pass a tuple `(zlow, zhigh)`; pass `z` to give the levels explicitly).

Data: swissALTIRegio topography, swissBEDROCK Quaternary thickness, swisstopo tectonic map (surface),
3D-ALPS layers & densities, Spada et al. (2013) Moho, 3D-ALPS-TR temperature, Morcles nappe model.

Returns `(; surface, model, z, tecto)`. `model` is a `CartData` with fields `Phases` (see
[`PHASES`](@ref)), `Vp`, `Vs` (km/s), `Rho` (kg/m³), `Temperature` (°C) and `MorclesUnit`
(index into [`MORCLES_UNITS`](@ref), 0 outside). In air cells (`Phases == 0`) Vp, Vs, Rho and
Temperature are set to `air_value` (default `NaN`, which ParaView can render transparent via the
color map's "NaN Opacity"; use e.g. `air_value = 0` for solvers that do not accept NaN).

The 3D model is [`valais_model_maps`](@ref) filled by [`valais_model_block`](@ref) over the whole grid; use these two
to build only a part of it.
"""
function build_valais_model(r::ValaisRegion = valais_region(); dx = 0.25, dx_surface = 0.1, z = nothing,
                            fine_band = :auto, dz_fine = 0.025, dz_max_up = 0.25, dz_max_down = 2.0,
                            growth = 1.08, velocities = DEFAULT_VELOCITIES, morcles = true, min_fill = 0.0,
                            temperature = true, air_value = NaN, glaciers = false, lakes = false,
                            outdir = datadir("output"), write = true)
    # ---------------- surface (high resolution) ----------------
    xs, ys, topo_s = fetch_topography(r; dx = dx_surface)
    _, _, bem_s, tmud_s = fetch_swissbedrock(r; dx = dx_surface)
    tecto = tectonic_map_on_grid(xs, ys)
    B_s = ice_water_thickness(r, topo_s; dx = dx_surface, glaciers, lakes)
    I_s = model_interfaces(xs, ys, topo_s, tmud_s; B_s...)
    Xs, Ys, Zs = xyz_grid(xs, ys, 0.0)
    Zs[:, :, 1] .= topo_s
    s3(a) = reshape(a, size(a)..., 1)
    surface = CartData(Xs, Ys, Zs, (
        Topography = s3(topo_s), BedrockElevation = s3(bem_s), QuaternaryThickness = s3(tmud_s),
        TectoLevel1 = s3(tecto.level1), TectoLevel2 = s3(tecto.level2), TectoUnit = s3(tecto.unit),
        QuaternaryMap = s3(Int.(tecto.quaternary)), IceThickness = s3(I_s.ice), WaterDepth = s3(I_s.water),
        BaseConsolidatedSediments = s3(I_s.base_cs), UpperLowerCrust = s3(I_s.conrad),
        MohoEuropean = s3(I_s.moho_eu), MohoAdriatic = s3(I_s.moho_ad),
    ))

    # ---------------- 3D grid ----------------
    maps = valais_model_maps(r; dx, z, fine_band, dz_fine, dz_max_up, dz_max_down, growth, morcles, min_fill, glaciers, lakes)
    z = maps.z
    nx, ny, nz = length(maps.x), length(maps.y), length(z)
    @info "3D grid: $(nx) × $(ny) × $(nz) = $(round(nx*ny*nz/1e6, digits=1)) M cells"
    fields = valais_model_block(maps; velocities, temperature, air_value)
    X, Y, Zg = xyz_grid(maps.x, maps.y, z)
    model = CartData(X, Y, Zg, fields)

    if write
        mkpath(outdir)
        write_paraview(surface, joinpath(outdir, "valais_surface"))
        write_paraview(model, joinpath(outdir, "valais_model"))
        _write_legend(joinpath(outdir, "legend.csv"), tecto)
        writedlm(joinpath(outdir, "z_levels_km.txt"), z)
    end
    return (; surface, model, z, tecto)
end

function _write_legend(file, t::TectonicMap)
    open(file, "w") do io
        println(io, "field,code,name")
        for (field, names) in (("TectoLevel1", t.names1), ("TectoLevel2", t.names2), ("TectoUnit", t.names_unit))
            for (i, n) in enumerate(names)
                println(io, field, ",", i, ",\"", n, "\"")
            end
        end
        for (n, p) in pairs(PHASES)
            println(io, "Phases,", p, ",", n)
        end
        for (i, n) in enumerate(MORCLES_UNITS)
            println(io, "MorclesUnit,", i, ",", n)
        end
    end
end

"""
    cartdata_to_geodata(d::CartData) -> GeoData

Convert a `CartData` in LV95 km coordinates to a `GeoData` (lon, lat, depth in km), e.g. to
combine it with seismic tomography or other GeophysicalModelGenerator datasets.
"""
function cartdata_to_geodata(d::CartData)
    lon, lat = lv95_to_lonlat(NumValue(d.x), NumValue(d.y))
    return GeoData(lon, lat, NumValue(d.z), d.fields)
end
