# Figures of the model on the rectilinear grid, and of the data sets combined in it.
#     julia --project=figures -t auto figures/rectilinear/make_figures.jl [figure ...]
# Figures: overview layers sections subsoil gravity density morcles (default: all). Output: figures/rectilinear/*.jpg (needs ImageMagick).
# density runs the density inversion of the Bouguer anomaly (invert_density: GMG voxel_grav / voxel_grav_adjoint + LSQR);
# GVW_BACKEND=metal runs it on the Apple GPU (Float32), otherwise on the CPU (Float64, slow). Its result is cached in
# data/output/density_inversion_figure.jls (delete it to rerun).
# Uses the cached data of data/raw/ (downloads what is missing, as the examples do).
using GeoValaisWallis3D, CairoMakie, Statistics, Serialization

const FIGDIR = @__DIR__
mkpath(FIGDIR)
const WHICH = isempty(ARGS) ? ["overview", "layers", "sections", "subsoil", "gravity", "density", "morcles"] : ARGS
const R = valais_region()
CairoMakie.activate!(type = "png", px_per_unit = 1.25)
set_theme!(fontsize = 13, Axis = (titlesize = 14, titlefont = :regular, xgridvisible = false, ygridvisible = false),
           Axis3 = (titlefont = :regular,))

# places for orientation (LV95, km)
const PLACES = [("Monthey", 2562.0, 1122.5), ("Martigny", 2571.5, 1105.5), ("Sion", 2594.0, 1120.0),
                ("Sierre", 2607.0, 1127.0), ("Visp", 2634.5, 1127.0), ("Brig", 2642.0, 1129.5),
                ("Zermatt", 2624.5, 1097.0), ("Leukerbad", 2614.5, 1136.5)]
# cross-sections of the "sections" figure
const SEC_A = (name = "A", E = 2594.0)          # N–S through Sion
const SEC_B = (name = "B", N = 1120.0)          # W–E along the Rhône valley at Sion

function hillshade(z, dx; az = 315.0, alt = 45.0, exag = 1.0)
    z = replace(z, NaN => 0.0) .* exag
    gx = similar(z); gy = similar(z)
    nx, ny = size(z)
    for j in 1:ny, i in 1:nx
        gx[i, j] = (z[min(i + 1, nx), j] - z[max(i - 1, 1), j]) / ((min(i + 1, nx) - max(i - 1, 1)) * dx)
        gy[i, j] = (z[i, min(j + 1, ny)] - z[i, max(j - 1, 1)]) / ((min(j + 1, ny) - max(j - 1, 1)) * dx)
    end
    a, e = deg2rad(az), deg2rad(alt)
    lx, ly, lz = cos(e) * sin(a), cos(e) * cos(a), sin(e)
    @. (-gx * lx - gy * ly + lz) / sqrt(gx^2 + gy^2 + 1)
end
shade!(ax, x, y, topo; kw...) = heatmap!(ax, x, y, hillshade(topo, x[2] - x[1]; kw...), colormap = :grays, colorrange = (0.2, 1.05))
function places!(ax; color = :black, names = true)
    for (n, e, nn) in PLACES
        scatter!(ax, [e], [nn], color = color, markersize = 6, strokecolor = :white, strokewidth = 1)
        names && text!(ax, e, nn, text = n, offset = (4, 3), fontsize = 12, color = color, glowcolor = (:white, 0.8), glowwidth = 2)
    end
end
mapaxis(f; kw...) = Axis(f; aspect = DataAspect(), xlabel = "E (km, LV95)", ylabel = "N (km)", kw...)
limits_region!(ax) = limits!(ax, R.E..., R.N...)
# JPEG (quality 85): the shaded relief makes PNGs of 3–8 MB; as JPEG they stay below ~1 MB for the README
function savefig(name, fig)
    png = joinpath(tempdir(), name * ".png"); save(png, fig)
    jpg = joinpath(FIGDIR, name * ".jpg")
    run(`magick $png -quality 85 $jpg`)            # ImageMagick
    println("written figures/rectilinear/$name.jpg ($(round(filesize(jpg) / 1e6, digits = 2)) MB)")
end
nanmask(a, m) = (b = Float64.(a); b[m] .= NaN; b)

# --------------------------------------------------------------------------------------------- overview
if "overview" in WHICH
    dx = 0.1
    x, y, topo = fetch_topography(R; dx)
    t = tectonic_map_on_grid(x, y)
    m = read_morcles()
    fig = Figure(size = (1500, 820))
    ax = mapaxis(fig[1, 1], title = "Topography (swissALTIRegio) and tectonic units (swisstopo 1:500 000, level 2)")
    shade!(ax, x, y, topo)
    used = sort(unique(t.level2[t.level2 .> 0]))
    cols = Makie.to_colormap(:tab20)
    code = zeros(size(t.level2)) .* NaN
    for (k, c) in enumerate(used)
        code[t.level2 .== c] .= k
    end
    heatmap!(ax, x, y, code, colormap = cgrad(cols[1:length(used)], length(used), categorical = true),
             colorrange = (0.5, length(used) + 0.5), alpha = 0.55, nan_color = :transparent)
    contour!(ax, x, y, topo, levels = [0.372], color = (:navy, 0.8), linewidth = 0.8)       # Lake Geneva shore
    mx, my = extrema(m.x), extrema(m.y)
    lines!(ax, [mx[1], mx[2], mx[2], mx[1], mx[1]], [my[1], my[1], my[2], my[2], my[1]], color = :black, linewidth = 2)
    text!(ax, mx[2], my[2], text = "Morcles 3D model", align = (:left, :top), fontsize = 11, offset = (3, -3), glowcolor = (:white, 0.8), glowwidth = 2)
    lines!(ax, [SEC_A.E, SEC_A.E], collect(R.N), color = :black, linestyle = :dash)
    lines!(ax, collect(R.E), [SEC_B.N, SEC_B.N], color = :black, linestyle = :dash)
    text!(ax, SEC_A.E, R.N[2], text = "A'", align = (:left, :top), offset = (4, -4)); text!(ax, SEC_A.E, R.N[1], text = "A", align = (:left, :bottom), offset = (4, 4))
    text!(ax, R.E[1], SEC_B.N, text = "B", align = (:left, :bottom), offset = (4, 4)); text!(ax, R.E[2], SEC_B.N, text = "B'", align = (:right, :bottom), offset = (-4, 4))
    places!(ax)
    limits_region!(ax)
    Legend(fig[1, 2], [PolyElement(color = (cols[k], 0.75)) for k in eachindex(used)], t.names2[used], "tectonic unit",
           framevisible = false, labelsize = 11, rowgap = 0)
    colsize!(fig.layout, 1, Relative(0.78))
    savefig("overview", fig)
end

# ------------------------------------------------------------------------------------- interfaces and layers
if "layers" in WHICH
    dx = 0.25
    x, y, topo = fetch_topography(R; dx)
    _, _, _, tmud = fetch_swissbedrock(R; dx)
    I = model_interfaces(x, y, topo, tmud)
    fig = Figure(size = (1500, 760))
    panels = [
        ("Topography (swissALTIRegio)", topo, :batlowW, (0.3, 4.2), "elevation (km)"),
        ("Quaternary thickness (swissBEDROCK)", nanmask(1e3 .* I.quat_thickness, I.quat_thickness .<= 0.002), Reverse(:lapaz), (0, 600), "thickness (m)"),
        ("Base of consolidated sediments (3D-ALPS)", I.base_cs, :batlow, (-12, 3), "elevation (km)"),
        ("Upper / lower crust boundary (3D-ALPS, Spada Moho)", I.conrad, :batlow, (-35, -10), "elevation (km)"),
        ("European Moho (Spada et al. 2013)", I.moho_eu, :batlow, (-60, -28), "elevation (km)"),
        ("Adriatic Moho (Spada et al. 2013)", I.moho_ad, :batlow, (-60, -10), "elevation (km)"),
    ]
    for (k, (title, A, cm, cr, lab)) in enumerate(panels)
        r, c = (k - 1) ÷ 3 + 1, (k - 1) % 3 + 1
        g = fig[r, c] = GridLayout()
        ax = mapaxis(g[1, 1]; title)
        shade!(ax, x, y, topo)
        hm = heatmap!(ax, x, y, A, colormap = cm, colorrange = cr, alpha = k == 1 ? 0.55 : 0.8, nan_color = :transparent)
        k >= 3 && contour!(ax, x, y, A, levels = range(cr..., step = k >= 5 ? 5 : 2.5), color = (:black, 0.5), linewidth = 0.6)
        places!(ax; names = false)
        limits_region!(ax)
        Colorbar(g[1, 2], hm, label = lab, height = Relative(0.7))
        r == 1 || hidexdecorations!(ax, label = false, ticklabels = false, ticks = false)
    end
    savefig("layers", fig)
end

# ------------------------------------------------------------------------------------------ cross-sections
if "sections" in WHICH
    maps = valais_model_maps(R; dx = 0.25, dz_fine = 0.025)
    ia = argmin(abs.(maps.x .- SEC_A.E)); jb = argmin(abs.(maps.y .- SEC_B.N))
    A = valais_model_block(maps, ia:ia, eachindex(maps.y); temperature = true)
    B = valais_model_block(maps, eachindex(maps.x), jb:jb; temperature = true)
    sq(a) = dropdims(a, dims = findfirst(==(1), size(a)))
    phcols = [:white, "#e9c46a", "#8fb8de", "#c0504d", "#7f6a93", "#3d5a40", "#2a9d8f"]
    phnames = ["air", "Quaternary", "consolidated sediments", "upper crust", "lower crust", "mantle", "Adriatic mantle"]
    phcm = cgrad(phcols, 7, categorical = true)
    fig = Figure(size = (1500, 1250))
    secs = [(A, maps.y, "A–A' (N–S, E = $(SEC_A.E) km)", "N (km)"), (B, maps.x, "B–B' (W–E, N = $(SEC_B.N) km)", "E (km)")]
    for (c, (S, h, name, xl)) in enumerate(secs)
        ph = Float64.(sq(S.Phases)); air = ph .== 0
        ax1 = Axis(fig[1, c], title = "$name: phases", ylabel = "elevation (km)")
        heatmap!(ax1, h, maps.z, ph, colormap = phcm, colorrange = (-0.5, 6.5))
        ax2 = Axis(fig[2, c], title = "$name: Vs (first-order velocity laws)", ylabel = "elevation (km)")
        h2 = heatmap!(ax2, h, maps.z, nanmask(sq(S.Vs), air), colormap = :lipari, colorrange = (0.3, 4.6), nan_color = :white)
        ax3 = Axis(fig[3, c], title = "$name: temperature (3D-ALPS-TR)", ylabel = "elevation (km)")
        h3 = heatmap!(ax3, h, maps.z, nanmask(sq(S.Temperature), air), colormap = :lajolla, colorrange = (0, 1200), nan_color = :white)
        contour!(ax3, h, maps.z, nanmask(sq(S.Temperature), air), levels = 200:200:1200, color = (:black, 0.5), linewidth = 0.6, labels = true, labelsize = 9)
        ax4 = Axis(fig[4, c], title = "$name: phases near the surface (vertically stretched)", ylabel = "elevation (km)", xlabel = xl)
        heatmap!(ax4, h, maps.z, ph, colormap = phcm, colorrange = (-0.5, 6.5))
        for ax in (ax1, ax2, ax3)
            ylims!(ax, -70, 5); xlims!(ax, extrema(h)...)
        end
        xlims!(ax4, extrema(h)...); ylims!(ax4, -1.0, 4.6)
        c == 2 && (Colorbar(fig[2, 3], h2, label = "Vs (km/s)"); Colorbar(fig[3, 3], h3, label = "T (°C)"))
    end
    Legend(fig[1, 3], [PolyElement(color = phcols[k], strokecolor = :grey60, strokewidth = 0.5) for k in 2:7], phnames[2:7], framevisible = false)
    Legend(fig[4, 3], [PolyElement(color = phcols[k], strokecolor = :grey60, strokewidth = 0.5) for k in 2:4], phnames[2:4], framevisible = false)
    savefig("sections", fig)
end

# ------------------------------------------------------------------------------ seismic subsoil classes (CREALP)
if "subsoil" in WHICH
    sm = subsoil_map(R; dx = 0.05)
    clcols = ["#3b5b92", "#5aa0c8", "#9ccf7a", "#f2c14e", "#e8743b", "#b5179e"]   # A (rock) … E (soft), F (special)
    clcm = cgrad(clcols, 6, categorical = true)
    cls = nanmask(sm.class, sm.class .== 0)
    ssm = Float64.(sm.ssm .> 0)
    fig = Figure(size = (1500, 1450))
    ax = mapaxis(fig[1, 1:2], title = "Seismic subsoil classes SIA 261 (BAFU national map; Rhône valley 1:25 000 maps of CREALP / canton of Valais), release $(sm.release)")
    shade!(ax, sm.x, sm.y, sm.topo)
    heatmap!(ax, sm.x, sm.y, cls, colormap = clcm, colorrange = (0.5, 6.5), nan_color = :transparent)
    contour!(ax, sm.x, sm.y, ssm, levels = [0.5], color = :black, linewidth = 1.2)
    places!(ax); limits_region!(ax)
    labels = ["A rock", "B", "C", "D", "E", "F special"]
    Legend(fig[1, 3], [[PolyElement(color = c) for c in clcols]; LineElement(color = :black)],
           [["$(SUBSOIL_CLASSES[k])" * (k == 1 ? " (rock)" : k == 5 ? " (thin soft layer on rock)" : k == 6 ? " (special)" : "") for k in 1:6];
            "spectral microzonation (classes replaced)"], "class", framevisible = false)
    # the subsoil map draped on the topography: Rhône valley from Martigny to Sierre
    ex, ny_ = (2566.0, 2612.0), (1100.0, 1132.0)
    ix = findall(e -> ex[1] <= e <= ex[2], sm.x); iy = findall(n -> ny_[1] <= n <= ny_[2], sm.y)
    xs, ys, zs = sm.x[ix], sm.y[iy], sm.topo[ix, iy]
    hs = clamp.(hillshade(zs, 0.05), 0, 1)
    rgb = map(CartesianIndices(zs)) do I
        g = 0.35 + 0.6hs[I]
        c = cls[ix[I[1]], iy[I[2]]]
        isnan(c) ? RGBf(g, g, g) : (k = Makie.to_color(clcols[Int(c)]); RGBf(k.r * (0.55 + 0.5g), k.g * (0.55 + 0.5g), k.b * (0.55 + 0.5g)))
    end
    ax3 = Axis3(fig[2, 1:3], title = "The same map draped on the topography: Rhône valley, Martigny–Sion–Sierre (vertical exaggeration 1.5, view from the south-west)",
                aspect = (ex[2] - ex[1], ny_[2] - ny_[1], 1.5 * 4.0), azimuth = -0.35π, elevation = 0.17π,
                xlabel = "E (km)", ylabel = "N (km)", zlabel = "z (km)", protrusions = 30)
    surface!(ax3, xs, ys, zs, color = rgb, shading = NoShading)
    for (n, e, nn) in PLACES
        (ex[1] <= e <= ex[2] && ny_[1] <= nn <= ny_[2]) || continue
        i = argmin(abs.(xs .- e)); j = argmin(abs.(ys .- nn))
        text!(ax3, Point3f(e, nn, zs[i, j] + 0.6), text = n, fontsize = 14, align = (:center, :bottom), color = :black, glowcolor = (:white, 0.85), glowwidth = 3)
    end
    zlims!(ax3, 0.3, 4.3)
    rowsize!(fig.layout, 1, Relative(0.5))
    savefig("subsoil", fig)
end

# --------------------------------------------------------------------------------------------- gravity
if "gravity" in WHICH
    gm = gravity_maps(R; dx = 0.25)
    fig = Figure(size = (1500, 900))
    panels = [
        ("Complete Bouguer anomaly, atlas 1:100 000 (points)", gm.bouguer, (-200, 20), Reverse(:RdBu), true),
        ("Complete Bouguer anomaly, map 1:500 000 (contours)", gm.bouguer500k, (-200, 20), Reverse(:RdBu), false),
        ("Attraction of the topography (2670 kg/m³)", gm.gtopo, (0, 300), :lajolla, false),
        ("Free-air anomaly = Bouguer + topography", gm.freeair, (-200, 200), Reverse(:RdBu), false),
    ]
    for (k, (title, A, cr, cm, pts)) in enumerate(panels)
        g = fig[(k - 1) ÷ 2 + 1, (k - 1) % 2 + 1] = GridLayout()
        local ax = mapaxis(g[1, 1]; title)
        shade!(ax, gm.x, gm.y, gm.topo)
        # diverging maps centred on zero (Bouguer: range not symmetric, so centre the colormap explicitly)
        crr = cm === :lajolla ? cr : (-1, 1) .* maximum(abs, cr)
        hm = heatmap!(ax, gm.x, gm.y, A, colormap = cm, colorrange = crr, alpha = 0.8, nan_color = :transparent)
        contour!(ax, gm.x, gm.y, A, levels = range(-300, 300, step = 20), color = (:black, 0.35), linewidth = 0.5)
        pts && scatter!(ax, gm.atlas.E, gm.atlas.N, markersize = 1.5, color = (:black, 0.5))
        k == 4 && scatter!(ax, gm.network.E, gm.network.N, marker = :utriangle, markersize = 8, color = :white, strokecolor = :black, strokewidth = 1)
        places!(ax; names = false); limits_region!(ax)
        Colorbar(g[1, 2], hm, label = "mGal", height = Relative(0.7))
    end
    Label(fig[0, 1:2], "swisstopo gravity: atlas points (dots), LSN2004 stations (triangles, free-air check: −3.6 ± 2.6 mGal)", fontsize = 14)
    savefig("gravity", fig)
end

# ------------------------------------------------------------------------- density inversion of the Bouguer anomaly
if "density" in WHICH
    const PAD, DXI, SIGMA = 20.0, 1.0, 2.0          # lateral padding and node spacing [km], data error [mGal]
    cache = datadir("output", "density_inversion_figure.jls")
    res = if isfile(cache)
        deserialize(cache)
    else
        if get(ENV, "GVW_BACKEND", "cpu") == "metal"
            @eval using Metal
            backend, ft = Base.invokelatest(() -> Metal.MetalBackend()), Float32
        else
            backend, ft = nothing, Float64
        end
        # starting model: the density of the 3D model (250 m), averaged onto the nodes of the inversion grid
        maps = valais_model_maps(R; dx = 0.25, dz_fine = 0.05, morcles = true)
        mb = valais_model_block(maps)
        air = mb.Phases .== 0
        padded = ValaisRegion((R.E[1] - PAD, R.E[2] + PAD), (R.N[1] - PAD, R.N[2] + PAD), R.z)
        xg, yg, topo = fetch_topography(padded; dx = DXI); topo = replace(topo, NaN => 0.0)
        zg = let zs = collect(-3.0:0.25:4.5), dz = 0.25     # 250 m down to −3 km, then growing by 25 % per level
            while zs[1] > R.z[1]
                dz = min(1.25dz, 4.0); pushfirst!(zs, max(zs[1] - dz, R.z[1]))
            end
            zs
        end
        g = density_on_grid(maps.x, maps.y, maps.z, mb.Rho, air, xg, yg, zg; topo)
        mb = nothing; GC.gc()
        c0 = bouguer_contrast(g.rho, g.rockfrac, zg)
        active = g.rockfrac .>= 0.5
        at = fetch_gravity_atlas(R; margin = PAD)
        sel = (xg[1] + 2 .<= at.E .<= xg[end] - 2) .& (yg[1] + 2 .<= at.N .<= yg[end] - 2) .& .!isnan.(at.bouguer)
        E, N, H, dobs = at.E[sel], at.N[sel], at.h[sel], at.bouguer[sel]
        inv = invert_density(xg, yg, zg, c0.contrast, active, topo, (E, N, H), dobs, SIGMA; z0 = 1.0, β = 2.0,
                             backend, float_type = ft)
        drho = zeros(size(inv.drho)); drho[active] .= inv.drho[active] ./ g.rockfrac[active]
        r = (; xg, yg, zg, topo, active, rho0 = g.rho, drho, E, N, dobs, dg0 = inv.dg0, dg = inv.dg, offset = inv.offset,
             iterations = inv.iterations, rms_history = inv.rms_history)
        serialize(cache, r); r
    end
    (; xg, yg, zg, topo, active, rho0, drho, E, N, dobs) = res
    inner = (R.E[1] .<= E .<= R.E[2]) .& (R.N[1] .<= N .<= R.N[2])
    rstart = dobs .- res.dg0 .- res.offset; rstart .-= mean(rstart[inner])     # starting model, same regional plane
    rinv = dobs .- res.dg
    rms(a) = sqrt(mean(abs2, a))
    box!(ax) = lines!(ax, [R.E[1], R.E[2], R.E[2], R.E[1], R.E[1]], [R.N[1], R.N[1], R.N[2], R.N[2], R.N[1]], color = :black, linestyle = :dash)
    fig = Figure(size = (1500, 1180))
    g1 = fig[1, 1] = GridLayout(); g2 = fig[2, 1] = GridLayout(); g3 = fig[3, 1] = GridLayout()
    dm = dobs .- mean(dobs[inner])
    for (n, (v, t, cr)) in enumerate([(dm, "observed Bouguer anomaly (atlas, mean removed)", 60),
                                      (rstart, "residual of the starting model: rms $(round(rms(rstart[inner]), digits = 1)) mGal", 40),
                                      (rinv, "residual of the inverted model: rms $(round(rms(rinv[inner]), digits = 1)) mGal", 8)])
        local ax = Axis(g1[1, 2n - 1], aspect = DataAspect(), title = t, titlesize = 13, xlabel = "E (km)", ylabel = n == 1 ? "N (km)" : "")
        sc = scatter!(ax, E, N, color = v, colormap = Reverse(:RdBu), colorrange = (-cr, cr), markersize = 3)
        contour!(ax, xg, yg, topo, levels = [1.0, 2.5], color = (:black, 0.3), linewidth = 0.5); box!(ax)
        Colorbar(g1[1, 2n], sc, label = "mGal", height = Relative(0.6))
    end
    local hm
    for (n, zz) in enumerate([0.0, -5.0, -15.0, -35.0])
        k = argmin(abs.(zg .- zz))
        local ax = Axis(g2[1, n], aspect = DataAspect(), title = "density change, z = $(round(zg[k], digits = 1)) km",
                        xlabel = "E (km)", ylabel = n == 1 ? "N (km)" : "")
        hm = heatmap!(ax, xg, yg, ifelse.(active[:, :, k], drho[:, :, k], NaN), colormap = Reverse(:RdBu), colorrange = (-200, 200), nan_color = :transparent)
        box!(ax); places!(ax; names = false)
    end
    Colorbar(g2[1, 5], hm, label = "Δρ (kg/m³)", height = Relative(0.6))
    i = argmin(abs.(xg .- SEC_A.E))
    ax = Axis(g3[1, 1], xlabel = "N (km)", ylabel = "elevation (km)", title = "A–A' (E = $(xg[i]) km): inverted density")
    h1 = heatmap!(ax, yg, zg, ifelse.(active[i, :, :], rho0[i, :, :] .+ drho[i, :, :], NaN), colormap = :batlow, colorrange = (2300, 3400), nan_color = :transparent)
    Colorbar(g3[1, 2], h1, label = "ρ (kg/m³)")
    ax = Axis(g3[1, 3], xlabel = "N (km)", title = "A–A': density change")
    h2 = heatmap!(ax, yg, zg, ifelse.(active[i, :, :], drho[i, :, :], NaN), colormap = Reverse(:RdBu), colorrange = (-200, 200), nan_color = :transparent)
    Colorbar(g3[1, 4], h2, label = "Δρ (kg/m³)")
    ax = Axis(g3[1, 5], xlabel = "LSQR iteration", ylabel = "rms misfit (mGal)", yscale = log10, yticks = [2, 5, 10, 20, 50, 100], title = "convergence (stops at σ = $(SIGMA) mGal)")
    lines!(ax, 0:length(res.rms_history) - 1, res.rms_history, color = "#2a78d6", linewidth = 2)
    hlines!(ax, [SIGMA], color = :gray, linestyle = :dash)
    colsize!(g3, 5, Relative(0.22))
    rowsize!(fig.layout, 1, Relative(0.3)); rowsize!(fig.layout, 2, Relative(0.27))
    Label(fig[0, 1], "Density inversion of the swisstopo Bouguer anomaly (GMG voxel_grav + voxel_grav_adjoint, LSQR, Li & Oldenburg depth weighting): " *
                     "$(length(E)) atlas points, $(count(active)) unknowns, $(res.iterations) iterations", fontsize = 14)
    savefig("density", fig)
end

# ------------------------------------------------------------------------------------------------ Morcles
if "morcles" in WHICH
    m = read_morcles()
    nx, ny, nz = size(m.unit)
    top = zeros(Int, nx, ny)                         # surface geology: the highest unit of each column
    for j in 1:ny, i in 1:nx
        k = findlast(!=(0), @view m.unit[i, j, :]); k === nothing || (top[i, j] = m.unit[i, j, k])
    end
    n = length(MORCLES_UNITS)
    cm = cgrad([Makie.to_color("#9e6b5a"); [get(cgrad(:roma), (k - 1) / (n - 2)) for k in 1:n-1]], n, categorical = true)
    fig = Figure(size = (1500, 700))
    ax = mapaxis(fig[1, 1], title = "Nappe de Morcles 3D model (Thornton et al. 2018): surface units")
    heatmap!(ax, m.x, m.y, nanmask(top, top .== 0), colormap = cm, colorrange = (0.5, n + 0.5), nan_color = :white)
    jc = argmin(abs.(m.y .- (m.y[1] + m.y[end]) / 2))
    lines!(ax, [m.x[1], m.x[end]], [m.y[jc], m.y[jc]], color = :black, linestyle = :dash)
    places!(ax)
    xlims!(ax, extrema(m.x)...); ylims!(ax, extrema(m.y)...)
    ax2 = Axis(fig[1, 2], title = "W–E section at N = $(round(m.y[jc], digits = 1)) km (dashed line)", xlabel = "E (km)", ylabel = "elevation (km)", aspect = DataAspect())
    heatmap!(ax2, m.x, m.z, nanmask(m.unit[:, jc, :], m.unit[:, jc, :] .== 0), colormap = cm, colorrange = (0.5, n + 0.5), nan_color = :white)
    xlims!(ax2, extrema(m.x)...)
    Legend(fig[1, 3], [PolyElement(color = cm[k]) for k in 1:n], MORCLES_UNITS, "formation (old → young)", framevisible = false, rowgap = 0, labelsize = 11)
    colsize!(fig.layout, 2, Relative(0.55))
    savefig("morcles", fig)
end

# ------------------------------------------------------------------------------- showcase (one 4:5 image, e.g. LinkedIn)
# Not in the default list: julia --project=figures -t auto figures/rectilinear/make_figures.jl showcase (needs a display: GLMakie; and the cache
# of the "density" figure). The 3D panels are rendered alone with GLMakie (depth buffer), trimmed with ImageMagick and placed
# as images; the rest is drawn with CairoMakie in the LaTeX font.
if "showcase" in WHICH
    cache = datadir("output", "density_inversion_figure.jls")
    isfile(cache) || error("run the density figure first (it caches the inversion in $cache)")
    dres = deserialize(cache)
    W = 1600                                          # logical width of the image (px_per_unit 2 -> 3200 × 4000 px)
    tex = (; regular = Makie.texfont(), bold = Makie.texfont(:bold), italic = Makie.texfont(:italic))
    ink, muted = RGBf(0.09, 0.09, 0.09), RGBf(0.33, 0.33, 0.31)
    shadecol(c, g) = (k = Makie.to_color(c); RGBf(k.r * (0.45 + 0.6g), k.g * (0.45 + 0.6g), k.b * (0.45 + 0.6g)))
    grey(g) = RGBf(0.3 + 0.65g, 0.3 + 0.65g, 0.3 + 0.65g)
    tmp = mktempdir()
    # render a 3D scene alone (GLMakie), trim the white border; returns the image file
    function render3d(name, sz, draw!)
        @eval using GLMakie
        Base.invokelatest() do
            GLMakie.activate!(px_per_unit = 2)
            f = Figure(size = sz, fonts = tex, figure_padding = 0)
            ax = Axis3(f[1, 1]; protrusions = 0, viewmode = :fit)
            hidedecorations!(ax); hidespines!(ax)
            draw!(ax)
            png = joinpath(tmp, name * ".png"); save(png, f)
            run(`magick $png -trim +repage $png`)
            CairoMakie.activate!(type = "png", px_per_unit = 2)
            png
        end
    end
    loadimg(f) = rotr90(Makie.FileIO.load(f))
    function imgaxis!(pos, f)
        img = loadimg(f)
        ax = Axis(pos, aspect = DataAspect(), backgroundcolor = :transparent)
        hidedecorations!(ax); hidespines!(ax)
        image!(ax, img, interpolate = true)
        ax, size(img)
    end

    # ---- 3D 1: subsoil classes draped on the Rhône valley
    sm = subsoil_map(R; dx = 0.05)
    clcols = ["#3b5b92", "#5aa0c8", "#9ccf7a", "#f2c14e", "#e8743b", "#b5179e"]
    ex, ny_ = (2566.0, 2612.0), (1100.0, 1132.0)
    ix = findall(e -> ex[1] <= e <= ex[2], sm.x); iy = findall(n -> ny_[1] <= n <= ny_[2], sm.y)
    zs = sm.topo[ix, iy]; hs = clamp.(hillshade(zs, 0.05), 0, 1)
    icecol, lakecol, quakecol = "#f2f9ff", "#0d2a5c", "#ffd23f"
    iceshade(hs) = (c = Makie.to_color(icecol); RGBf(c.r * (0.72 + 0.3hs), c.g * (0.72 + 0.3hs), c.b * (0.72 + 0.3hs)))
    BWh = ice_water_thickness(R, sm.topo; dx = 0.05)         # glaciers (Grab et al. 2021) and lakes (swissBATHY3D)
    rgb = [sm.class[ix[I[1]], iy[I[2]]] == 0 ? (BWh.ice[ix[I[1]], iy[I[2]]] > 0.01 ? iceshade(hs[I]) : grey(hs[I])) :
           shadecol(clcols[sm.class[ix[I[1]], iy[I[2]]]], hs[I]) for I in CartesianIndices(zs)]
    hero = render3d("hero", (W, 900), ax -> begin
        ax.aspect = (ex[2] - ex[1], ny_[2] - ny_[1], 1.5 * 4.0); ax.azimuth = -0.35π; ax.elevation = 0.2π
        surface!(ax, sm.x[ix], sm.y[iy], zs, color = rgb, shading = NoShading)
        for (n, e, nn) in PLACES
            (ex[1] <= e <= ex[2] && ny_[1] <= nn <= ny_[2]) || continue
            i = argmin(abs.(sm.x[ix] .- e)); j = argmin(abs.(sm.y[iy] .- nn))
            text!(ax, Point3f(e, nn, zs[i, j] + 0.7), text = n, fontsize = 30, align = (:center, :bottom), color = ink,
                  glowcolor = (:white, 0.9), glowwidth = 4)
        end
        zlims!(ax, 0.3, 4.3)
    end)

    # ---- 3D 2: fence diagram of the model
    dx = 0.25
    x, y, topo = fetch_topography(R; dx)
    t = tectonic_map_on_grid(x, y)
    used = sort(unique(t.level2[t.level2 .> 0])); tcols = Makie.to_colormap(:tab20)
    tcol = Dict(c => tcols[k] for (k, c) in enumerate(used))
    hs2 = clamp.(hillshade(topo, dx), 0, 1)
    maps = valais_model_maps(R; dx, dz_fine = 0.05)
    ia = argmin(abs.(maps.x .- SEC_A.E)); jb = argmin(abs.(maps.y .- SEC_B.N))
    A = valais_model_block(maps, ia:ia, eachindex(maps.y)); B = valais_model_block(maps, eachindex(maps.x), jb:jb)
    phcols = ["#ffffff", "#e9c46a", "#8fb8de", "#c0504d", "#7f6a93", "#3d5a40", "#2a9d8f"]
    phrgba(p) = p == 0 ? RGBAf(1, 1, 1, 0) : (c = Makie.to_color(phcols[p + 1]); RGBAf(c.r, c.g, c.b, 1))
    zsurf = replace(topo, NaN => 0.0); zsurf[(x .> SEC_A.E) .& (y' .< SEC_B.N)] .= NaN   # cut the corner facing the viewer
    srgb = [t.level2[I] > 0 ? shadecol(tcol[t.level2[I]], hs2[I]) : grey(hs2[I]) for I in CartesianIndices(topo)]
    BW = ice_water_thickness(R, topo; dx)
    srgb = [BW.water[I] > 0.002 ? shadecol(lakecol, hs2[I]) : BW.ice[I] > 0.01 ? iceshade(hs2[I]) : srgb[I] for I in CartesianIndices(topo)]
    # hypocentres since 1975 (M ≥ 1) within 5 km of the cut faces, projected onto them; historical earthquakes of Mw ≥ 5.5 on the surface
    SQ = fetch_seismicity()
    q = (SQ.time .>= "1975") .& .!isnan.(SQ.z) .& (SQ.mag .>= 1.0)
    ka = findall(q .& (abs.(SQ.E .- SEC_A.E) .< 5) .& (SQ.N .<= SEC_B.N)); kb = findall(q .& (abs.(SQ.N .- SEC_B.N) .< 5) .& (SQ.E .>= SEC_A.E))
    qpts = [Point3f.(SEC_A.E + 0.4, SQ.N[ka], SQ.z[ka]); Point3f.(SQ.E[kb], SEC_B.N - 0.4, SQ.z[kb])]; qmag = [SQ.mag[ka]; SQ.mag[kb]]
    kh = findall((SQ.time .< "1975") .& (SQ.mag .>= 5.5))
    topoat(e, n) = replace(topo, NaN => 0.0)[argmin(abs.(x .- e)), argmin(abs.(y .- n))]
    za = maps.z; jA = findall(<(SEC_B.N), maps.y); iB = findall(>(SEC_A.E), maps.x)
    fence = render3d("fence", (1000, 900), ax -> begin
        ax.aspect = (130, 95, 50); ax.azimuth = -0.33π; ax.elevation = 0.13π
        surface!(ax, x, y, zsurf, color = srgb, shading = NoShading)
        scatter!(ax, qpts, color = quakecol, markersize = 2.5 .+ 2.2 .* qmag,
                 strokecolor = (:black, 0.6), strokewidth = 0.5)
        scatter!(ax, [Point3f(SQ.E[k], SQ.N[k], topoat(SQ.E[k], SQ.N[k]) + 1.0) for k in kh], marker = :star5,
                 markersize = 7 .* SQ.mag[kh] .- 18, color = :red, strokecolor = :black, strokewidth = 0.8)
        surface!(ax, fill(maps.x[ia], length(jA), length(za)), repeat(maps.y[jA], 1, length(za)), repeat(za', length(jA), 1),
                 color = phrgba.(A.Phases[1, jA, :]), shading = NoShading)
        surface!(ax, repeat(maps.x[iB], 1, length(za)), fill(maps.y[jb], length(iB), length(za)), repeat(za', length(iB), 1),
                 color = phrgba.(B.Phases[iB, 1, :]), shading = NoShading)
        zlims!(ax, -60, 5)
    end)

    # ---- composition (CairoMakie, LaTeX font)
    fig = Figure(size = (W, 3000), fonts = tex, backgroundcolor = :white, figure_padding = (36, 36, 24, 24))
    Label(fig[1, 1:2], "GeoValaisWallis3D.jl", fontsize = 64, font = :bold, color = ink, halign = :left, tellwidth = false)
    Label(fig[2, 1:2], "An open 3D geological and geophysical model of Valais/Wallis, built in Julia with GeophysicalModelGenerator.jl",
          fontsize = 27, color = muted, halign = :left, tellwidth = false)
    cap(pos, s) = Label(pos, s, fontsize = 27, color = ink, halign = :left, justification = :left, tellwidth = false, word_wrap = true)

    cap(fig[3, 1:2], "Seismic subsoil classes of the Rhône valley (CREALP / BAFU) on the swissALTIRegio relief, Martigny – Sion – Sierre")
    _, hsz = imgaxis!(fig[4, 1:2], hero)
    Legend(fig[5, 1:2], [PolyElement(color = c) for c in clcols], ["A rock", "B", "C", "D", "E thin soft layer", "F special"],
           "SIA 261 class:", orientation = :horizontal, titleposition = :left, framevisible = false, labelsize = 24, titlesize = 24,
           titlefont = :regular, colgap = 14, tellwidth = false)

    cap(fig[6, 1], "The 3D model: tectonic map, glaciers and lakes on the relief, crust and Moho to 60 km, earthquakes of the SED catalogue")
    _, fsz = imgaxis!(fig[7, 1], fence)
    Legend(fig[8, 1], [[PolyElement(color = phcols[k]) for k in 3:6]; PolyElement(color = icecol, strokecolor = :grey60, strokewidth = 0.5);
                       PolyElement(color = lakecol); MarkerElement(marker = :circle, color = quakecol, strokecolor = :black, strokewidth = 0.5, markersize = 16);
                       MarkerElement(marker = :star5, color = :red, strokecolor = :black, strokewidth = 0.5, markersize = 20)],
           ["sedimentary nappes", "upper crust", "lower crust", "mantle", "glaciers", "lakes", "earthquakes since 1975 (within 5 km)", "historical Mw ≥ 5.5"],
           orientation = :horizontal, framevisible = false, labelsize = 24, colgap = 14, nbanks = 2, tellwidth = false)

    # gravity inversion
    k = argmin(abs.(dres.zg .+ 16.0))
    inner = (R.E[1] .<= dres.E .<= R.E[2]) .& (R.N[1] .<= dres.N .<= R.N[2])
    rs = dres.dobs .- dres.dg0 .- dres.offset; rs .-= mean(rs[inner]); rinv = dres.dobs .- dres.dg
    rmsv(a) = round(sqrt(mean(abs2, a)), digits = 1)
    cap(fig[6, 2], "Gravity inversion (GMG adjoint + LSQR)")
    g = fig[7:8, 2] = GridLayout()
    ax1 = Axis(g[1, 1], aspect = DataAspect(), titlesize = 24, titlefont = :regular, titlecolor = ink,
               title = "density change at $(round(Int, -dres.zg[k])) km depth\nBouguer misfit $(rmsv(rs[inner])) → $(rmsv(rinv[inner])) mGal")
    hidedecorations!(ax1)
    hm = heatmap!(ax1, dres.xg, dres.yg, ifelse.(dres.active[:, :, k], dres.drho[:, :, k], NaN), colormap = Reverse(:RdBu), colorrange = (-200, 200))
    lines!(ax1, [R.E[1], R.E[2], R.E[2], R.E[1], R.E[1]], [R.N[1], R.N[1], R.N[2], R.N[2], R.N[1]], color = ink, linestyle = :dash)
    for (n, e, nn) in PLACES[[3, 6]]
        scatter!(ax1, [e], [nn], color = ink, markersize = 9)
        text!(ax1, e, nn, text = n, fontsize = 22, offset = (6, 4), color = ink, glowcolor = (:white, 0.8), glowwidth = 3)
    end
    text!(ax1, 2648.0, 1068.0, text = "Ivrea body", fontsize = 22, align = (:center, :top), color = ink, glowcolor = (:white, 0.8), glowwidth = 3)
    Colorbar(g[2, 1], hm, label = "Δρ (kg/m³)", vertical = false, flipaxis = false, width = Relative(0.8), labelsize = 22, ticklabelsize = 20)

    # Morcles section
    m = read_morcles()
    jc = argmin(abs.(m.y .- (m.y[1] + m.y[end]) / 2))
    nU = length(MORCLES_UNITS)
    mcm = cgrad([Makie.to_color("#9e6b5a"); [get(cgrad(:roma), (k - 1) / (nU - 2)) for k in 1:nU-1]], nU, categorical = true)
    cap(fig[9, 1:2], "Nappe de Morcles 3D geological model (Thornton et al. 2018), 18 formations: W–E section at N = $(round(m.y[jc], digits = 1)) km")
    axm = Axis(fig[10, 1:2], aspect = DataAspect(), xlabel = "E (km)", ylabel = "z (km)", xlabelsize = 22, ylabelsize = 22,
               xticklabelsize = 20, yticklabelsize = 20, xgridvisible = false, ygridvisible = false)
    heatmap!(axm, m.x, m.z, nanmask(m.unit[:, jc, :], m.unit[:, jc, :] .== 0), colormap = mcm, colorrange = (0.5, nU + 0.5))
    xlims!(axm, extrema(m.x)...); ylims!(axm, 0, 2.8)
    grp = [("crystalline basement", 1), ("Triassic–Jurassic", 6), ("Cretaceous", 12), ("Eocene–Oligocene", 17)]
    axislegend(axm, [PolyElement(color = mcm[k]) for (_, k) in grp], first.(grp), position = :lb, framevisible = true,
               framecolor = (:white, 0), backgroundcolor = (:white, 0.85), labelsize = 22, rowgap = 2, margin = (8, 8, 8, 8))

    Label(fig[11, 1:2], "Data: swisstopo (swissALTIRegio, swissBEDROCK, swissBATHY3D, Tectonic Map 1:500 000, Gravimetric Atlas), BAFU / CREALP, 3D-ALPS " *
                        "(Spooner et al. 2019), Spada et al. (2013), Thornton et al. (2018), Grab et al. (2021), Swiss Seismological Service.    github.com/tduretz/GeoValaisWallis3D.jl",
          fontsize = 18, color = muted, halign = :left, justification = :left, tellwidth = false, word_wrap = true)
    colsize!(fig.layout, 1, Relative(0.6))
    cw = W - 72                                   # content width (figure padding)
    rowsize!(fig.layout, 4, Fixed(cw * hsz[2] / hsz[1]))
    rowsize!(fig.layout, 7, Fixed(0.6cw * fsz[2] / fsz[1]))
    rowsize!(fig.layout, 10, Fixed((cw - 70) * 2.8 / (m.x[end] - m.x[1]) + 60))
    rowgap!(fig.layout, 6)
    for r in (2, 4, 7, 9); rowgap!(fig.layout, r, 18); end
    resize_to_layout!(fig)
    savefig("showcase", fig)
end
