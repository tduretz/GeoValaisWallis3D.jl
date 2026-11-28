# Figures of the model on the curvilinear mesh of GeophysicalModelGenerator (curvilinear_grid): the grid levels follow the
# topography, the base of the consolidated sediments, the Conrad and the European Moho.
#     julia --project=figures -t auto figures/curvilinear/make_figures.jl [figure ...]
# Figures: overview mesh layers sections subsoil gravity density glaciers morcles showcase (default: all).
# Output: figures/curvilinear/*.jpg (needs ImageMagick). overview and showcase render their 3D panels with GLMakie (needs a
# display). gravity and density use GMG voxel_grav on curvilinear grids; GVW_BACKEND=metal runs them on the Apple GPU
# (Float32), otherwise on the CPU (Float64, slow). Their results are cached in data/output/curvilinear_*.jls (delete to rerun).
# The figures of the rectilinear grid are made by figures/rectilinear/make_figures.jl and are not changed by this script.
using GeoValaisWallis3D, CairoMakie, Statistics, Serialization

const FIGDIR = @__DIR__
mkpath(FIGDIR)
const ALL = ["overview", "mesh", "layers", "sections", "subsoil", "gravity", "density", "glaciers", "morcles", "showcase"]
const WHICH = isempty(ARGS) ? ALL : ARGS
const R = valais_region()
CairoMakie.activate!(type = "png", px_per_unit = 1.25)
set_theme!(fontsize = 13, Axis = (titlesize = 14, titlefont = :regular, xgridvisible = false, ygridvisible = false),
           Axis3 = (titlefont = :regular,))
const METAL = get(ENV, "GVW_BACKEND", "cpu") == "metal"
METAL && @eval using Metal
backend_ft() = METAL ? (Base.invokelatest(() -> Metal.MetalBackend()), Float32) : (nothing, Float64)

const PLACES = [("Monthey", 2562.0, 1122.5), ("Martigny", 2571.5, 1105.5), ("Sion", 2594.0, 1120.0),
                ("Sierre", 2607.0, 1127.0), ("Visp", 2634.5, 1127.0), ("Brig", 2642.0, 1129.5),
                ("Zermatt", 2624.5, 1097.0), ("Leukerbad", 2614.5, 1136.5)]
const SEC_A = (name = "A", E = 2594.0)          # N–S through Sion
const SEC_B = (name = "B", N = 1120.0)          # W–E along the Rhône valley at Sion
const PHCOLS = ["#ffffff", "#e9c46a", "#8fb8de", "#c0504d", "#7f6a93", "#3d5a40", "#2a9d8f"]
const PHNAMES = ["air", "Quaternary", "consolidated sediments", "upper crust", "lower crust", "mantle", "Adriatic mantle"]
const PHCM = cgrad(PHCOLS, 7, categorical = true)
# the meshed surfaces: label and line colour
const SURF = (topography = ("topography", :black), fine_base = ("base of the 25 m layer", "#6b6b6b"),
              base_cs = ("base of consolidated sediments", "#1f4e79"), conrad = ("Conrad", "#7a1f1f"),
              moho = ("European Moho", "#14361d"), bottom = ("bottom", :black))

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
function savefig(name, fig)
    png = joinpath(tempdir(), "curvilinear_" * name * ".png"); save(png, fig)
    jpg = joinpath(FIGDIR, name * ".jpg")
    run(`magick $png -quality 85 $jpg`)
    println("written figures/curvilinear/$name.jpg ($(round(filesize(jpg) / 1e6, digits = 2)) MB)")
end
nanmask(a, m) = (b = Float64.(a); b[m] .= NaN; b)
sq(a) = dropdims(a, dims = findfirst(==(1), size(a)))

# A field V at the nodes of a curvilinear section (horizontal positions h, heights Z[i, k]) drawn as the dual cells of the
# nodes, each in one colour: corners midway between the nodes, so that cells follow the mesh exactly.
function curvheatmap!(ax, h, Z, V; kw...)
    n, m = size(Z)
    hc = vcat(h[1], (h[1:end-1] .+ h[2:end]) ./ 2, h[end])
    zc = [(Z[clamp(i - 1, 1, n), clamp(k - 1, 1, m)] + Z[clamp(i, 1, n), clamp(k - 1, 1, m)] +
           Z[clamp(i - 1, 1, n), clamp(k, 1, m)] + Z[clamp(i, 1, n), clamp(k, 1, m)]) / 4 for i in 1:n+1, k in 1:m+1]
    pts = Point2f[]; faces = Makie.GLTriangleFace[]; col = Float32[]
    for k in 1:m, i in 1:n
        v = V[i, k]; isnan(v) && continue
        b = length(pts)
        push!(pts, Point2f(hc[i], zc[i, k]), Point2f(hc[i + 1], zc[i + 1, k]), Point2f(hc[i + 1], zc[i + 1, k + 1]), Point2f(hc[i], zc[i, k + 1]))
        push!(faces, Makie.GLTriangleFace(b + 1, b + 2, b + 3), Makie.GLTriangleFace(b + 1, b + 3, b + 4))
        append!(col, (v, v, v, v))
    end
    mesh!(ax, pts, faces; color = col, shading = NoShading, kw...)
end
# grid lines of a curvilinear section: every `every`-th level and column (thin), the meshed surfaces (bold)
function meshlines!(ax, h, Z, levels; every = 1, cols = every, color = (:black, 0.35), linewidth = 0.4, surfaces = true)
    for k in 1:every:size(Z, 2)
        lines!(ax, h, Z[:, k]; color, linewidth)
    end
    for i in 1:cols:length(h)
        lines!(ax, [h[i], h[i]], [Z[i, 1], Z[i, end]]; color, linewidth)
    end
    surfaces || return
    for (s, l) in pairs(levels)
        s in (:bottom,) && continue
        lines!(ax, h, Z[:, l], color = SURF[s][2], linewidth = s == :fine_base ? 1.0 : 1.8, linestyle = s == :fine_base ? :dash : :solid)
    end
end
# horizontal slice at height zq of a field F at the nodes of a curvilinear grid (heights Z), linear along the columns
function zslice(Z, F, zq)
    out = fill(NaN, size(Z, 1), size(Z, 2))
    for j in axes(Z, 2), i in axes(Z, 1)
        c = @view Z[i, j, :]
        (zq < c[1] || zq > c[end]) && continue
        k = clamp(searchsortedlast(c, zq), 1, length(c) - 1)
        w = (zq - c[k]) / (c[k + 1] - c[k])
        out[i, j] = (1 - w) * F[i, j, k] + w * F[i, j, k + 1]
    end
    out
end
surflegend(levels; kw...) = ([LineElement(color = SURF[s][2], linewidth = 2, linestyle = s == :fine_base ? :dash : :solid) for s in keys(levels) if s != :bottom],
                             [SURF[s][1] for s in keys(levels) if s != :bottom])

# the model on the curvilinear mesh (250 m columns); built once, used by most figures
const NEEDMAPS = any(in(WHICH), ["overview", "mesh", "layers", "sections", "subsoil", "showcase"])
const MAPS = NEEDMAPS ? valais_curvilinear_maps(R; dx = 0.25) : nothing
col_index(v, a) = argmin(abs.(v .- a))
section_A(maps; kw...) = (i = col_index(maps.x, SEC_A.E); (; i, h = maps.y, Z = maps.z[i, :, :], F = valais_model_block(maps, i:i, eachindex(maps.y); kw...)))
section_B(maps; kw...) = (j = col_index(maps.y, SEC_B.N); (; j, h = maps.x, Z = maps.z[:, j, :], F = valais_model_block(maps, eachindex(maps.x), j:j; kw...)))

# --------------------------------------------------------------------------------------------- 3D rendering (GLMakie)
const TEX = (; regular = Makie.texfont(), bold = Makie.texfont(:bold), italic = Makie.texfont(:italic))
shadecol(c, g) = (k = Makie.to_color(c); RGBf(k.r * (0.45 + 0.6g), k.g * (0.45 + 0.6g), k.b * (0.45 + 0.6g)))
grey(g) = RGBf(0.3 + 0.65g, 0.3 + 0.65g, 0.3 + 0.65g)
const TMP = mktempdir()
function render3d(name, sz, draw!)
    @eval using GLMakie
    Base.invokelatest() do
        GLMakie.activate!(px_per_unit = 2)
        f = Figure(size = sz, fonts = TEX, figure_padding = 0)
        ax = Axis3(f[1, 1]; protrusions = 0, viewmode = :fit)
        hidedecorations!(ax); hidespines!(ax)
        draw!(ax)
        png = joinpath(TMP, name * ".png"); save(png, f)
        run(`magick $png -trim +repage $png`)
        CairoMakie.activate!(type = "png", px_per_unit = 1.25)
        png
    end
end
function imgaxis!(pos, f)
    img = rotr90(Makie.FileIO.load(f))
    ax = Axis(pos, aspect = DataAspect(), backgroundcolor = :transparent)
    hidedecorations!(ax); hidespines!(ax)
    image!(ax, img, interpolate = true)
    ax, size(img)
end
# fence diagram of the curvilinear model: tectonic map on the relief (SE quadrant cut away), the phases on the two cut faces
# with the mesh lines (every `every`-th level) and the meshed surfaces
function render_fence(name, sz; every = 6, linesurf = true)
    maps = MAPS
    x, y = maps.x, maps.y
    topo = maps.z[:, :, end]
    t = tectonic_map_on_grid(x, y)
    used = sort(unique(t.level2[t.level2 .> 0])); tcols = Makie.to_colormap(:tab20)
    tcol = Dict(c => tcols[k] for (k, c) in enumerate(used))
    hs = clamp.(hillshade(topo, 0.25), 0, 1)
    A = section_A(maps); B = section_B(maps)
    ia, jb = A.i, B.j
    jA = findall(<=(y[jb]), y); iB = findall(>=(x[ia]), x)
    phrgba(p) = (c = Makie.to_color(PHCOLS[p + 1]); RGBAf(c.r, c.g, c.b, 1))
    zsurf = copy(topo); zsurf[(x .> x[ia]) .& (y' .< y[jb])] .= NaN
    srgb = [t.level2[I] > 0 ? shadecol(tcol[t.level2[I]], hs[I]) : grey(hs[I]) for I in CartesianIndices(topo)]
    render3d(name, sz, ax -> begin
        ax.aspect = (130, 95, 50); ax.azimuth = -0.33π; ax.elevation = 0.13π
        surface!(ax, x, y, zsurf, color = srgb, shading = NoShading)
        nz = size(maps.z, 3)
        # face A (x = x[ia], south part), face B (y = y[jb], east part): the mesh nodes are the surface vertices
        surface!(ax, fill(x[ia], length(jA), nz), repeat(y[jA], 1, nz), A.Z[jA, :], color = phrgba.(A.F.Phases[1, jA, :]), shading = NoShading)
        surface!(ax, repeat(x[iB], 1, nz), fill(y[jb], length(iB), nz), B.Z[iB, :], color = phrgba.(B.F.Phases[iB, 1, :]), shading = NoShading)
        ox, oy = 0.15, -0.15                                     # lines slightly in front of the faces
        for k in 1:every:nz
            lines!(ax, fill(x[ia] + ox, length(jA)), y[jA], A.Z[jA, k], color = (:black, 0.28), linewidth = 0.6)
            lines!(ax, x[iB], fill(y[jb] + oy, length(iB)), B.Z[iB, k], color = (:black, 0.28), linewidth = 0.6)
        end
        if linesurf
            for (s, l) in pairs(maps.levels)
                s in (:bottom, :topography, :fine_base) && continue
                lines!(ax, fill(x[ia] + ox, length(jA)), y[jA], A.Z[jA, l], color = SURF[s][2], linewidth = 2.5)
                lines!(ax, x[iB], fill(y[jb] + oy, length(iB)), B.Z[iB, l], color = SURF[s][2], linewidth = 2.5)
            end
        end
        zlims!(ax, -70, 5)
    end)
end

# --------------------------------------------------------------------------------------------- overview (3D)
if "overview" in WHICH
    fence = render_fence("overview_fence", (1400, 1150); every = 4)
    fig = Figure(size = (1500, 1100), fonts = TEX, figure_padding = (24, 24, 16, 16))
    Label(fig[1, 1], "The Valais model on a curvilinear mesh (GeophysicalModelGenerator curvilinear_grid): " *
                     "$(length(MAPS.x)) × $(length(MAPS.y)) × $(size(MAPS.z, 3)) nodes, every 4th level drawn",
          fontsize = 20, halign = :left, tellwidth = false)
    _, fsz = imgaxis!(fig[2, 1], fence)
    els, labs = surflegend((; base_cs = 0, conrad = 0, moho = 0))
    Legend(fig[3, 1], [[PolyElement(color = PHCOLS[k]) for k in 2:6]; els], [PHNAMES[2:6]; labs], orientation = :horizontal,
           framevisible = false, labelsize = 17, colgap = 12, nbanks = 2, tellwidth = false)
    Label(fig[4, 1], "Grid levels follow the topography, the base of the consolidated sediments, the Conrad and the European Moho; " *
                     "the Quaternary fill lies in the top 1 km with 25 m spacing. Surface: tectonic map (swisstopo 1:500 000) on the relief.",
          fontsize = 15, color = :grey30, halign = :left, justification = :left, word_wrap = true, tellwidth = false)
    rowsize!(fig.layout, 2, Fixed((1500 - 48) * fsz[2] / fsz[1]))
    resize_to_layout!(fig)
    savefig("overview", fig)
end

# --------------------------------------------------------------------------------------------- the mesh (sections)
if "mesh" in WHICH
    maps = MAPS
    B = section_B(maps); A = section_A(maps)
    fig = Figure(size = (1500, 1300))
    ax1 = Axis(fig[1, 1:2], title = "B–B' (W–E, N = $(SEC_B.N) km): the curvilinear mesh (every 4th level and column) on the phases",
               ylabel = "elevation (km)", xlabel = "E (km)")
    curvheatmap!(ax1, B.h, B.Z, Float64.(sq(B.F.Phases)), colormap = PHCM, colorrange = (-0.5, 6.5))
    meshlines!(ax1, B.h, B.Z, maps.levels; every = 4)
    xlims!(ax1, extrema(B.h)...); ylims!(ax1, -70, 5)
    # zoom: Rhône valley at Sion, every level and column
    zx = (2586.0, 2602.0); ix = findall(e -> zx[1] <= e <= zx[2], B.h)
    ax2 = Axis(fig[2, 1], title = "zoom: Rhône valley at Sion, every level and column (vertical exaggeration 2)",
               ylabel = "elevation (km)", xlabel = "E (km)")
    curvheatmap!(ax2, B.h[ix], B.Z[ix, :], Float64.(sq(B.F.Phases)[ix, :]), colormap = PHCM, colorrange = (-0.5, 6.5))
    meshlines!(ax2, B.h[ix], B.Z[ix, :], maps.levels; every = 1, linewidth = 0.3)
    xlims!(ax2, zx...); ylims!(ax2, -1.2, 1.6)
    ax2.aspect = AxisAspect((zx[2] - zx[1]) / (2 * 2.8))
    # vertical spacing in a column (Sion)
    i = col_index(B.h, 2594.0); z = B.Z[i, :]; dz = diff(z); zm = (z[1:end-1] .+ z[2:end]) ./ 2
    ax3 = Axis(fig[2, 2], title = "vertical spacing in the column at Sion", xlabel = "cell height (m)", ylabel = "elevation (km)",
               xscale = log10, xticks = [10, 25, 100, 250, 500, 1000])
    lines!(ax3, 1e3 .* dz, zm, color = :black, linewidth = 1.5)
    for (s, l) in pairs(maps.levels)
        s == :bottom && continue
        hlines!(ax3, [z[l]], color = SURF[s][2], linestyle = s == :fine_base ? :dash : :solid, linewidth = 1.5)
    end
    ylims!(ax3, -70, 5)
    els, labs = surflegend(maps.levels)
    Legend(fig[1, 3], [[PolyElement(color = PHCOLS[k]) for k in 2:6]; els], [PHNAMES[2:6]; labs], framevisible = false)
    Label(fig[0, 1:3], "Curvilinear mesh of GeophysicalModelGenerator: $(length(maps.x)) × $(length(maps.y)) × $(size(maps.z, 3)) nodes, " *
                       "cells per layer $(join(maps.nz, ", ")) (top to bottom: 25 m, 100 m, 500 m, 500 m, 1 km where the layer is thickest)", fontsize = 14)
    colsize!(fig.layout, 2, Relative(0.25))
    savefig("mesh", fig)
end

# ------------------------------------------------------------------------------------- the meshed surfaces (maps)
if "layers" in WHICH
    maps = MAPS; x, y = maps.x, maps.y; I = maps.interfaces; L = maps.levels; Z = maps.z
    topo = Z[:, :, L.topography]
    moved(s) = abs.(Z[:, :, L[s]] .- maps.surfaces[s]) .> 1e-3
    wedge = nanmask(I.moho_ad, .!(I.moho_ad .> I.moho_eu))
    fig = Figure(size = (1500, 760))
    panels = [
        ("Mesh level $(L.topography): topography (swissALTIRegio)", topo, :batlowW, (0.3, 4.2), "elevation (km)", nothing),
        ("Quaternary thickness (swissBEDROCK), in the 25 m layer", nanmask(1e3 .* I.quat_thickness, I.quat_thickness .<= 0.002), Reverse(:lapaz), (0, 600), "thickness (m)", nothing),
        ("Mesh level $(L.base_cs): base of consolidated sediments", Z[:, :, L.base_cs], :batlow, (-12, 3), "elevation (km)", :base_cs),
        ("Mesh level $(L.conrad): Conrad", Z[:, :, L.conrad], :batlow, (-35, -10), "elevation (km)", nothing),
        ("Mesh level $(L.moho): European Moho (Spada et al. 2013)", Z[:, :, L.moho], :batlow, (-60, -28), "elevation (km)", nothing),
        ("Adriatic Moho: not meshed (crosses the Conrad), a material wedge", wedge, :batlow, (-25, -15), "elevation (km)", nothing),
    ]
    for (k, (title, A, cm, cr, lab, mv)) in enumerate(panels)
        r, c = (k - 1) ÷ 3 + 1, (k - 1) % 3 + 1
        g = fig[r, c] = GridLayout()
        ax = mapaxis(g[1, 1]; title)
        shade!(ax, x, y, topo)
        hm = heatmap!(ax, x, y, A, colormap = cm, colorrange = cr, alpha = k == 1 ? 0.55 : 0.8, nan_color = :transparent)
        k >= 3 && contour!(ax, x, y, A, levels = range(cr..., step = k >= 5 ? 5 : 2.5), color = (:black, 0.5), linewidth = 0.6)
        isnothing(mv) || contour!(ax, x, y, Float64.(moved(mv)), levels = [0.5], color = :red, linewidth = 1.0)
        places!(ax; names = false)
        limits_region!(ax)
        Colorbar(g[1, 2], hm, label = lab, height = Relative(0.7))
        r == 1 || hidexdecorations!(ax, label = false, ticklabels = false, ticks = false)
    end
    nmv = count(moved(:base_cs))
    Label(fig[0, 1:3], "Surfaces of the curvilinear mesh: topography, Conrad and Moho are grid levels exactly; the base of the consolidated sediments too, " *
                       "except where they pinch out (red outline, $(round(100nmv / length(topo), digits = 1)) % of the columns: level ≤ $(round(Int, 1e3maximum(abs.(Z[:, :, L.base_cs] .- maps.surfaces.base_cs)))) m below)", fontsize = 14)
    savefig("layers", fig)
end

# ------------------------------------------------------------------------------------------ cross-sections
if "sections" in WHICH
    maps = MAPS
    A = section_A(maps; temperature = true); B = section_B(maps; temperature = true)
    fig = Figure(size = (1500, 1250))
    secs = [(A, "A–A' (N–S, E = $(SEC_A.E) km)", "N (km)"), (B, "B–B' (W–E, N = $(SEC_B.N) km)", "E (km)")]
    local h2, h3
    for (c, (S, name, xl)) in enumerate(secs)
        ph = Float64.(sq(S.F.Phases)); h = S.h; Zs = S.Z
        ax1 = Axis(fig[1, c], title = "$name: phases and meshed surfaces", ylabel = "elevation (km)")
        curvheatmap!(ax1, h, Zs, ph, colormap = PHCM, colorrange = (-0.5, 6.5))
        meshlines!(ax1, h, Zs, maps.levels; every = 1000, cols = 10^6)
        ax2 = Axis(fig[2, c], title = "$name: Vs (first-order velocity laws)", ylabel = "elevation (km)")
        h2 = curvheatmap!(ax2, h, Zs, Float64.(sq(S.F.Vs)), colormap = :lipari, colorrange = (0.3, 4.6))
        ax3 = Axis(fig[3, c], title = "$name: temperature (3D-ALPS-TR)", ylabel = "elevation (km)")
        T = Float64.(sq(S.F.Temperature))
        h3 = curvheatmap!(ax3, h, Zs, T, colormap = :lajolla, colorrange = (0, 1200))
        ax4 = Axis(fig[4, c], title = "$name: near the surface, every 4th level (vertically stretched)", ylabel = "elevation (km)", xlabel = xl)
        curvheatmap!(ax4, h, Zs, ph, colormap = PHCM, colorrange = (-0.5, 6.5))
        meshlines!(ax4, h, Zs, maps.levels; every = 4, cols = 10^6, color = (:black, 0.12), linewidth = 0.25)
        for ax in (ax1, ax2, ax3)
            ylims!(ax, -70, 5); xlims!(ax, extrema(h)...)
        end
        xlims!(ax4, extrema(h)...); ylims!(ax4, -1.0, 4.6)
    end
    Colorbar(fig[2, 3], h2, label = "Vs (km/s)"); Colorbar(fig[3, 3], h3, label = "T (°C)")
    els, labs = surflegend(maps.levels)
    Legend(fig[1, 3], [[PolyElement(color = PHCOLS[k], strokecolor = :grey60, strokewidth = 0.5) for k in 2:7]; els], [PHNAMES[2:7]; labs], framevisible = false, labelsize = 11, rowgap = 0)
    Legend(fig[4, 3], [PolyElement(color = PHCOLS[k], strokecolor = :grey60, strokewidth = 0.5) for k in 2:4], PHNAMES[2:4], framevisible = false)
    savefig("sections", fig)
end

# ------------------------------------------------------------------------------ seismic subsoil classes (CREALP)
if "subsoil" in WHICH
    maps = MAPS
    sm = subsoil_map(R; dx = 0.05)
    clcols = ["#3b5b92", "#5aa0c8", "#9ccf7a", "#f2c14e", "#e8743b", "#b5179e"]
    clcm = cgrad(clcols, 6, categorical = true)
    cls = nanmask(sm.class, sm.class .== 0)
    fig = Figure(size = (1500, 1500))
    ax = mapaxis(fig[1, 1:3], title = "Seismic subsoil classes SIA 261 (BAFU; Rhône valley 1:25 000 maps of CREALP / canton of Valais), release $(sm.release); " *
                                       "N–S sections of the curvilinear model below (dashed)")
    shade!(ax, sm.x, sm.y, sm.topo)
    heatmap!(ax, sm.x, sm.y, cls, colormap = clcm, colorrange = (0.5, 6.5), nan_color = :transparent)
    contour!(ax, sm.x, sm.y, Float64.(sm.ssm .> 0), levels = [0.5], color = :black, linewidth = 1.2)
    secs = [("Martigny", 2571.5, (1100.0, 1112.0)), ("Sion", 2594.0, (1113.0, 1127.0)), ("Visp", 2634.5, (1121.0, 1133.0))]
    for (n, e, yr) in secs
        lines!(ax, [e, e], collect(yr), color = :black, linestyle = :dash, linewidth = 2)
    end
    places!(ax); limits_region!(ax)
    Legend(fig[1, 4], [[PolyElement(color = c) for c in clcols]; LineElement(color = :black)],
           [["$(SUBSOIL_CLASSES[k])" * (k == 1 ? " (rock)" : k == 5 ? " (thin soft layer on rock)" : k == 6 ? " (special)" : "") for k in 1:6];
            "spectral microzonation (classes replaced)"], "class", framevisible = false)
    local hv
    for (c, (n, e, yr)) in enumerate(secs)
        i = col_index(maps.x, e); j = findall(v -> yr[1] <= v <= yr[2], maps.y)
        F = valais_model_block(maps, i:i, j)
        Zs = maps.z[i, j, :]; h = maps.y[j]
        axs = Axis(fig[2, c], title = "$n (E = $(maps.x[i]) km): Vs on the mesh, every level", xlabel = "N (km)", ylabel = c == 1 ? "elevation (km)" : "")
        hv = curvheatmap!(axs, h, Zs, Float64.(sq(F.Vs)), colormap = :lipari, colorrange = (0.3, 3.5))
        meshlines!(axs, h, Zs, maps.levels; every = 1, cols = 10^6, color = (:white, 0.25), linewidth = 0.25)
        # Quaternary base (not meshed: a material boundary inside the 25 m layer)
        bed = maps.interfaces.bed[i, j]; q = maps.interfaces.quat_thickness[i, j] .> 0
        lines!(axs, h, ifelse.(q, bed, NaN), color = :white, linewidth = 1.5)
        # surface subsoil class along the section
        ks = [col_index(sm.x, maps.x[i])]; js = [col_index(sm.y, v) for v in h]
        sc = [sm.class[ks[1], jj] for jj in js]
        scatter!(axs, h, Zs[:, end] .+ 0.12, color = [s == 0 ? RGBAf(0, 0, 0, 0) : Makie.to_color(clcols[s]) for s in sc], marker = :rect, markersize = 9)
        xlims!(axs, yr...); ylims!(axs, -0.8, 3.2)
    end
    Colorbar(fig[2, 4], hv, label = "Vs (km/s)", height = Relative(0.8))
    Label(fig[3, 1:4], "The Quaternary fill (white line: its base, a material boundary) lies in the 25 m layer of the mesh that follows the topography; " *
                       "squares above the surface: subsoil class of the map at the section.", fontsize = 13, halign = :left, tellwidth = false)
    rowsize!(fig.layout, 1, Relative(0.6))
    savefig("subsoil", fig)
end

# --------------------------------------------------------------------------------------------- gravity
# attraction of the topographic masses (sea level to topography, 2670 kg/m³) with GMG voxel_grav on curvilinear grids whose
# top level is the topography: 250 m columns out to 10 km around the region, 1 km columns (hollow inside) out to 100 km
if "gravity" in WHICH
    cache = datadir("output", "curvilinear_gravity_topo.jls")
    res = if isfile(cache)
        deserialize(cache)
    else
        backend, ft = backend_ft()
        grow(m) = ValaisRegion((R.E[1] - m, R.E[2] + m), (R.N[1] - m, R.N[2] + m), R.z)
        near = grow(10.0)
        # stations every 1 km, on the 250 m topography of the mesh (bilinear between its columns)
        xs, ys, _ = fetch_topography(R; dx = 1.0)
        xn, yn, tn = fetch_topography(near; dx = 0.25); tn = max.(replace(tn, NaN => 0.0), 0.005)
        function bil(e, n)
            i = clamp(searchsortedlast(xn, e), 1, length(xn) - 1); j = clamp(searchsortedlast(yn, n), 1, length(yn) - 1)
            a = (e - xn[i]) / (xn[i + 1] - xn[i]); b = (n - yn[j]) / (yn[j + 1] - yn[j])
            (1 - a) * (1 - b) * tn[i, j] + a * (1 - b) * tn[i + 1, j] + (1 - a) * b * tn[i, j + 1] + a * b * tn[i + 1, j + 1]
        end
        Xs, Ys = [a for a in xs, _ in ys], [b for _ in xs, b in ys]
        ts = bil.(Xs, Ys)
        St = CartData(reshape(Xs, :, 1, 1), reshape(Ys, :, 1, 1), reshape(ts, :, 1, 1), (h = reshape(ts, :, 1, 1),))
        function layer_gravity(r, dx, nz, hole)
            x, y, t = fetch_topography(r; dx); t = max.(replace(t, NaN => 0.0), 0.005)
            G = curvilinear_grid(x, y, (GeoValaisWallis3D._surface(x, y, t), 0.0), (nz,))
            X, Y, _ = xyz_grid(x, y, 1:nz + 1)
            Grid = CartData(X, Y, G.z.val, (z = G.z.val,))
            rho = fill(2670.0, size(X))
            isnothing(hole) || (rho[(hole.E[1] .< X .< hole.E[2]) .& (hole.N[1] .< Y .< hole.N[2])] .= 0.0)
            @info "topographic layer: $(size(X)) nodes"
            vec(voxel_grav(Grid, rho; stations = St, refMod = zeros(nz + 1), rhoTol = 0, backend, float_type = ft)[1])
        end
        t0 = time()
        g = layer_gravity(near, 0.25, 4, nothing) .+ layer_gravity(grow(100.0), 1.0, 2, (E = near.E, N = near.N))
        @info "GMG topographic attraction in $(round(time() - t0, digits = 1)) s"
        gp = topographic_attraction(vec(Xs), vec(Ys), vec(ts), R)
        r = (; x = xs, y = ys, topo = ts, ggmg = reshape(g, size(Xs)), gprism = reshape(gp, size(Xs)))
        serialize(cache, r); r
    end
    gm = gravity_maps(R; dx = 1.0)
    fig = Figure(size = (1500, 900))
    d = res.ggmg .- res.gprism
    panels = [
        ("Complete Bouguer anomaly, atlas 1:100 000 (points)", gm.bouguer, (-200, 200), Reverse(:RdBu), "mGal"),
        ("Attraction of the topography: GMG voxel_grav on curvilinear layers (2670 kg/m³)", res.ggmg, (0, 300), :lajolla, "mGal"),
        ("GMG (flat Earth) minus prisms and line masses of the package (with Earth curvature)", d, (-1, 1) .* max(5.0, maximum(abs, d)), Reverse(:RdBu), "mGal"),
        ("Free-air anomaly = Bouguer + GMG topography", gm.bouguer .+ res.ggmg, (-200, 200), Reverse(:RdBu), "mGal"),
    ]
    for (k, (title, A, cr, cm, lab)) in enumerate(panels)
        g = fig[(k - 1) ÷ 2 + 1, (k - 1) % 2 + 1] = GridLayout()
        local ax = mapaxis(g[1, 1]; title, titlesize = 13)
        shade!(ax, res.x, res.y, res.topo)
        hm = heatmap!(ax, res.x, res.y, A, colormap = cm, colorrange = cr, alpha = 0.8, nan_color = :transparent)
        contour!(ax, res.x, res.y, A, levels = k == 3 ? range(-30, 30, step = 2) : range(-300, 300, step = 20), color = (:black, 0.35), linewidth = 0.5)
        k == 1 && scatter!(ax, gm.atlas.E, gm.atlas.N, markersize = 1.5, color = (:black, 0.5))
        places!(ax; names = false); limits_region!(ax)
        Colorbar(g[1, 2], hm, label = lab, height = Relative(0.7))
    end
    Label(fig[0, 1:2], "Gravity on 1 km stations at the topography. GMG − package: mean $(round(mean(filter(isfinite, d)), digits = 1)) mGal, " *
                       "std $(round(std(filter(isfinite, d)), digits = 1)) mGal (the GMG grids are flat; the package includes the Earth's curvature)", fontsize = 14)
    savefig("gravity", fig)
end

# ------------------------------------------------------------------------- density inversion of the Bouguer anomaly
if "density" in WHICH || "showcase" in WHICH
    const SIGMA = 2.0
    cache = datadir("output", "curvilinear_density_inversion.jls")
    DRES = if isfile(cache)
        deserialize(cache)
    else
        backend, ft = backend_ft()
        d = curvilinear_density_model(R; pad = 20.0, dx = 1.0)
        c = curvilinear_bouguer_contrast(d.z, d.rho)
        X, Y, _ = xyz_grid(d.x, d.y, 1:size(d.z, 3))
        Grid = CartData(X, Y, d.z, (rho = d.rho,))
        at = fetch_gravity_atlas(R; margin = 20.0)
        sel = (d.x[1] + 2 .<= at.E .<= d.x[end] - 2) .& (d.y[1] + 2 .<= at.N .<= d.y[end] - 2) .& .!isnan.(at.bouguer)
        E, N, H, dobs = at.E[sel], at.N[sel], at.h[sel], at.bouguer[sel]
        active = trues(size(d.z))
        t0 = time()
        inv = invert_density(Grid, c.contrast, active, (E, N, H), dobs, SIGMA; z0 = 1.0, β = 2.0, backend, float_type = ft)
        r = (; x = d.x, y = d.y, z = d.z, levels = d.levels, nz = d.nz, topo = d.topo, rho0 = d.rho, drho = inv.drho, E, N, dobs,
             dg0 = inv.dg0, dg = inv.dg, offset = inv.offset, iterations = inv.iterations, rms_history = inv.rms_history,
             seconds = time() - t0)
        serialize(cache, r); r
    end
end
if "density" in WHICH
    (; x, y, z, levels, topo, rho0, drho, E, N, dobs) = DRES
    inner = (R.E[1] .<= E .<= R.E[2]) .& (R.N[1] .<= N .<= R.N[2])
    rstart = dobs .- DRES.dg0 .- DRES.offset; rstart .-= mean(rstart[inner])
    rinv = dobs .- DRES.dg
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
        contour!(ax, x, y, topo, levels = [1.0, 2.5], color = (:black, 0.3), linewidth = 0.5); box!(ax)
        Colorbar(g1[1, 2n], sc, label = "mGal", height = Relative(0.6))
    end
    local hm
    for (n, zz) in enumerate([0.0, -5.0, -15.0, -35.0])
        local ax = Axis(g2[1, n], aspect = DataAspect(), title = "density change, z = $(zz) km", xlabel = "E (km)", ylabel = n == 1 ? "N (km)" : "")
        hm = heatmap!(ax, x, y, zslice(z, drho, zz), colormap = Reverse(:RdBu), colorrange = (-200, 200), nan_color = :transparent)
        box!(ax); places!(ax; names = false)
    end
    Colorbar(g2[1, 5], hm, label = "Δρ (kg/m³)", height = Relative(0.6))
    i = col_index(x, SEC_A.E); Zs = z[i, :, :]
    ax = Axis(g3[1, 1], xlabel = "N (km)", ylabel = "elevation (km)", title = "A–A' (E = $(x[i]) km): inverted density, meshed surfaces")
    h1 = curvheatmap!(ax, y, Zs, rho0[i, :, :] .+ drho[i, :, :], colormap = :batlow, colorrange = (2300, 3400))
    meshlines!(ax, y, Zs, levels; every = 1000, cols = 10^6)
    Colorbar(g3[1, 2], h1, label = "ρ (kg/m³)")
    ax = Axis(g3[1, 3], xlabel = "N (km)", title = "A–A': density change, every level of the mesh")
    h2 = curvheatmap!(ax, y, Zs, drho[i, :, :], colormap = Reverse(:RdBu), colorrange = (-200, 200))
    meshlines!(ax, y, Zs, levels; every = 1, cols = 10^6, color = (:black, 0.15), linewidth = 0.3)
    Colorbar(g3[1, 4], h2, label = "Δρ (kg/m³)")
    ax = Axis(g3[1, 5], xlabel = "LSQR iteration", ylabel = "rms misfit (mGal)", yscale = log10, yticks = [2, 5, 10, 20, 50, 100], title = "convergence (stops at σ = $(SIGMA) mGal)")
    lines!(ax, 0:length(DRES.rms_history) - 1, DRES.rms_history, color = "#2a78d6", linewidth = 2)
    hlines!(ax, [SIGMA], color = :gray, linestyle = :dash)
    colsize!(g3, 5, Relative(0.22))
    rowsize!(fig.layout, 1, Relative(0.3)); rowsize!(fig.layout, 2, Relative(0.27))
    Label(fig[0, 1], "Density inversion of the Bouguer anomaly on a curvilinear mesh (GMG voxel_grav + voxel_grav_adjoint, LSQR): $(length(E)) atlas points, " *
                     "$(length(z)) unknowns ($(join(size(z), " × ")), no air nodes), $(DRES.iterations) iterations", fontsize = 14)
    savefig("density", fig)
end

# ------------------------------------------------------------------------------------------- glaciers and lakes
# ice (Swiss Glacier Thickness R2020) and lake water (swissBATHY3D) between the topography and the solid surface: maps, their
# gravity effect at the atlas points (GMG voxel_grav, curvilinear density meshes with and without them), sections
if "glaciers" in WHICH
    PH9 = [PHCOLS; "#bfe3f2"; "#2f6db5"]; PH9CM = cgrad(PH9, 9, categorical = true)
    PH9N = [PHNAMES; "glacier ice"; "lake water"]
    cache = datadir("output", "curvilinear_glaciers_gravity.jls")
    GR = if isfile(cache)
        deserialize(cache)
    else
        backend, ft = backend_ft()
        at = fetch_gravity_atlas(R; margin = 20.0)
        function gravity_of(; kw...)
            d = curvilinear_density_model(R; pad = 20.0, dx = 1.0, kw...)
            sel = (d.x[1] + 2 .<= at.E .<= d.x[end] - 2) .& (d.y[1] + 2 .<= at.N .<= d.y[end] - 2) .& .!isnan.(at.bouguer)
            X, Y, _ = xyz_grid(d.x, d.y, 1:size(d.z, 3))
            St = CartData(reshape(at.E[sel], :, 1, 1), reshape(at.N[sel], :, 1, 1), reshape(at.h[sel], :, 1, 1), (h = reshape(at.h[sel], :, 1, 1),))
            # densities relative to 2670 kg/m³: both meshes have the same outer shape, so the reference cancels in the difference
            g = vec(voxel_grav(CartData(X, Y, d.z, (rho = d.rho,)), d.rho .- 2670.0; stations = St, refMod = zeros(size(d.z, 3)),
                               rhoTol = 0, backend, float_type = ft)[1])
            (; g, sel)
        end
        t0 = time()
        g0 = gravity_of(); g1 = gravity_of(glaciers = true, lakes = true)
        @info "gravity of the ice and water in $(round(time() - t0, digits = 1)) s"
        r = (; E = at.E[g0.sel], N = at.N[g0.sel], dg = g1.g .- g0.g)
        serialize(cache, r); r
    end
    # maps at 100 m
    xm, ym, tm = fetch_topography(R; dx = 0.1)
    BW = ice_water_thickness(R, tm; dx = 0.1)
    MG = valais_curvilinear_maps(R; dx = 0.25, glaciers = true, lakes = true)
    I = MG.interfaces
    ki = argmax(I.ice); kl = argmax(I.water)
    vol(a, dx) = sum(a) * dx^2
    fig = Figure(size = (1500, 1250))
    g1 = fig[1, 1] = GridLayout()
    ax = mapaxis(g1[1, 1]; title = "Ice thickness (Grab et al. 2021) and lake depth (swissBATHY3D) at 100 m")
    shade!(ax, xm, ym, tm)
    hi = heatmap!(ax, xm, ym, nanmask(1e3 .* BW.ice, BW.ice .<= 0.002), colormap = Reverse(:devon), colorrange = (0, 600), nan_color = :transparent)
    hw = heatmap!(ax, xm, ym, nanmask(1e3 .* BW.water, BW.water .<= 0.002), colormap = Reverse(:deep), colorrange = (0, 300), nan_color = :transparent)
    for (k, n) in ((ki, "A"), (kl, "B"))
        e = MG.x[Tuple(k)[1]]
        lines!(ax, [e, e], collect(R.N), color = :black, linewidth = 1.2, linestyle = :dash)
        text!(ax, e, R.N[2], text = "$(n)–$(n)'", align = (:center, :bottom), offset = (0, 2), fontsize = 13)
    end
    places!(ax); limits_region!(ax)
    Colorbar(g1[1, 2], hi, label = "ice thickness (m)", height = Relative(0.7))
    Colorbar(g1[1, 3], hw, label = "water depth (m)", height = Relative(0.7))
    # gravity effect, and what it does to the residual of the starting model of the density inversion
    g2 = fig[1, 2] = GridLayout()
    inner = (R.E[1] .<= GR.E .<= R.E[2]) .& (R.N[1] .<= GR.N .<= R.N[2])
    ax = mapaxis(g2[1, 1]; title = "Gravity effect of ice and water at the atlas points (ρ 917 and 1000 instead of 2670 / 2000 kg/m³)")
    shade!(ax, xm, ym, tm)
    contour!(ax, xm, ym, BW.ice, levels = [0.002], color = (:black, 0.5), linewidth = 0.5)
    sc = scatter!(ax, GR.E, GR.N, color = GR.dg, colormap = Reverse(:RdBu), colorrange = (-8, 8), markersize = 5,
                  strokecolor = (:black, 0.3), strokewidth = 0.2)
    limits_region!(ax)
    Colorbar(g2[1, 2], sc, label = "Δg (mGal)", height = Relative(0.7))
    note = "in the region: Δg from $(round(minimum(GR.dg[inner]), digits = 1)) to $(round(maximum(GR.dg[inner]), digits = 1)) mGal"
    dcache = datadir("output", "curvilinear_density_inversion.jls")
    if isfile(dcache)
        D = deserialize(dcache)
        if D.E == GR.E && D.N == GR.N
            rms(a) = sqrt(mean(abs2, a))
            r0 = D.dobs .- D.dg0 .- D.offset; r0 .-= mean(r0[inner])
            r1 = r0 .- GR.dg; r1 .-= mean(r1[inner])
            s = inner .& (GR.dg .< -1)
            note *= "\nresidual (observed − model) of the starting model of the density inversion: rms $(round(rms(r0[inner]), digits = 1)) → " *
                    "$(round(rms(r1[inner]), digits = 1)) mGal with ice and water.\nAt the $(count(s)) stations with Δg < −1 mGal (mean $(round(mean(GR.dg[s]), digits = 1)) mGal) " *
                    "the residual is already $(round(mean(r0[s]), digits = 1)) mGal (elsewhere $(round(mean(r0[inner .& .!s]), digits = 1))):\n" *
                    "the model is lighter than observed there, and adding the ice makes it worse. The atlas reduction may already treat\n" *
                    "glaciers and lakes, or the crystalline massifs under them are denser than in 3D-ALPS; the gravity figures keep them off."
        end
    end
    Label(g2[2, 1:2], note, fontsize = 12, halign = :left, tellwidth = false)
    # sections through the thickest ice and the deepest lake: every level and column of the 250 m curvilinear mesh
    for (c, k, n, nl, zl, ttl) in ((1, ki, "A", (1135.0, 1162.0), (0.8, 4.3), "glacier ice above the bed (Aletsch, Konkordiaplatz)"),
                                   (2, kl, "B", (1128.0, 1150.0), (-1.2, 1.0), "Lac Léman: water, Quaternary below the lake floor"))
        i = Tuple(k)[1]
        jj = findall(v -> nl[1] <= v <= nl[2], MG.y)
        F = valais_model_block(MG, i:i, jj)
        Zs = MG.z[i, jj, :]
        local ax = Axis(fig[2, c], title = "$(n)–$(n)' (E = $(MG.x[i]) km): $ttl", xlabel = "N (km)", ylabel = "elevation (km)")
        curvheatmap!(ax, MG.y[jj], Zs, Float64.(sq(F.Phases)), colormap = PH9CM, colorrange = (-0.5, 8.5))
        meshlines!(ax, MG.y[jj], Zs, MG.levels; every = 1, cols = 1, color = (:black, 0.2), linewidth = 0.3)
        lines!(ax, MG.y[jj], I.solid[i, jj], color = :black, linewidth = 1.2)
        xlims!(ax, nl...); ylims!(ax, zl...)
    end
    Legend(fig[3, 1:2], [PolyElement(color = PH9[k]) for k in (2, 3, 4, 8, 9)], PH9N[[2, 3, 4, 8, 9]], orientation = :horizontal, framevisible = false)
    rowsize!(fig.layout, 1, Relative(0.5))
    Label(fig[0, 1:2], "Glaciers and lakes in the model (glaciers = true, lakes = true): $(round(vol(I.ice, 0.25), digits = 1)) km³ of ice " *
                       "(up to $(round(Int, 1e3maximum(I.ice))) m) and $(round(vol(I.water, 0.25), digits = 1)) km³ of lake water " *
                       "(up to $(round(Int, 1e3maximum(I.water))) m) in the region; black line: the solid surface", fontsize = 14)
    savefig("glaciers", fig)
end

# ------------------------------------------------------------------------------------------------ Morcles
# the Morcles voxel model (50 m) sampled on a 2D curvilinear mesh (GMG curvilinear_grid with one y) that follows the topography
if "morcles" in WHICH
    m = read_morcles()
    n = length(MORCLES_UNITS)
    cm = cgrad([Makie.to_color("#9e6b5a"); [get(cgrad(:roma), (k - 1) / (n - 2)) for k in 1:n-1]], n, categorical = true)
    jc = argmin(abs.(m.y .- (m.y[1] + m.y[end]) / 2)); yc = m.y[jc]
    xt, yt, tt = fetch_topography(ValaisRegion((floor(m.x[1]), ceil(m.x[end])), (floor(yc) - 1, floor(yc) + 1), R.z); dx = 0.05)
    jt = col_index(yt, yc)
    topoline(e) = (i = clamp(searchsortedlast(xt, e), 1, length(xt) - 1); w = clamp((e - xt[i]) / (xt[i + 1] - xt[i]), 0, 1); (1 - w) * tt[i, jt] + w * tt[i + 1, jt])
    xs = collect(m.x[1]:0.05:m.x[end])
    zb = m.z[1]
    M = curvilinear_grid(xs, yc, (topoline, zb), (round(Int, (maximum(topoline.(xs)) - zb) / 0.025),))
    Z = M.z.val[:, 1, :]
    U = [Float64(morcles_unit_at(m, xs[i], yc, Z[i, k])) for i in eachindex(xs), k in axes(Z, 2)]
    U[U .== 0] .= NaN
    fig = Figure(size = (1500, 1000))
    ax1 = Axis(fig[1, 1], title = "Morcles voxel model (Thornton et al. 2018, 50 m): W–E section at N = $(round(yc, digits = 1)) km",
               ylabel = "elevation (km)", aspect = DataAspect())
    heatmap!(ax1, m.x, m.z, nanmask(m.unit[:, jc, :], m.unit[:, jc, :] .== 0), colormap = cm, colorrange = (0.5, n + 0.5), nan_color = :white)
    lines!(ax1, xs, topoline.(xs), color = :black, linewidth = 1)
    ax2 = Axis(fig[2, 1], title = "the same, sampled on a 2D curvilinear mesh (GMG curvilinear_grid, 50 m columns, $(size(Z, 2)) levels from the topography to $(zb) km)",
               xlabel = "E (km)", ylabel = "elevation (km)", aspect = DataAspect())
    curvheatmap!(ax2, xs, Z, U, colormap = cm, colorrange = (0.5, n + 0.5))
    lines!(ax2, xs, Z[:, end], color = :black, linewidth = 1)
    for ax in (ax1, ax2)
        xlims!(ax, extrema(m.x)...); ylims!(ax, zb, 3.3)
    end
    # zoom on a ridge: staircase of the voxels against the terrain-following mesh
    zx = (2572.0, 2575.0)
    for (c, (t, mesh)) in enumerate((("voxels", false), ("curvilinear mesh, every level", true)))
        ax = Axis(fig[3, 1][1, c], title = "zoom: $t", xlabel = "E (km)", ylabel = c == 1 ? "elevation (km)" : "", aspect = DataAspect())
        if mesh
            curvheatmap!(ax, xs, Z, U, colormap = cm, colorrange = (0.5, n + 0.5))
            ii = findall(e -> zx[1] - 0.1 <= e <= zx[2] + 0.1, xs)
            meshlines!(ax, xs[ii], Z[ii, :], (;); every = 1, color = (:black, 0.3), linewidth = 0.3, surfaces = false)
        else
            heatmap!(ax, m.x, m.z, nanmask(m.unit[:, jc, :], m.unit[:, jc, :] .== 0), colormap = cm, colorrange = (0.5, n + 0.5), nan_color = :white)
        end
        lines!(ax, xs, topoline.(xs), color = :black, linewidth = 1.5)
        ii = findall(e -> zx[1] <= e <= zx[2], xs); tz = topoline.(xs[ii])
        xlims!(ax, zx...); ylims!(ax, minimum(tz) - 0.6, maximum(tz) + 0.1)
    end
    Legend(fig[1:3, 2], [PolyElement(color = cm[k]) for k in 1:n], MORCLES_UNITS, "formation (old → young)", framevisible = false, rowgap = 0, labelsize = 11)
    savefig("morcles", fig)
end

# ------------------------------------------------------------------------------- showcase (one image, e.g. LinkedIn)
if "showcase" in WHICH
    W = 1600
    ink, muted = RGBf(0.09, 0.09, 0.09), RGBf(0.33, 0.33, 0.31)
    fence = render_fence("showcase_fence", (1400, 1150); every = 6)
    maps = MAPS
    B = section_B(maps)
    fig = Figure(size = (W, 3000), fonts = TEX, backgroundcolor = :white, figure_padding = (36, 36, 24, 24))
    Label(fig[1, 1:2], "GeoValaisWallis3D.jl", fontsize = 64, font = :bold, color = ink, halign = :left, tellwidth = false)
    Label(fig[2, 1:2], "Valais/Wallis in 3D on a curvilinear mesh of GeophysicalModelGenerator.jl: grid levels follow the topography, the sediments, the Conrad and the Moho",
          fontsize = 27, color = muted, halign = :left, justification = :left, word_wrap = true, tellwidth = false)
    cap(pos, s) = Label(pos, s, fontsize = 27, color = ink, halign = :left, justification = :left, tellwidth = false, word_wrap = true)

    cap(fig[3, 1:2], "The 3D model: tectonic map on the relief, phases and mesh down to 70 km")
    _, fsz = imgaxis!(fig[4, 1:2], fence)
    els, labs = surflegend((; base_cs = 0, conrad = 0, moho = 0))
    Legend(fig[5, 1:2], [[PolyElement(color = PHCOLS[k]) for k in 3:6]; els], [["sedimentary nappes", "upper crust", "lower crust", "mantle"]; labs],
           orientation = :horizontal, framevisible = false, labelsize = 24, colgap = 14, nbanks = 2, tellwidth = false)

    # zoom on the mesh in the Rhône valley
    cap(fig[6, 1], "Terrain-following mesh at Sion: 25 m levels hold the Quaternary fill")
    zx = (2586.0, 2602.0); ix = findall(e -> zx[1] <= e <= zx[2], B.h)
    axz = Axis(fig[7, 1], xlabel = "E (km)", ylabel = "z (km)", xlabelsize = 22, ylabelsize = 22, xticklabelsize = 20, yticklabelsize = 20)
    curvheatmap!(axz, B.h[ix], B.Z[ix, :], Float64.(sq(B.F.Phases)[ix, :]), colormap = PHCM, colorrange = (-0.5, 6.5))
    meshlines!(axz, B.h[ix], B.Z[ix, :], maps.levels; every = 1, linewidth = 0.3)
    xlims!(axz, zx...); ylims!(axz, -1.2, 1.6)
    # gravity inversion
    dres = DRES
    zq = -16.0
    inner = (R.E[1] .<= dres.E .<= R.E[2]) .& (R.N[1] .<= dres.N .<= R.N[2])
    rs = dres.dobs .- dres.dg0 .- dres.offset; rs .-= mean(rs[inner]); rinv = dres.dobs .- dres.dg
    rmsv(a) = round(sqrt(mean(abs2, a)), digits = 1)
    cap(fig[6, 2], "Gravity inversion on the mesh")
    g = fig[7, 2] = GridLayout()
    ax1 = Axis(g[1, 1], aspect = DataAspect(), titlesize = 24, titlefont = :regular, titlecolor = ink,
               title = "density change at $(round(Int, -zq)) km depth\nBouguer misfit $(rmsv(rs[inner])) → $(rmsv(rinv[inner])) mGal")
    hidedecorations!(ax1)
    hm = heatmap!(ax1, dres.x, dres.y, zslice(dres.z, dres.drho, zq), colormap = Reverse(:RdBu), colorrange = (-200, 200))
    lines!(ax1, [R.E[1], R.E[2], R.E[2], R.E[1], R.E[1]], [R.N[1], R.N[1], R.N[2], R.N[2], R.N[1]], color = ink, linestyle = :dash)
    for (n, e, nn) in PLACES[[3, 6]]
        scatter!(ax1, [e], [nn], color = ink, markersize = 9)
        text!(ax1, e, nn, text = n, fontsize = 22, offset = (6, 4), color = ink, glowcolor = (:white, 0.8), glowwidth = 3)
    end
    Colorbar(g[2, 1], hm, label = "Δρ (kg/m³)", vertical = false, flipaxis = false, width = Relative(0.8), labelsize = 22, ticklabelsize = 20)

    Label(fig[8, 1:2], "Data: swisstopo (swissALTIRegio, swissBEDROCK, Tectonic Map 1:500 000, Gravimetric Atlas), 3D-ALPS (Spooner et al. 2019), " *
                       "Spada et al. (2013).    github.com/tduretz/GeoValaisWallis3D.jl",
          fontsize = 18, color = muted, halign = :left, justification = :left, tellwidth = false, word_wrap = true)
    colsize!(fig.layout, 1, Relative(0.55))
    cw = W - 72
    rowsize!(fig.layout, 4, Fixed(cw * fsz[2] / fsz[1]))
    rowsize!(fig.layout, 7, Fixed(0.55cw * 0.62))
    rowgap!(fig.layout, 6)
    for r in (2, 5, 7); rowgap!(fig.layout, r, 18); end
    resize_to_layout!(fig)
    savefig("showcase", fig)
end
