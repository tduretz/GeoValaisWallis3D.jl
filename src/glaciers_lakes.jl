# Glaciers (ice thickness) and lakes (water depth): the bodies above the solid surface that the topography includes.
#
# swissALTIRegio is the surface of the ice and of the lakes, and swissBEDROCK counts the Quaternary thickness from the
# glacier bed and the lake floor. Without these two data sets the model fills glaciers with rock and starts the Quaternary
# at the lake surface (so the bedrock under a lake is too high by the water depth).

# Swiss Glacier Thickness – Release 2020 (Grab et al. 2021, J. Glaciol.; ETH Research Collection, CC BY 4.0):
# glacier bed elevation of the glaciers of the Swiss Glacier Inventory, one 10 m LV95 GeoTIFF per glacier.
const GLACIER_BED_URL = "https://www.research-collection.ethz.ch/bitstreams/91676987-d124-44b8-98d0-89b942174252/download"
const GLACIER_BED_MD5 = "d887deda439979c6c4d1ac54f5202070"
const GLACIER_BED_ZIP = "07_GlacierBed_SGI.zip"
# swissBATHY3D (swisstopo, OGD): lake-floor elevation of the larger Swiss lakes, 1 km tiles of ESRI ASCII grids
const SWISSBATHY3D_ITEMS = "https://data.geo.admin.ch/api/stac/v1/collections/ch.swisstopo.swissbathy3d/items?limit=100"

"""
    download_glacier_beds(; force = false) -> directory

Download and unpack the glacier bed elevations of *Swiss Glacier Thickness – Release 2020* (Grab, M. et al., ETH Zurich,
doi:10.3929/ethz-b-000434697, CC BY 4.0; 122 MB) into `data/raw/glaciers/`. The archive is checked against its MD5
sum: the ETH server answers too many requests with an HTML page, which is detected and reported.
"""
function download_glacier_beds(; force = false)
    dir = datadir("raw", "glaciers")
    out = joinpath(dir, "07_GlacierBed_SGI")
    (isdir(out) && !force) && return out
    mkpath(dir)
    zip = joinpath(dir, GLACIER_BED_ZIP)
    if force || !isfile(zip) || _md5(zip) != GLACIER_BED_MD5
        @info "Downloading the glacier beds of Swiss Glacier Thickness R2020 (122 MB) ..."
        Downloads.download(GLACIER_BED_URL, zip)
        if _md5(zip) != GLACIER_BED_MD5
            msg = filesize(zip) < 100_000 && occursin("Rate Limited", read(zip, String)) ?
                  "the ETH Research Collection rate-limited the download; try again in a few minutes" :
                  "MD5 mismatch for $(zip)"
            rm(zip; force = true)
            error("download_glacier_beds: " * msg)
        end
    end
    run(`unzip -q -o $zip -d $dir`)
    return out
end

_md5(file) = first(split(readchomp(Sys.isapple() ? `md5 -q $file` : `md5sum $file`)))

# Accumulate the valid pixels of an LV95 raster (metres) onto the cells of the grid with lower-left corner (E0, N0) km,
# spacing dx km and nx × ny cells: S += value (m), C += pixel area / cell area.
function _accumulate_raster!(S, C, ds, E0, N0, dx)
    gt = AG.getgeotransform(ds)
    b = AG.getband(ds, 1)
    nd = AG.getnodatavalue(b)
    A = AG.read(b)
    w = abs(gt[2] * gt[6]) / (dx * 1e3)^2
    nx, ny = size(S)
    @inbounds for jp in axes(A, 2), ip in axes(A, 1)
        v = A[ip, jp]
        (isnothing(nd) || v != nd) && isfinite(v) || continue
        E = (gt[1] + (ip - 0.5) * gt[2]) / 1e3
        N = (gt[4] + (jp - 0.5) * gt[6]) / 1e3
        i = floor(Int, (E - E0) / dx) + 1
        j = floor(Int, (N - N0) / dx) + 1
        (1 <= i <= nx && 1 <= j <= ny) || continue
        S[i, j] += v * w; C[i, j] += w
    end
    return nothing
end

# (x, y) cell centres, mean value over the covered part of each cell (m → km; NaN where uncovered) and covered fraction
function _cover_grid(r, dx, files; open = AG.read)
    nx, ny = round(Int, (r.E[2] - r.E[1]) / dx), round(Int, (r.N[2] - r.N[1]) / dx)
    S = zeros(nx, ny); C = zeros(nx, ny)
    for f in files
        open(f) do ds
            _accumulate_raster!(S, C, ds, r.E[1], r.N[1], dx)
        end
    end
    x = [r.E[1] + (i - 0.5) * dx for i in 1:nx]
    y = [r.N[1] + (j - 0.5) * dx for j in 1:ny]
    mean = [c > 0 ? s / c / 1e3 : NaN for (s, c) in zip(S, C)]
    return x, y, mean, min.(C, 1.0)
end

# raster bounding box (km) without reading the data
function _raster_bbox(ds)
    gt = AG.getgeotransform(ds)
    E = (gt[1], gt[1] + AG.width(ds) * gt[2]) ./ 1e3
    N = (gt[4] + AG.height(ds) * gt[6], gt[4]) ./ 1e3
    return E, N
end

_overlaps(E, N, r) = E[1] < r.E[2] && E[2] > r.E[1] && N[1] < r.N[2] && N[2] > r.N[1]

function _write_cover(file, x, y, mean, frac)
    open(file, "w") do io
        serialize(io, (; x, y, mean, frac))
    end
end
_read_cover(file) = (c = open(deserialize, file); (c.x, c.y, c.mean, c.frac))

"""
    fetch_glacier_bed(r = valais_region(); dx = 0.25, cache = true) -> (x, y, bed, fraction)

Glacier bed elevation (km a.s.l.) of *Swiss Glacier Thickness – Release 2020* (Grab et al. 2021) on the grid of
[`fetch_topography`](@ref) (cell centres `x`, `y`, spacing `dx`): the mean over the glacier pixels of each cell (NaN
without glacier), and the fraction of the cell covered by glaciers. Cached in `data/raw/glacierbed_<region>.jls`.
"""
function fetch_glacier_bed(r::ValaisRegion = valais_region(); dx = 0.25, cache = true)
    file = datadir("raw", "glacierbed_$(_region_tag(r, dx)).jls")
    cache && isfile(file) && return _read_cover(file)
    dir = download_glacier_beds()
    files = filter(f -> endswith(f, ".tif") && AG.read(ds -> _overlaps(_raster_bbox(ds)..., r), f),
                   readdir(dir; join = true))
    @info "Glacier beds: $(length(files)) glaciers in the region"
    x, y, bed, frac = _cover_grid(r, dx, files)
    _write_cover(file, x, y, bed, frac)
    return x, y, bed, frac
end

"""
    fetch_lake_bathymetry(r = valais_region(); dx = 0.25, cache = true) -> (x, y, floor, fraction)

Lake-floor elevation (km a.s.l.) of swisstopo **swissBATHY3D** on the grid of [`fetch_topography`](@ref): the mean over
the lake pixels of each cell (NaN without lake) and the fraction of the cell covered by a lake. The lakes whose extent
meets the region are downloaded (in Valais: Lac Léman, 400 MB) to `data/raw/swissbathy3d/`; only their tiles in the
region are read. Cached in `data/raw/lakefloor_<region>.jls`.
"""
function fetch_lake_bathymetry(r::ValaisRegion = valais_region(); dx = 0.25, cache = true)
    file = datadir("raw", "lakefloor_$(_region_tag(r, dx)).jls")
    cache && isfile(file) && return _read_cover(file)
    dir = datadir("raw", "swissbathy3d"); mkpath(dir)
    # the STAC items (one per lake) whose extent meets the region: their ESRI ASCII grid archives
    lon, lat = lv95_to_lonlat([r.E[1], r.E[2], r.E[1], r.E[2]], [r.N[1], r.N[1], r.N[2], r.N[2]])
    url = SWISSBATHY3D_ITEMS * "&bbox=" * join(round.((minimum(lon), minimum(lat), maximum(lon), maximum(lat)), digits = 5), ",")
    json = String(take!(Downloads.download(url, IOBuffer())))
    hrefs = unique(String(m[1]) for m in eachmatch(r"\"href\"\s*:\s*\"([^\"]+\.esriasciigrid\.zip)\"", json))
    tiles = String[]
    for href in hrefs
        zip = joinpath(dir, basename(href))
        if !isfile(zip)
            @info "Downloading $(basename(href)) (swissBATHY3D) ..."
            Downloads.download(href, zip * ".part"); mv(zip * ".part", zip)
        end
        # tiles named ..._<E km>_<N km>.asc by their lower-left corner (1 km)
        names = String.(split(readchomp(`unzip -Z1 $zip`)))
        want = filter(names) do n
            m = match(r"_(\d{4})_(\d{4})\.asc$", n)
            m !== nothing && _overlaps((parse(Int, m[1]), parse(Int, m[1]) + 1), (parse(Int, m[2]), parse(Int, m[2]) + 1), r)
        end
        isempty(want) && continue
        tdir = joinpath(dir, splitext(basename(zip))[1]); mkpath(tdir)
        missing_tiles = filter(n -> !isfile(joinpath(tdir, n)), want)
        isempty(missing_tiles) || run(`unzip -q -o $zip $missing_tiles -d $tdir`)
        append!(tiles, joinpath.(tdir, want))
    end
    @info "swissBATHY3D: $(length(tiles)) lake tiles in the region"
    x, y, fl, frac = _cover_grid(r, dx, tiles)
    _write_cover(file, x, y, fl, frac)
    return x, y, fl, frac
end

"""
    ice_water_thickness(r, topo; dx, glaciers = true, lakes = true) -> (; ice, water)

Ice thickness and water depth (km) on the grid of `fetch_topography(r; dx)` with surface `topo` (km): in each cell the
covered fraction times the height of the surface above the glacier bed or the lake floor (see [`fetch_glacier_bed`](@ref),
[`fetch_lake_bathymetry`](@ref)), so that the mean over a partly covered cell is kept. Where the surface lies below the
glacier bed (glaciers that have melted since the inventory) the ice thickness is 0. Zeros for a body that is not asked for.
"""
function ice_water_thickness(r::ValaisRegion, topo; dx, glaciers = true, lakes = true)
    function body(f)
        _, _, b, c = f(r; dx)
        size(b) == size(topo) || throw(DimensionMismatch("topo is not on the grid of fetch_topography(r; dx = $dx)"))
        return [c[k] > 0 && !isnan(topo[k]) ? c[k] * max(topo[k] - b[k], 0.0) : 0.0 for k in eachindex(topo)]
    end
    ice = glaciers ? reshape(body(fetch_glacier_bed), size(topo)) : zeros(size(topo))
    water = lakes ? reshape(body(fetch_lake_bathymetry), size(topo)) : zeros(size(topo))
    return (; ice, water)
end
