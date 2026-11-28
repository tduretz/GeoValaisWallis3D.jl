# Seismicity: the earthquake catalogue of the Swiss Seismological Service (SED) and focal mechanisms.

# SED FDSN event service: the instrumental catalogue (from 1975) and the historical events of ECOS (from 250 AD)
const SED_FDSN_EVENT = "http://eida.ethz.ch/fdsnws/event/1/query"
# Truttmann, Diehl & Herwegh (2023), data of "Hypocenter-based 3D imaging of active faults" (JGR 128, e2023JB026352):
# focal mechanisms and double-difference relocations of the St-Léonard and Anzère sequences; Zenodo, CC BY 4.0
const TRUTTMANN2023_URL = "https://zenodo.org/api/records/7855039/files/Truttmann_etal_2023.zip/content"
const TRUTTMANN2023_MD5 = "1ace08887161daec1660dc3602634bf3"

# lon/lat bounding box of the region extended by `margin` km
function _lonlat_box(r, margin)
    E = (r.E[1] - margin, r.E[2] + margin); N = (r.N[1] - margin, r.N[2] + margin)
    lon, lat = lv95_to_lonlat([E[1], E[2], E[1], E[2]], [N[1], N[1], N[2], N[2]])
    return extrema(lon), extrema(lat)
end
_inside(E, N, r, margin) = r.E[1] - margin <= E <= r.E[2] + margin && r.N[1] - margin <= N <= r.N[2] + margin

"""
    fetch_seismicity(r = valais_region(); margin = 0, eventtype = "earthquake", minmagnitude = -Inf, cache = true)
        -> NamedTuple of vectors (time, E, N, z, depth, mag, magtype, name, type, id)

The earthquake catalogue of the **Swiss Seismological Service** (SED, FDSN event service `eida.ethz.ch`) in the region
extended by `margin` km: the instrumental catalogue since 1975 and the historical events of ECOS (ECOS-09, with moment
magnitudes Mw). `E`, `N` are LV95 km; `depth` (km) as reported by the SED, `z = -depth` (NaN where the SED gives no depth,
e.g. many historical events); `mag` with its type `magtype` (MLhc, MLh, Ml: local magnitudes; Mw); `time` an ISO string.

`eventtype` keeps one SED event type (`"earthquake"`; `nothing` keeps all, including quarry blasts, explosions,
landslides and ice quakes). The whole catalogue of the box (all types and magnitudes) is downloaded once and cached in
`data/raw/sed_catalogue_<region>.txt` (delete it to update).
"""
function fetch_seismicity(r::ValaisRegion = valais_region(); margin = 0, eventtype = "earthquake", minmagnitude = -Inf,
                          cache = true)
    mkpath(datadir("raw"))
    file = datadir("raw", "sed_catalogue_$(_region_tag(r, 1.0))_m$(round(Int, margin)).txt")
    if !(cache && isfile(file))
        (lo, la) = _lonlat_box(r, margin + 1.0)
        q = "?starttime=0001-01-01&minlatitude=$(round(la[1], digits = 4))&maxlatitude=$(round(la[2], digits = 4))" *
            "&minlongitude=$(round(lo[1], digits = 4))&maxlongitude=$(round(lo[2], digits = 4))&format=text"
        @info "Downloading the SED earthquake catalogue (FDSN event service) ..."
        txt = String(take!(Downloads.download(SED_FDSN_EVENT * q, IOBuffer())))
        startswith(txt, "#EventID") || error("fetch_seismicity: unexpected answer of $(SED_FDSN_EVENT): $(first(txt, 200))")
        write(file, txt)
        @info "Saved $(file) ($(count(==('\n'), txt) - 1) events)"
    end
    return _parse_sed_text(readlines(file), r, margin; eventtype, minmagnitude)
end

function _parse_sed_text(lines, r, margin; eventtype = "earthquake", minmagnitude = -Inf)
    time = String[]; lon = Float64[]; lat = Float64[]; depth = Float64[]; mag = Float64[]
    magtype = String[]; name = String[]; type = String[]; id = String[]
    num(s) = (v = tryparse(Float64, s); v === nothing ? NaN : v)
    for l in lines
        (isempty(l) || startswith(l, "#")) && continue
        c = split(l, '|')
        length(c) < 14 && continue
        t = String(c[14])
        (eventtype === nothing || t == eventtype) || continue
        m = num(c[11])
        (isnan(m) ? -Inf : m) >= minmagnitude || continue
        push!(id, c[1]); push!(time, c[2]); push!(lat, num(c[3])); push!(lon, num(c[4])); push!(depth, num(c[5]))
        push!(magtype, c[10]); push!(mag, m); push!(name, strip(c[13])); push!(type, t)
    end
    E, N = isempty(lon) ? (Float64[], Float64[]) : lonlat_to_lv95(lon, lat)
    keep = [_inside(E[k], N[k], r, margin) for k in eachindex(E)]
    s(v) = v[keep]
    return (; time = s(time), E = s(E), N = s(N), z = s(-depth), depth = s(depth), mag = s(mag), magtype = s(magtype),
            name = s(name), type = s(type), id = s(id))
end

"""
    download_truttmann2023(; force = false) -> directory

Download the data of Truttmann, Diehl & Herwegh (2023, JGR Solid Earth, doi:10.1029/2023JB026352; Zenodo
doi:10.5281/zenodo.7855039, CC BY 4.0, 35 MB) to `data/raw/truttmann2023/`: focal mechanisms and double-difference
relocations of the St-Léonard (2014–2016) and Anzère (2019) earthquake sequences north of Sion.
"""
function download_truttmann2023(; force = false)
    dir = datadir("raw", "truttmann2023")
    out = joinpath(dir, "JGR_data_Zenodo")
    (isdir(out) && !force) && return out
    mkpath(dir)
    zip = joinpath(dir, "Truttmann_etal_2023.zip")
    if force || !isfile(zip) || _md5(zip) != TRUTTMANN2023_MD5
        @info "Downloading Truttmann et al. (2023) from Zenodo (35 MB) ..."
        Downloads.download(TRUTTMANN2023_URL, zip)
        _md5(zip) == TRUTTMANN2023_MD5 || (rm(zip; force = true); error("download_truttmann2023: MD5 mismatch"))
    end
    run(`unzip -q -o $zip -x "__MACOSX/*" -d $dir`)
    return out
end

"""
    read_focal_mechanisms(files...) -> NamedTuple of vectors

Focal mechanisms in the format of the SED compilations (semicolon-separated, a header line `Yr;Mo;Dy;Hr:Mi;Lat;Lon;Z;
Mag;A;Strike1;Dip1;Rake1;Strike2;Dip2;Rake2;Pazim;Pdip;Tazim;Tdip;Q;Type;Loc;Lit`), as published by Truttmann et al.
(2023). Returns `time` (ISO, to the minute), `E`, `N` (LV95 km), `z = -depth` (km), `mag` (ML), `active` (the active
nodal plane, 1 or 2; 0 unknown), the nodal planes `strike1, dip1, rake1, strike2, dip2, rake2`, the P and T axes
(`p_azimuth, p_plunge, t_azimuth, t_plunge`, degrees), `quality`, `style` (e.g. SS strike-slip, NF normal, TF thrust),
`locality` and `reference`.
"""
function read_focal_mechanisms(files::AbstractString...)
    T = (time = String[], lon = Float64[], lat = Float64[], z = Float64[], mag = Float64[], active = Int[],
         strike1 = Float64[], dip1 = Float64[], rake1 = Float64[], strike2 = Float64[], dip2 = Float64[], rake2 = Float64[],
         p_azimuth = Float64[], p_plunge = Float64[], t_azimuth = Float64[], t_plunge = Float64[], quality = Int[],
         style = String[], locality = String[], reference = String[])
    for f in files
        lines = filter(l -> !isempty(strip(l)) && !startswith(l, "#"), readlines(f))
        h = strip.(split(lines[1], ';'))
        nh = count(!isempty, h)
        ix(n) = findfirst(==(n), h)
        hm = ix("Hr:Mi")
        for l in lines[2:end]
            c = strip.(split(l, ';'))
            v(n) = parse(Float64, c[ix(n)])
            hh, mm = split(c[hm], ':')
            push!(T.time, "$(c[ix("Yr")])-$(lpad(c[ix("Mo")], 2, '0'))-$(lpad(c[ix("Dy")], 2, '0'))T$(lpad(hh, 2, '0')):$(lpad(mm, 2, '0'))")
            push!(T.lat, v("Lat")); push!(T.lon, v("Lon")); push!(T.z, -v("Z")); push!(T.mag, v("Mag"))
            push!(T.active, parse(Int, c[ix("A")]))
            for (k, n) in ((:strike1, "Strike1"), (:dip1, "Dip1"), (:rake1, "Rake1"), (:strike2, "Strike2"), (:dip2, "Dip2"),
                           (:rake2, "Rake2"), (:p_azimuth, "Pazim"), (:p_plunge, "Pdip"), (:t_azimuth, "Tazim"), (:t_plunge, "Tdip"))
                push!(T[k], v(n))
            end
            push!(T.quality, parse(Int, c[ix("Q")])); push!(T.style, c[ix("Type")])
            # the Truttmann files have one more column than their header (a sequence code before the locality)
            extra = count(!isempty, c) > nh
            push!(T.locality, strip(c[ix("Loc") + extra])); push!(T.reference, strip(c[ix("Lit") + extra]))
        end
    end
    E, N = isempty(T.lon) ? (Float64[], Float64[]) : lonlat_to_lv95(T.lon, T.lat)
    return merge((; E, N), T)
end

"""
    fetch_focal_mechanisms(r = valais_region(); files = nothing, margin = 0) -> NamedTuple of vectors

The focal mechanisms in the region (see [`read_focal_mechanisms`](@ref)): by default those of Truttmann et al. (2023)
(St-Léonard and Anzère sequences, [`download_truttmann2023`](@ref)); `files` reads other files in the SED format instead.
"""
function fetch_focal_mechanisms(r::ValaisRegion = valais_region(); files = nothing, margin = 0)
    if files === nothing
        d = joinpath(download_truttmann2023(), "catalogs_focalmechanisms")
        files = [joinpath(d, f) for f in readdir(d) if startswith(f, "FocalMechanisms_") && !occursin("HEADER", f)]
    end
    F = read_focal_mechanisms(files...)
    keep = [_inside(F.E[k], F.N[k], r, margin) for k in eachindex(F.E)]
    return map(v -> v[keep], F)
end

"""
    read_relocations_truttmann2023() -> NamedTuple of vectors (time, E, N, z, mag, cluster, sequence)

The double-difference relocations (hypoDD) of the St-Léonard and Anzère sequences of Truttmann et al. (2023).
"""
function read_relocations_truttmann2023()
    d = joinpath(download_truttmann2023(), "catalogs_hypoDD")
    time = String[]; lon = Float64[]; lat = Float64[]; z = Float64[]; mag = Float64[]; cluster = Int[]; sequence = String[]
    for (f, s) in (("hypoDD_StLeonard.txt", "St-Léonard"), ("hypoDD_Anzere.txt", "Anzère"))
        A = readdlm(joinpath(d, f), '\t', Any; skipstart = 1)
        for k in axes(A, 1)
            push!(lat, A[k, 2]); push!(lon, A[k, 3]); push!(z, -A[k, 4]); push!(mag, A[k, 17]); push!(cluster, A[k, 24])
            push!(time, "$(A[k, 11])-$(lpad(A[k, 12], 2, '0'))-$(lpad(A[k, 13], 2, '0'))T$(lpad(A[k, 14], 2, '0')):$(lpad(A[k, 15], 2, '0'))")
            push!(sequence, s)
        end
    end
    E, N = lonlat_to_lv95(lon, lat)
    return (; time, E, N, z, mag, cluster, sequence)
end

"""
    fault_vectors(strike, dip, rake) -> (n, u)

Unit normal `n` of the fault plane and slip direction `u` (of the hanging wall) in (north, east, down) coordinates, from
strike, dip and rake in degrees (Aki & Richards convention).
"""
function fault_vectors(strike, dip, rake)
    s, d, r = deg2rad(strike), deg2rad(dip), deg2rad(rake)
    n = (-sin(d) * sin(s), sin(d) * cos(s), -cos(d))
    u = (cos(r) * cos(s) + cos(d) * sin(r) * sin(s), cos(r) * sin(s) - cos(d) * sin(r) * cos(s), -sin(r) * sin(d))
    return n, u
end

"""
    pt_axes(strike, dip, rake) -> ((p_azimuth, p_plunge), (t_azimuth, t_plunge))

Azimuth and plunge (degrees) of the P (pressure) and T (tension) axes of a double couple.
"""
function pt_axes(strike, dip, rake)
    n, u = fault_vectors(strike, dip, rake)
    function azpl(v)
        v = v[3] < 0 ? .-v : v                         # pointing down
        (mod(rad2deg(atan(v[2], v[1])), 360.0), rad2deg(asin(clamp(v[3] / sqrt(sum(abs2, v)), -1, 1))))
    end
    return azpl(n .- u), azpl(n .+ u)
end

"""
    beachball(strike, dip, rake; n = 61) -> (x, y, polarity)

Lower-hemisphere, equal-area projection of the P-wave polarities of a double couple on an `n × n` grid of the unit disc
(`x` east, `y` north): `polarity` is 1 for compression, 0 for dilatation and NaN outside the disc. For drawing focal
mechanisms (e.g. `heatmap!(e .+ s .* x, n .+ s .* y, polarity)` at the epicentre `e`, `n`, with size `s`).
"""
function beachball(strike, dip, rake; n = 61)
    nv, u = fault_vectors(strike, dip, rake)
    x = range(-1, 1, length = n); y = range(-1, 1, length = n)
    pol = fill(NaN, n, n)
    for j in 1:n, i in 1:n
        ρ = hypot(x[i], y[j]); ρ > 1 && continue
        inc = 2asin(ρ / sqrt(2))                        # take-off angle from the vertical (down)
        φ = atan(x[i], y[j])                            # azimuth from north
        γ = (sin(inc) * cos(φ), sin(inc) * sin(φ), cos(inc))
        pol[i, j] = sum(u .* γ) * sum(nv .* γ) >= 0 ? 1.0 : 0.0
    end
    return collect(x), collect(y), pol
end
