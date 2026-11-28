"""
    GeoValaisWallis3D

Tools to build a 3D geological/geophysical model of the Valais/Wallis region (Swiss Alps)
from the surface down to the Moho, based on GeophysicalModelGenerator.jl.

Data sources
- Topography: swisstopo swissALTIRegio (10 m DEM, LV95), https://map.geo.admin.ch
- Surface tectonics: swisstopo Tectonic Map of Switzerland 1:500 000
- Quaternary deposits: swisstopo swissBEDROCK (bedrock elevation / thickness of unconsolidated deposits)
- Crustal layers & densities: 3D-ALPS (Spooner et al. 2019); temperature: 3D-ALPS-TR (Spooner et al. 2020)
- Moho: Spada et al. (2013, GJI) Western Alps Moho maps (European + Adriatic Moho)
- Nappe de Morcles 3D geological model (Thornton et al. 2018, Sci. Data)
- Gravity: swisstopo Gravimetric Atlas 1:100 000 (Bouguer anomalies), Gravimetric Map 1:500 000,
  National Gravity Network LSN2004 (absolute gravity)
- Seismic subsoil classes (SIA 261): BAFU national map, including the Valais 1:25 000 maps (CREALP)
- Seismicity: SED earthquake catalogue (FDSN; ECOS historical events); focal mechanisms of Truttmann et al. (2023)

The model is built in Swiss LV95 coordinates (EPSG:2056), expressed in km, on a rectilinear grid
(`build_valais_model`) or on a curvilinear mesh of GeophysicalModelGenerator whose levels follow the topography,
the base of the consolidated sediments, the Conrad and the Moho (`build_valais_curvilinear_model`).
"""
module GeoValaisWallis3D

using Reexport
@reexport using GeophysicalModelGenerator
import ArchGDAL as AG
using DelimitedFiles
using Downloads
using Serialization: serialize, deserialize
using LinearAlgebra, Statistics
using Krylov: lsqr
using LinearMaps: LinearMap

export ValaisRegion, valais_region, datadir
export lv95_to_lonlat, lonlat_to_lv95
export fetch_topography, fetch_swissbedrock
export fetch_seismicity, fetch_focal_mechanisms, read_focal_mechanisms, read_relocations_truttmann2023, download_truttmann2023
export fault_vectors, pt_axes, beachball
export download_glacier_beds, fetch_glacier_bed, fetch_lake_bathymetry, ice_water_thickness
export download_tectonic_map, tectonic_map_on_grid, TectonicMap
export download_moho_spada, read_moho_spada, interpolate_scattered
export download_3dalps, download_3dalps_tr, alps3d_on_grid, temperature_3dalps
export read_morcles, MorclesModel, MORCLES_UNITS, morcles_unit_at
export refined_levels, fine_band_from_bedrock, drop_thin_deposits
export LayerVelocity, DEFAULT_VELOCITIES
export build_valais_model, valais_model_maps, valais_model_block, model_interfaces, PHASES, cartdata_to_geodata
export fetch_gravity_atlas, fetch_bouguer_500k, fetch_gravity_network
export topographic_attraction, gravity_maps, write_gravity_paraview, model_gravity, upward_continue
export density_on_grid, bouguer_contrast, invert_density, interpolate_nodes
export SUBSOIL_CLASSES, download_subsoil_classes, subsoil_release, subsoil_classes_on_grid, subsoil_map, write_subsoil_paraview
export download_spectral_microzonation, spectral_microzonation_on_grid
export CURVILINEAR_SURFACES, meshing_surfaces, curvilinear_levels, valais_curvilinear_maps, build_valais_curvilinear_model
export curvilinear_density_model, curvilinear_bouguer_contrast, curvilinear_heights, layer_cells, interpolate_curvilinear

"""
    datadir(parts...)

Path inside the data directory (downloads go to `data/raw`, outputs to `data/output`). This is the
package's `data/` folder, or `\$GEOVALAIS_DATADIR` if that environment variable is set (e.g. a scratch
file system on a cluster).
"""
datadir(parts...) = joinpath(get(ENV, "GEOVALAIS_DATADIR", joinpath(dirname(@__DIR__), "data")), parts...)

include("region.jl")
include("topography.jl")
include("bedrock.jl")
include("glaciers_lakes.jl")
include("tectonics.jl")
include("moho.jl")
include("alps3d.jl")
include("morcles.jl")
include("grid.jl")
include("velocity.jl")
include("model.jl")
include("gravity.jl")
include("subsoil.jl")
include("curvilinear.jl")
include("seismicity.jl")

end # module GeoValaisWallis3D
