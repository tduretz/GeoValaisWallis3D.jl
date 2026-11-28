"""
    LayerVelocity(; Vp0, dVp = 0, Vs0, dVs = 0, Vpmax = Inf, Vsmax = Inf, rho = NaN)

Seismic velocities of one model layer as a linear function of the depth `d` (km) below
the topographic surface: `Vp = min(Vp0 + dVp*d, Vpmax)`, `Vs = min(Vs0 + dVs*d, Vsmax)` (km/s).
`rho` (kg/m³) overrides the 3D-ALPS density when not `NaN`.
"""
Base.@kwdef struct LayerVelocity
    Vp0::Float64
    dVp::Float64 = 0.0
    Vs0::Float64
    dVs::Float64 = 0.0
    Vpmax::Float64 = Inf
    Vsmax::Float64 = Inf
    rho::Float64 = NaN
end

velocity(l::LayerVelocity, d) = (min(l.Vp0 + l.dVp * d, l.Vpmax), min(l.Vs0 + l.dVs * d, l.Vsmax))

"""
Default (first-order, to be calibrated) velocity laws per phase.

- Quaternary fill: Vs increasing from ~0.35 km/s at the surface to ~0.8 km/s at 500 m depth,
  in the range reported for the Rhône valley (e.g. Roten et al. 2006, 2008; Frischknecht & Wagner 2004).
- Consolidated (Alpine sedimentary nappes) / upper crust / lower crust / mantle: typical values of
  the western Alpine crustal model (Wagner et al. 2012; Diehl et al. 2009).
- Adriatic mantle wedge (Ivrea body): Vp ≈ 7.4 km/s.
- Glacier ice: Vp 3.8, Vs 1.9 km/s, 917 kg/m³ (typical of temperate glacier ice, Vp ≈ 3.6–3.9 km/s). Lake water: Vp 1.48 km/s, Vs 0 (a fluid: check that the wave solver accepts it), 1000 kg/m³.
"""
const DEFAULT_VELOCITIES = (
    quaternary             = LayerVelocity(Vp0 = 1.60, dVp = 1.6, Vs0 = 0.35, dVs = 0.9, Vpmax = 2.6, Vsmax = 1.1, rho = 2000.0),
    consolidated_sediments = LayerVelocity(Vp0 = 4.80, dVp = 0.08, Vs0 = 2.70, dVs = 0.045, Vpmax = 5.9, Vsmax = 3.4),
    upper_crust            = LayerVelocity(Vp0 = 5.80, dVp = 0.02, Vs0 = 3.35, dVs = 0.012, Vpmax = 6.3, Vsmax = 3.65),
    lower_crust            = LayerVelocity(Vp0 = 6.50, dVp = 0.01, Vs0 = 3.70, dVs = 0.006, Vpmax = 7.0, Vsmax = 4.0),
    mantle                 = LayerVelocity(Vp0 = 8.00, dVp = 0.003, Vs0 = 4.50, dVs = 0.002),
    adriatic_mantle        = LayerVelocity(Vp0 = 7.40, dVp = 0.02, Vs0 = 4.20, dVs = 0.01, Vpmax = 8.0, Vsmax = 4.5),
    ice                    = LayerVelocity(Vp0 = 3.80, Vs0 = 1.90, rho = 917.0),
    water                  = LayerVelocity(Vp0 = 1.48, Vs0 = 0.0, rho = 1000.0),
)
