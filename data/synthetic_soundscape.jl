# synthetic_soundscape.jl
# Synthetic time–frequency toy soundscape for explaining NMF.
#
# Convention: matrices are (F × T) = (frequency bins × time frames)
using Random

"""
    toy_spectrogram(F, T; kwargs...) -> V, parts

Construct a toy spectrogram `V(f,t)` directly in time–frequency space as a sum of
simple, interpretable components:

    V = V_whistle + V_clicks + V_ship + V_bg

Math (conceptual):
- Whistle: narrowband ridge around a smooth contour fw(t)
    V_whistle(f,t) = Aw * exp(-(f - fw(t))^2 / (2σw^2))
  where fw(t) = f0 + Δf*sin(2π t / Tmod)

- Clicks: sparse broadband vertical bursts at times tᵢ
    V_clicks(f,t) = Ac * Σᵢ exp(-(t - tᵢ)^2 / (2σt^2)) * g_click(f)
  where g_click(f) is a (nearly) flat spectral shape.

- Ship: persistent low-frequency band

- Background: weak diffuse noise floor
    V_bg(f,t) = Ab * noise(f,t)

Ecological / signal rationale:
These components mimic common soundscape morphologies (tonal whistles, impulsive
clicks, persistent low-frequency shipping, and diffuse background) to make the
NMF idea visually intuitive. This is not a physical audio simulation.

Returns:
- `V::Matrix{Float64}`: total toy spectrogram (F×T)
- `parts::NamedTuple`: individual component matrices for inspection/tuning
"""
function toy_spectrogram(
    F::Int, T::Int;
    # Whistle params (in bin units)
    Aw::Float64=1.0, f0::Float64=0.55F, Δf::Float64=0.15F, Tmod::Float64=0.50T, σw::Float64=6.0,
    # Click params
    Ac::Float64=0.9, nclicks::Int=10, σt::Float64=1.5,
    # Ship params
    As::Float64=0.8, σs::Float64=35.0, ship_mod_amp::Float64=0.15, ship_mod_T::Float64=0.80T,
    # Background params
    Ab::Float64=0.06, rng_seed::Int=1
)
    # Axes (index space)
    f = collect(1:F)
    t = collect(1:T)

    # --- Whistle: ridge along a smooth contour ---
    fw = f0 .+ Δf .* sin.(2π .* t ./ Tmod)              # length T
    V_whistle = Aw .* exp.(-((f .- fw').^2) ./ (2σw^2)) # F×T

    # --- Clicks: sparse vertical bursts (broadband) ---
    # choose click times away from edges for nicer visuals

    rng = MersenneTwister(rng_seed)
    click_times = sort(rand(rng, 10:(T-10), nclicks))

    # time-envelope for clicks (1×T), sum of Gaussians in time
    env_t = zeros(Float64, T)
    for ti in click_times
        env_t .+= exp.(-((t .- ti).^2) ./ (2σt^2))
    end
    env_t ./= maximum(env_t) > 0 ? maximum(env_t) : 1.0

    # nearly-flat spectral shape for clicks (F×1)
    g_click = ones(Float64, F)

    V_clicks = Ac .* (g_click .* env_t')               # F×T

    # --- Ship: low-frequency band, persistent with slow modulation ---
    g_ship = exp.(-(f.^2) ./ (2σs^2))                  # F
    ship_env = 1 .+ ship_mod_amp .* sin.(2π .* t ./ ship_mod_T)  # T
    V_ship = As .* (g_ship .* ship_env')               # F×T

    # --- Background: weak diffuse noise floor ---
    V_bg = Ab .* rand(rng, Float64, F, T)

    V = V_whistle .+ V_clicks .+ V_ship .+ V_bg

    parts = (whistle=V_whistle, clicks=V_clicks, ship=V_ship, bg=V_bg,
             fw=fw, click_times=click_times)

    return V, parts
end


using CairoMakie

fig = Figure()
ax  = Axis(fig[1,1], xlabel="time", ylabel="frequency")
heatmap!(ax, V; interpolate=false)
fig
