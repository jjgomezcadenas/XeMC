"""
Event-level processing: multi-gamma events through fast kernel + FV stack.

Orchestrates the two-stage transport:
1. Fast kernel (`transport_gamma_fastkernel` from tracking_fast.jl)
   for rapid veto rejection outside the FV
2. Full stack (`propagate_gamma` from tracking.jl) for detailed
   deposit collection inside the FV

Returns raw FV deposits for offline SS/MS classification.
"""


"""
    EventProcessingResult

Result of processing a multi-gamma event.

`status` values:
- `:fv` — at least one gamma reached the FV; `deposits` contains the
  full deposit stack for offline SS/MS classification
- `:vetoed` — a deposit in Skin or TPC active exceeded the veto threshold
- `:no_fv` — no gamma reached the FV
"""
struct EventProcessingResult
    status::Symbol          # :fv, :vetoed, :no_fv
    has_tpc_veto::Bool
    has_skin_veto::Bool
    n_processed::Int
    deposits::Vector{Deposit}  # FV deposits for :fv events; empty otherwise
end


"""
    process_event(gammas, fk, vol::PhysicalVolume, cfg, rng) -> EventProcessingResult

Process a multi-gamma event through the fast kernel and full stack.
Returns `:fv` with the raw deposit stack for offline SS/MS classification,
`:vetoed` if a veto threshold is exceeded, or `:no_fv` if no gamma
reaches the fiducial volume.
"""
function process_event(gammas, fk::FastKernelGeometry, vol::PhysicalVolume,
                       cfg::SimConfig, rng::AbstractRNG)
    has_tpc_veto = false
    has_skin_veto = false
    n_processed = 0
    fv_deposits = Deposit[]

    empty_deps = Deposit[]

    for gamma in gammas
        result = transport_gamma_fastkernel(gamma, fk, cfg, rng)
        n_processed += 1

        if result.status == :handoff_fv
            deps = propagate_gamma(
                result.energy_MeV,
                vol,
                cfg;
                position=(result.position[1], result.position[2], result.position[3]),
                direction=(result.direction[1], result.direction[2], result.direction[3]),
                rng=rng,
            )
            append!(fv_deposits, deps)
        elseif result.status == :vetoed_tpc
            has_tpc_veto = true
            return EventProcessingResult(:vetoed, has_tpc_veto, has_skin_veto, n_processed, empty_deps)
        elseif result.status == :vetoed_skin
            has_skin_veto = true
            return EventProcessingResult(:vetoed, has_tpc_veto, has_skin_veto, n_processed, empty_deps)
        end
    end

    if !isempty(fv_deposits)
        return EventProcessingResult(:fv, has_tpc_veto, has_skin_veto, n_processed, fv_deposits)
    end
    EventProcessingResult(:no_fv, has_tpc_veto, has_skin_veto, n_processed, empty_deps)
end
