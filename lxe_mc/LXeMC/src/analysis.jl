"""
Offline analysis utilities for FV deposit stacks.

Functions for SS/MS classification and deposit clustering, used
after event processing (not during tracking).
"""


"""
    cluster_deposits_in_z(deposits, dz_cm; E_min=0.0)

Group deposits into clusters whose z-extent is within `dz_cm`.
Clusters below `E_min` are discarded.
Returns `(z_centroid, total_energy)` tuples.
"""
function cluster_deposits_in_z(deposits::Vector{Deposit},
                               dz_cm::Float64;
                               E_min::Float64=0.0)::Vector{Tuple{Float64,Float64}}
    isempty(deposits) && return Tuple{Float64,Float64}[]

    sorted = sort(deposits, by=d -> d.position[3])

    all_clusters = Tuple{Float64,Float64}[]
    cur_zs = Float64[sorted[1].position[3]]
    cur_es = Float64[sorted[1].energy]

    @inbounds for i in 2:length(sorted)
        z = sorted[i].position[3]
        e = sorted[i].energy
        if z - cur_zs[end] < dz_cm
            push!(cur_zs, z)
            push!(cur_es, e)
        else
            E_tot = sum(cur_es)
            z_cent = sum(z * e for (z, e) in zip(cur_zs, cur_es)) / E_tot
            push!(all_clusters, (z_cent, E_tot))
            cur_zs = Float64[z]
            cur_es = Float64[e]
        end
    end

    E_tot = sum(cur_es)
    z_cent = sum(z * e for (z, e) in zip(cur_zs, cur_es)) / E_tot
    push!(all_clusters, (z_cent, E_tot))

    E_min > 0.0 ? filter(c -> c[2] >= E_min, all_clusters) : all_clusters
end


"""
    is_single_site(deposits, dz_cm; E_min=0.0) -> Bool

True if deposits form at most one cluster above `E_min`.
"""
function is_single_site(deposits::Vector{Deposit}, dz_cm::Float64;
                        E_min::Float64=0.0)::Bool
    length(cluster_deposits_in_z(deposits, dz_cm; E_min=E_min)) <= 1
end
