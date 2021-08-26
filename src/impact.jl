@inline get_r_type(sde::SDE_VI_Oscillator1D) = get_f_type(sde.wall[1].r)

function Wall(_r::Union{<:Number,<:Function},d)
    r = Scalar_Or_Function(_r)
    Wall{typeof(r),typeof(d)}(r,d)
end
function Wall(_r::Union{<:Number,<:Function},d::Integer)
    Wall(_r,float(d))
end

function _walls(_d1,_d2, _r1, _r2)
    d1, d2 = minmax(_d1, _d2)
    if d1 == _d1
        r1, r2 = _r1, _r2
    else
        r1, r2 = _r2, _r1
    end
    Wall(r1,d1), Wall(r2,d2)
end
function symmetricwalls(d,r)
    @assert d>0 "d<0!"
    _walls(-d, d, r, r)
end

## Integration kernel functions
function get_integ_limits(IK::IntegrationKernel{sdeT}) where sdeT<:SDE_VI_Oscillator1D{wT,oscT} where {oscT,wT<:Union{Tuple{wT1},Tuple{wT1,wT2}}} where {wT1<:Wall, wT2<:Wall}
    # assuming a single impact to the closer wall and r < 1
    Δt = get_Δt(IK)
    x = IK.pdgrid.xs[1][IK.idx₁[1]]
    d = IK.sde.wall[IK.impactinterval.wallID[1]].pos 
    r = IK.sde.wall[IK.impactinterval.wallID[1]].r
    vmin, vmax = IK.xs[1], IK.xs[end]
    return get_integ_limits!(IK,vmin, vmax, x, Δt, d, r)
end

# Can change the IK 
@inline function get_integ_limits!(IK,vmin, vmax, x, Δt, d, r)
    vᵢ = (x-d)/Δt # v_{0,I} - 0 and 1 solution for ξ
    if !(vmin < vᵢ < vmax)
        zero_out_impactinterval!(IK.impactinterval)
        return vmin, vmax
    end

    if isapprox(vᵢ,zero(vᵢ),atol = eps(typeof(vᵢ)))
        res = IK.impactinterval.wallID[1] == 1 ? (vmin, zero(vᵢ)) : (zero(vᵢ), vmax)
        zero_out_impactinterval!(IK.impactinterval, Q_atwall = true, Q_zero_wallID = false)
        return res
    end
    
    # TODO in case of variable r use a solver!
    vᵢᵢ = get_vII(x, d, Δt, r)  # v_{0,II} - 1 and 2 solution for ξ
    if vᵢᵢ < vmin
        zero_out_impactinterval!(IK.impactinterval)
        return (vmin, vᵢ)
    elseif vmax < vᵢᵢ
        zero_out_impactinterval!(IK.impactinterval)
        return (vᵢ, vmax)
    end

    if vᵢ<vᵢᵢ
        set_impactinterval!(IK.impactinterval, vᵢᵢ, vmax)
        return (vᵢ, vᵢᵢ, vmax)
    else
        set_impactinterval!(IK.impactinterval, vmin, vᵢᵢ)
        return (vmin, vᵢᵢ, vᵢ)
    end
end

function set_impactinterval!(ii::ImpactInterval, _min, _max)
    ii.Q_atwall[1] = false
    ii.lims[1] = _min
    ii.lims[2] = _max
end

@inline function zero_out_impactinterval!(ii::ImpactInterval; Q_atwall::Bool = false, Q_zero_wallID = true)
    ii.lims .= zero(eltype(ii.lims))
    if Q_zero_wallID
        ii.wallID[1] = zero(eltype(ii.wallID))
    end
    ii.Q_atwall[1] = Q_atwall
    ii
end

@inline get_r(IK::IntegrationKernel) = nothing
function get_r(IK::IntegrationKernel{sdeT}) where sdeT<:SDE_VI_Oscillator1D
    wallID = IK.impactinterval.wallID[1]
    if wallID != 0
        IK.sde.wall[wallID].r
    else
        IK.sde.wall[1].r
    end
end

@inline function get_vII(x,d,Δt, r::Scalar_Or_Function{rT}) where rT<:Number
    (d-x)/(r.f*Δt)
end

function is_impact(v₀,IK::IntegrationKernel)
    false
end
function is_impact(v₀,IK::IntegrationKernel{sdeT}) where sdeT<:SDE_VI_Oscillator1D
    !(IK.impactinterval.Q_atwall[1]) && IK.impactinterval.wallID[1] != 0 && is_impact(v₀,IK.impactinterval)
end
is_impact(v₀,ii::ImpactInterval) = ii.lims[1] < v₀ < ii.lims[2] 

@inline function update_idx1!(IK::IntegrationKernel{sdeT,iT0,iT1},idx₁::NTuple{M,eT}) where {sdeT<:SDE_VI_Oscillator1D,iT0,iT1<:AbstractVector{eT}} where eT<:Number where M
    for m in 1:M
        IK.idx₁[m] = idx₁[m]
    end
    x = IK.pdgrid.xs[1][IK.idx₁[1]]
    set_wallID!(IK,x)
    nothing
end
@inline function set_wallID!(IK::IntegrationKernel{sdeT},x) where sdeT<:SDE_VI_Oscillator1D{wT} where {wT<:Union{Wall, Tuple{wT0}} where wT0<:Wall}
    IK.impactinterval.wallID[1] = 1
end
@inline function set_wallID!(IK::IntegrationKernel{sdeT},x) where sdeT<:SDE_VI_Oscillator1D{wT} where {wT<:Union{Vector{wT0}, Tuple{wT1,wT2}}} where {wT0<:Wall, wT1<:Wall, wT2<:Wall}
    if wT isa Vector && length(IK.sde.walls)==1
        IK.impactinterval.wallID[1] = 1
    else
        IK.impactinterval.wallID[1] = (abs(IK.sde.wall[1].pos - x) < abs(IK.sde.wall[2].pos - x)) ? 1 : 2
    end
end

## VI problem initialize
function create_symmetric_VI_PDGrid(sde::SDE_Oscillator1D, d, r, v_ax::aT, Nₓ::Integer; x_interpolation = :chebyshev, kwargs...) where aT<:Axis
    walls = symmetricwalls(d,r);
    x_ax = Axis(-d, d, Nₓ, interpolation = x_interpolation)
    vi_sde = SDE_VI_Oscillator1D(sde,walls)
    vi_sde, PDGrid(vi_sde, x_ax, v_ax; kwargs...)
end