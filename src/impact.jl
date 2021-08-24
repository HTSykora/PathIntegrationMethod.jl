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
function get_integ_limits(IK::IntegrationKernel{sdeT}) where sdeT<:SDE_VI_Oscillator1D{wT,oscT} where {oscT,wT<:Tuple{wT1,wT2}} where {wT1<:Wall, wT2<:Wall}
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
    if isapprox(vᵢ,zero(vᵢ),atol = eps(typeof(vᵢ)))
        res = IK.impactinterval.wallID[1] == 1 ? (vmin, vᵢ) : (vᵢ, vmax)
        zero_out_impactinterval!(IK.impactinterval, Q_atwall = true, zero_wallID = false)
        return res
    elseif vmin < vᵢ < vmax
        # TODO in case of variable r use a solver!
        vᵢᵢ = get_vII(x,d,Δt, r)  # v_{0,II} - 1 and 2 solution for ξ
        if vmin < vᵢᵢ < vmax
            if vᵢ<vᵢᵢ
                IK.impactinterval.lims[1] = vᵢᵢ
                IK.impactinterval.lims[2] = vmax
                return (vᵢ, vᵢᵢ, vmax)
            else
                IK.impactinterval.lims[1] = vmin
                IK.impactinterval.lims[2] = vᵢᵢ
                return (vmin, vᵢᵢ, vᵢ)
            end
        else
            zero_out_impactinterval!(IK.impactinterval)
            return vᵢ<vᵢᵢ ? (vᵢ, vmax) : (vmin, vᵢ)
        end
    else
        zero_out_impactinterval!(IK.impactinterval)
        return vmin, vmax
    end
end
@inline function zero_out_impactinterval!(ii::ImpactInterval; Q_atwall::Bool = false, zero_wallID = true)
    ii.lims .= zero(eltype(ii.lims))
    if zero_wallID
        ii.wallID[1] = zero(eltype(ii.wallID))
    end
    if Q_atwall && !ii.Q_atwall[1]
        ii.Q_atwall[1] = true
    elseif ii.Q_atwall[1]
        ii.Q_atwall[1] = false
    end
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


## VI problem initialize
function create_symmetric_VI_PDGrid(sde::SDE_Oscillator1D, d, r, v_ax::aT, Nₓ::Integer; x_interpolation = :chebyshev, kwargs...) where aT<:Axis
    walls = symmetricwalls(d,r);
    x_ax = Axis(-d, d, Nₓ, interpolation = x_interpolation)
    vi_sde = SDE_VI_Oscillator1D(sde,walls)
    vi_sde, PDGrid(vi_sde, x_ax, v_ax; kwargs...)
end