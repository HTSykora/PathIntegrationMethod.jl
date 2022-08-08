(mpdf::MarginalPDF)(x...) = mpdf.pdf(x...)

function initialise_mPDF(pdf,IDs)
    Tuple(initialise_mPDF(pdf,_IDs) for _IDs in IDs)
end
function initialise_mPDF(pdf::InterpolatedFunction{T,N},IDs::Union{Integer,NTuple{n,<:Integer}}) where {T, N, n}
    dims = get_mpdf_intdims(pdf,IDs)
    
    mpdf = get_mpdf0(pdf,IDs)
    p0 = reshape(mpdf.p, (i in dims ? 1 : size(pdf.p,i) for i in 1:N)...)
    temp = similar(pdf.p);
    wMX = get_wMX(pdf,dims);

    MarginalPDF(mpdf,IDs,wMX, temp, p0,dims)
end

function get_mpdf_intdims(pdf::InterpolatedFunction{T,N},ID::Integer) where {T, N}
    @assert ID ≤ N "ID > $N"
    @assert 0 < ID "ID < 0"
    
    Tuple(i for i in 1:N if i != ID)
end
function get_mpdf_intdims(pdf::InterpolatedFunction{T,N},IDs::NTuple{n,Integer}) where {n,T,N}
    @assert maximum(IDs) ≤ N "maximum(IDs) > $N"
    @assert 0 < minimum(IDs) "minimum(IDs) < 0"
    
    Tuple(i for i in 1:N if !(i ∈ IDs))
end

function get_mpdf0(pdf,ID::Integer)
    InterpolatedFunction(eltype(pdf.p), pdf.axes[ID])
end

function get_mpdf0(pdf,IDs::NTuple{n,Integer}) where n
    InterpolatedFunction(eltype(pdf.p), getindex.(Ref(pdf.axes),IDs)...)
end

function get_wMX(pdf::InterpolatedFunction{T,N},dims) where {T,N}
    wMX = similar(pdf.p)
    fill!(wMX,one(T))

    for dim in dims
        multiply_weightMX!(wMX,pdf.axes[dim].wts,dim)
    end
    wMX
end

function multiply_weightMX!(res::AbstractArray{T,N}, ws::AbstractVector{wT},dim::Integer) where {T,T1,wT,N,N1}
    sl = Slicer(dim,N)

    @inbounds for (i,w) in enumerate(ws)
        update_sliceridx!(sl,i)
        vres = view(res,sl)
        vres .= vres .* w
    end

    res
end


function update_mPDF!(mpdf::MarginalPDF,pdf::InterpolatedFunction; detached = false, kwargs...)
    mpdf.temp .= mpdf.wMX .* pdf.p

    mpdf.p0 .= sum(mpdf.temp, dims=mpdf.dims)
    if detached # In case the structure is reconstructed using JLD2 we should use `detached = true`
        for i in eachindex(mpdf.p0)
            @inbounds mpdf.pdf.p[i] = mpdf.p0[i]
        end
    end
    nothing
end

# Utils
function Slicer(n,N, idx0 = 1; size_t = Int)
    slv = [size_t(idx0)];
    slicer = ((Colon() for _ in 1:n-1)...,slv,(Colon() for _ in n+1:N)...)
    Slicer{n,N,size_t,typeof(slicer)}(slicer)
end
function update_sliceridx!(sl::Slicer{n,N,idT}, idx::idT) where {n,N,idT}
    sl.slicer[n][1] = idx
    nothing
end
Base.view(p::AbstractArray,sl::Slicer) = Base.view(p,sl.slicer...)

