function get_marginalpdf(pdf; dims::Integer = 1)
    if dims isa Integer
        get_marginalpdf(pdf, (dims,))
    elseif dims isa NTuple{N,<:Integer} where N
        get_marginalpdf(pdf, dims)
    end
end
get_marginalpdf(pdf{T,N}, dims::NTuple{N,<:Integer}) where N = pdf
function get_marginalpdf(pdf{T,Np}, _dims::NTuple{Nd,<:Integer}) where {Np,Nd}
    sl = size(pdf.p)
    steps = Np - Nd
    
    dims = [i for i in 1:Np if !(i in _dims)]
    
    


end

function get_marginalpdf(pdf; dims = (1,))