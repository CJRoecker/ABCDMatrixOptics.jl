export propagate, trace

using DataFrames
"""
    propagate(e::Union{Element, Vector{<:Element}}, b)

Propagate a beam `b` either by a single element `e::Element` or `Vector{<:Element}`.

Return is the final beam.
Also available as `e * b`.

## Example
```jldoctest
julia> beam = FreeSpace(10) * GeometricBeam(w=10.0, k=1.0)
GeometricBeam{Float64}(20.0, 1.0, 10.0)

julia> beam = [ThinLens(10), FreeSpace(10)] * GeometricBeam(w=10.0, k=0.0)
GeometricBeam{Float64}(0.0, -1.0, 10.0)

julia> beam.w ≈ 0
true
```
"""
function propagate(m::Matrix, b::Vector{T}) where T
    v = m * b
    return v
end

function propagate(e::Element, b::Vector{T}) where T
    return transfer_matrix(e) * b
end

function propagate(me::MAElement, b::Vector{T}) where T
    return transfer_matrix(me.e) * b + me.v
end

function propagate(m::Matrix, b::GeometricBeam{T}) where T
    e = userDefinedElement(m)
    w, k = transfer_matrix(e) * [b.w, b.k]
    return GeometricBeam{T}(w, k, b.zpos + dz(e))
end

function propagate(e::Element, b::GeometricBeam{T}) where T
    w, k = transfer_matrix(e) * [b.w, b.k]
    return GeometricBeam{T}(w, k, b.zpos + dz(e))
end

function propagate(me::MAElement, b::GeometricBeam{T}) where T
    w, k = transfer_matrix(me.e) * [b.w, b.k] + me.v
    return GeometricBeam{T}(w, k, b.zpos + dz(me.e))
end

function propagate(e::Element, b::GaussianBeam{T}) where T
    ABCD = transfer_matrix(e)
    A, B, C, D = ABCD[1,1], ABCD[1,2], ABCD[2,1], ABCD[2,2]
    q_new = (A * b.q + B) / (C * b.q + D)
    return _n(T, e, b, q_new)
end

_n(::Type{T}, e::Element, b, q_new) where T = GaussianBeam{T}(q_new, b.zpos + dz(e), b.n, b.λ) 
_n(::Type{T}, e::Interface, b, q_new) where T = GaussianBeam{T}(q_new, b.zpos + dz(e), e.n2, b.λ) 



function propagate(es::Vector{<:Element}, b::AbstractBeam)
    return reduce((a,b) -> propagate(b, a), es, init=b)
end

function propagate(es::Vector{<:Any}, b::AbstractBeam)
    return reduce((a,b) -> propagate(b, a), es, init=b)
end

function propagate(es::Vector{<:Any}, b::Vector)
    return reduce((a,b) -> propagate(b, a), es, init=b)
end

# Vector{Union{<:Element, MAElement}}
# Base.:*(e::Union{Matrix, Element, Vector{<:Element}}, b::AbstractBeam) = propagate(e, b)
Base.:*(e::Union{Matrix, Element, MAElement, Vector{<:Element}}, b::AbstractBeam) = propagate(e, b)
Base.:*(e::Union{Matrix, Element, MAElement, Vector{<:Any}}, b::AbstractBeam) = propagate(e, b)
Base.:*(a::Element, b::Element) = transfer_matrix(a) * transfer_matrix(b)
Base.:*(a::Matrix, b::Element) = a * transfer_matrix(b)
Base.:*(a::Vector{<:Element}, b::Vector) = transfer_matrix(a) * b 
Base.:*(a::Element, b::Vector) = transfer_matrix(a) * b 
Base.:*(a::MAElement, b::Vector) = propagate(a, b)
Base.:*(a::Vector{<:Any}, b::Vector{T} where T) = propagate(a,b) 
Base.:*(a::DataFrame, b::AbstractBeam) = propagate(removeMA(a.Element),b) 





# overwritten standard functions for misaligned elements
Base.:+(a::GeometricBeam, b::GeometricBeam) = GeometricBeam(w = (a.w + b.w), k = (a.k + b.k))
Base.:+(a::GeometricBeam, b::MAVector) = GeometricBeam(w = (a.w + b.s), k = (a.k + b.σ))
Base.:+(a::MAVector, b::GeometricBeam) = GeometricBeam(w = (b.w + a.s), k = (b.k + a.σ))
Base.:+(a::Vector, b::MAVector) = [a[1]+b.s, a[2]+b.σ]
Base.:+(a::MAVector, b::MAVector) = MAVector(s = (a.s + b.s), σ = (a.σ + b.σ))
Base.:*(a::Matrix, b::MAVector) = MAVector(s = (a[1,1]*b.s + a[1,2]*b.σ), σ = (a[2,1]*b.s + a[2,2]*b.σ))
Base.:*(a::Element, b::MAVector) = transfer_matrix(a) * b
Base.:*(a::Element, mb::MAElement) = MAElement(transfer_matrix(a) * transfer_matrix(mb.e), transfer_matrix(a) * mb.v) #Element * Misaligned Element, Ergebnis noch verifizieren
Base.:*(a::Matrix, mb::MAElement) = MAElement(a * transfer_matrix(mb.e), a * mb.v) #Ergebnis noch verifizieren
Base.:*(ma::MAElement, b::Element) = MAElement(transfer_matrix(ma.e) * transfer_matrix(b), ma.v) #Misaligned Element * Element, Ergebnis noch verifizieren, 
Base.:*(ma::MAElement, mb::MAElement) = MAElement(transfer_matrix(ma.e) * transfer_matrix(mb.e), ma.e*mb.v + ma.v) #Misaligned Element * Misaligned Element, Ergebnis noch verifizieren, 
# Base.:*(ma::MAElement, b::GeometricBeam) = propagate(ma, b)


"""
    trace(elems::Vector{<:Element}, b0::AbstractBeam)

Trace a beam `b0` through a vector of elements `elems`.
All intermediate states of the beam will be recorded.

Return is a `Vector` of states where the last entry is the final beam.
Final beam is equivalent to `propagate(elems, b0)`.


## Example
```jldoctest
julia> trace([ThinLens(10), FreeSpace(10)], GeometricBeam(w=10.0, k=0.0))
3-element Vector{GeometricBeam{Float64}}:
 GeometricBeam{Float64}(10.0, 0.0, 0.0)
 GeometricBeam{Float64}(10.0, -1.0, 0.0)
 GeometricBeam{Float64}(0.0, -1.0, 10.0)
```
"""
# function trace(elems::Vector{<:Element}, b0::B) where B
#     bs = Vector{B}(undef, 1)
#     bs[1] = b0
#     for (idx, e) in enumerate(elems)
#         push!(bs, trace(e, bs[end])...)
#     end
#     return bs
# end

function trace(elems::Vector{<:Any}, b0::B) where B
    bs = Vector{B}(undef, 1)
    bs[1] = b0
    for (idx, e) in enumerate(elems)
        push!(bs, trace(e, bs[end])...)
    end
    return bs
end

function trace(elems::Element, b0::B) where B
    return [propagate(elems, b0)]
end

function trace(elems::MAElement, b0::B) where B
    return [propagate(elems, b0)]
end
