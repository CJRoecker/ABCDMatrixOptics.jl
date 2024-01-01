export Element, FreeSpace, Interface, ThinLens, MAThinLens, ThickLens, Mirror, MAMirror
export transfer_matrix

abstract type Element{T} end

"""
    FreeSpace(dz)

Construct a free space propagation over distance `dz`.
"""
@with_kw_noshow struct FreeSpace{T<:Number} <: Element{T}
    dz::T
end

"""
    ThickLens(;  R1, R2, t, n_lens=1.5 n1=1.0, n2=1.0)

Construct a thick lens with the keywords:

* `R1` radius of curvature of first surface
* `R2` radius of curvature of second surface
* `t`: thickness of lens
* `n_lens=1.5` refractive index of lens
* `n1=1`: refractive index of the medium of the incoming side
* `n2=1`: refractive index of the medium of the exiting side

"""
struct ThickLens{T<:Number} <: Element{T}
    R1::T
    R2::T
    t::T
    n_lens::T
    n1::T
    n2::T
end

function ThickLens(; R1, R2, t, n_lens=1.5, n1=1.0, n2=1.0)
    R1, R2, t, n_lens, n1, n2 = promote(R1, R2, t, n_lens, n1, n2)
    return ThickLens{typeof(R1)}(R1, R2, t, n_lens, n1, n2)
end

struct Interface{T<:Number} <: Element{T}
    n1::T
    n2::T
    R::T
end

"""
    Interface(n1, n2)

Creates a flat interface with refractive index `n1` on the entering side
and `n2` on the new medium.
"""
Interface(n1::Integer, n2::Integer) = Interface{Float64}(promote(n1, n2, Inf)...)
Interface(n1, n2, R=Inf) = Interface{promote_type(typeof(n1), typeof(n2), typeof(R))}(promote(n1, n2, R)...)

"""
    Interface(; n1, n2, R=Inf)

Creates a curved interface with radius `R` and with refractive index `n1` on the entering side
and `n2` on the new medium.
"""
Interface(; n1, n2, R=Inf) = Interface{promote_type(typeof(n1), typeof(n2), typeof(R))}(promote(n1, n2, R)...)



@with_kw_noshow struct ThinLens{T<:Number} <: Element{T}
    f::T
    θ::T = typeof(f)(0)   # misalignment: angle
    Δ::T = typeof(f)(0)   # misalignment: offset 
end

"""
    ThinLens(f)

Creates a thin lens with focal length `f`.
"""
ThinLens(f::T) where {T} = ThinLens{T}(f, 0, 0)


"""
    ThinLens(R1, R2, n_lens, n)

Creates a thin lens defined by the first radius of curvature `R1`, the second `R2`.
The lens refractive index is `n_lens` and the outer refractive index is `n`.
"""
ThinLens(R1, R2, n_lens=1.5, n=1.0) = ThinLens(inv((n_lens - n) / n * (1 / R1 - 1 / R2)))

"""
    Mirror(R=Inf)

Mirror with radius of curvature `R`.
Per default `Inf`, so a flat mirror.
"""
@with_kw_noshow struct Mirror{T} <: Element{T}
    R::T = Inf
    θ::T = 0   # misalignment: angle
    Δ::T = 0   # misalignment: offset 
end

"""
    Mirror(R)

Creates a Mirror with Radius of curvature `R`.
"""
Mirror(R::T) where {T} = Mirror{T}(R, 0, 0)

struct userDefinedElement{T<:Number} <: Element{T}
    m::Matrix{T}
end

"""
    MAVector(T)

Construct a Vector containing information on the misalignment of a standard element.
"""
struct MAVector{T<:Number}
    s::T
    σ::T
end

MAVector(e::ThinLens) = MAVector{typeof(e.Δ/e.f)}(0, e.Δ/e.f)
MAVector(e::Mirror) = MAVector{typeof(e.Δ/e.R)}(0, 2*(e.Δ/e.R+e.θ))
MAVector(;s, σ) = MAVector(s, σ) 


"""
    MAElement(T)

Construct a misaligned Element containing a standard element and a misalignment Vector.
"""
struct MAElement{T<:Number}
    e::Element{T}
    v::MAVector{T}
end

MAElement(e::ThinLens) = MAElement{promote_type(typeof(e.f), typeof(e.Δ/e.f))}(e, MAVector(e))
MAElement(e::Mirror) = MAElement{promote_type(typeof(e.R), typeof(e.θ), typeof(e.Δ), typeof(e.Δ/e.R))}(e, MAVector(e))
MAElement(e::Matrix, v::MAVector) = MAElement(userDefinedElement(e),v)


"""
    MAThinLens(f)

Creates a thin lens with focal length `f`, angular misalignment 'θ', and displacement 'Δ'.
"""
MAThinLens(f, θ, Δ) = MAElement(ThinLens{promote_type(typeof(f), typeof(θ), typeof(Δ), typeof(Δ/f))}(promote(f, θ, Δ)...))

"""
    Mirror(R, θ, Δ)

Creates a misaligned Mirror with Radius of curvature `R`, angular misalignment 'θ', and displacement 'Δ'..
"""
MAMirror(R, θ, Δ) = MAElement(Mirror{promote_type(typeof(R), typeof(θ), typeof(Δ), typeof(Δ/R))}(promote(R, θ, Δ)...))

# definitions of dz
"""
    dz(element::Element)

Returns how much an element changes the optical distance `z`.
"""
dz(e::FreeSpace) = e.dz
dz(e::Interface{T}) where {T} = zero(T)
dz(e::ThinLens{T}) where {T} = zero(T)
dz(e::Mirror{T}) where {T} = zero(T)
dz(e::ThickLens) = e.t
# dz(e::Matrix) = Inf
dz(e::Matrix) = e[1,2]  # total propagation length corresponds to Element B 
dz(e::userDefinedElement) = dz(e.m)



"""
    transfer_matrix(element::Element) 

Returns the Ray Transfer (ABCD) matrix associated with the given, optical element.
"""
transfer_matrix(e::Matrix) = e
transfer_matrix(e::userDefinedElement) = transfer_matrix(e.m)
transfer_matrix(e::Interface) = [1 0; ((e.n1-e.n2)/(e.R*e.n2)) (e.n1/e.n2)]
transfer_matrix(e::ThinLens) = [1 0; -1/e.f 1]
transfer_matrix(e::Mirror) = [1 0; -2/e.R 1]
transfer_matrix(e::FreeSpace) = [1 e.dz; 0 1]
transfer_matrix(e::ThickLens) = transfer_matrix([Interface(n1=e.n1, n2=e.n_lens, R=e.R1),
    FreeSpace(e.t),
    Interface(n1=e.n_lens, n2=e.n2, R=e.R2)])
transfer_matrix(e::MAElement) = transfer_matrix(removeMA(e))




"""
    transfer_matrix(elements)

Returns the Ray Transfer (ABCD) matrix associated with
an optical system described by a collection (e.g. a vector or
iteration) of optical elements.
"""
function transfer_matrix(elements::Vector{<:Element})
    return mapfoldr(transfer_matrix, (a, b) -> b * a, elements)
end

function transfer_matrix(elements::Vector{Union{<:Element, MAElement}})
    return mapfoldr(transfer_matrix, (a, b) -> b * a, elements)
end

function transfer_matrix(elements::Vector{Any})
    return mapfoldr(transfer_matrix, (a, b) -> b * a, elements)
end


"""
    discretize(e)

Discretizes the elements for plots. Nothing is done expect for FreeSpace, which is split up
"""

discretize(e::FreeSpace, N::Int) = fill(FreeSpace(e.dz / N), N)
discretize(e::Element, N::Int) = e
discretize(els::Vector{<:Element}, N::Int) = vcat(discretize.(els, Ref(N))...)

"""
    Base.isapprox(a::Vector{<:Element}, b::Vector{<:Element})

Compare two vectors of elements using Base.isapprox for each element's
ray matrix (ABCD entries). Does consequently not consider one
discretization of element FreeSpace different from another, or one
realization of an imaging system from another as long as both achieve
(within tolerances) the same imaging.

!!! note
    The `atol` (absolute tolerance) parameter can be used but is
    typically nonsensical as it will be used for each of the
    ray matrix entries ABCD which usually differ vastly in magnitude.

"""
Base.isapprox(
    a::Union{Element,Vector{<:Element}}, b::Union{Element,Vector{<:Element}}; kwargs...
) = isapprox(transfer_matrix(a), transfer_matrix(b); kwargs...)

Base.isapprox(a::Matrix, b::Union{Element,Vector{<:Element}}; kwargs...) = isapprox(a, transfer_matrix(b); kwargs...)
Base.isapprox(a::Union{Element,Vector{<:Element}}, b::Matrix; kwargs...) = Base.isapprox(b, a)
