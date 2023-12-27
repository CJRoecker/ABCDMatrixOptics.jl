export removeMA, analyzeSystemGeometrically, create_rays, trace_all, getGeoRayCoordinates, getStabilityParameter, getqParameter
using RecipesBase

"""
    removeMA(vma::Vector{<:Any})

Converts a misaligned System into an aligned system for comfortable plotting
    of gaussian beam propagation with standard functions. 
Removes misalignment vectors and does type conversion
"""
function removeMA(vma::Vector{<:Any})
    v = []
    for element in vma
        push!(v, sortElement(element))    
    end
    return convert(Vector{Element},v)
end

"""
    sortElement(el::MAElement)

returns the Matrix of the (mis)aligned optical element
"""
function sortElement(el::MAElement)
    return el.e 
end

function sortElement(el::Element)
    return el
end

# General workflow: use create rays to generate beams. 
# trace the generated beam through the System

#exampe usage
# par_coordinates, div_coordinates = analyzeSystemGeometrically(SystemMatrix=M)
# plot(par_coordinates) # automatically generates a quiver plot for the ray bundle
# LEGACY
# quiver(div_coordinates.z,div_coordinates.w, quiver=(div_coordinates.dz,div_coordinates.dw))
# quiver(par_coordinates.z,par_coordinates.w, quiver=(par_coordinates.dz,par_coordinates.dw))

"""
    create_rays(;n = 5, xmin = -10, xmax = 10, kmin = -0.1, kmax = 0.1)

create a number of n parallel(k=0) and diverging(w=0) rays ranging
from xmin to xmax and kmin and kmax, respectively.
"""
function create_rays(;n = 5, xmin = -10e-3, xmax = 10e-3, kmin = -0.1e-3, kmax = 0.1e-3)
    parallel_rays = Vector{GeometricBeam{Float64}}(undef, 0)
    divergent_rays = Vector{GeometricBeam{Float64}}(undef, 0)

    for w in range(xmin, xmax, n)
        push!(parallel_rays, GeometricBeam(w=w, k=0)) 
    end
    for k in range(kmin, kmax, n)
        push!(divergent_rays, GeometricBeam(w=0, k=k)) 
    end

    return parallel_rays, divergent_rays
end

"""
    trace_all(rays::Vector{GeometricBeam{T}} where T, M)

This function traces n rays through a system M with m elements
    Tracing one ray through a system M with m elements yields m new rays
    The function returns a Vector of length n containing each m Vectors  
"""
function trace_all(rays::Vector{GeometricBeam{T}} where T, M)
    collection = Vector{Vector{GeometricBeam{Float64}}}(undef,length(rays))
    for i in eachindex(collection)
        collection[i] = trace(M, rays[i])
    end
    return collection
end

"""
    getGeoRayCoordinates(vect::Vector{GeometricBeam{T}} where T; x =[], y=[], u=[], v=[])

Extract the coordinates of one GeometricRay
"""
function getGeoRayCoordinates(vect::Vector{GeometricBeam{T}} where T; 
    x =Vector{Float64}(undef,0), y=Vector{Float64}(undef,0), 
    u=Vector{Float64}(undef,0), v=Vector{Float64}(undef,0))
    for i in eachindex(vect)
        if i==1 
            continue 
        else
            push!(x, vect[i-1].zpos)
            push!(y, vect[i-1].w)   
            push!(u, vect[i].zpos-vect[i-1].zpos)
            push!(v, vect[i].w-vect[i-1].w)
        end
    end
    return [x,y,u,v]    
end

"""
    getGeoRayCoordinates(vectvect::Vector{Vector{GeometricBeam{T}}} where T, coordinateVect=[])

Extract all the coordinates of a Vector of Geometric Rays
"""
function getGeoRayCoordinates(vectvect::Vector{Vector{GeometricBeam{T}}} where T, coordinateVect=[])
    for (i,vect) in enumerate(vectvect)
        if i ==1
            coordinateVect = getGeoRayCoordinates(vect)
        else
            coordinateVect = getGeoRayCoordinates(vect; x=coordinateVect[1], y=coordinateVect[2], u=coordinateVect[3], v=coordinateVect[4])
        end
    end
    # return coordinateVect
    return GeoRayCoords(coordinateVect[1], coordinateVect[2], coordinateVect[3], coordinateVect[4])
end

"""
    analyzeSystemGeometrically(;SystemMatrix, n=10)

convenience wrapper for standard workflow: 
    1. create rays for calculation
    2. Trace rays through the system
    3. Extract and return the coordinates of the vectors
"""
function analyzeSystemGeometrically(SystemMatrix; n::Integer=10, x::Number=10e-3, k::Number=0.1)
    par_rays, div_rays = ABCDMatrixOptics.create_rays(n=n,xmin=-x,xmax=x,kmin=-k,kmax=k)
    div_collection = ABCDMatrixOptics.trace_all(div_rays, SystemMatrix);  # collection of traces of all divergent rays
    par_collection = ABCDMatrixOptics.trace_all(par_rays, SystemMatrix);  # collection of traces of all parallel rays
    div_coordinates = ABCDMatrixOptics.getGeoRayCoordinates(div_collection);
    par_coordinates = ABCDMatrixOptics.getGeoRayCoordinates(par_collection);    
    return par_coordinates, div_coordinates
end

function analyzeSystemGeometrically(SystemMatrix, ray::GeometricBeam{T} where T)
    rays = [ray]
    traced = ABCDMatrixOptics.trace_all(rays, SystemMatrix);  # collection of traces of all divergent rays
    coordinates = ABCDMatrixOptics.getGeoRayCoordinates(traced);
    return coordinates
end

function analyzeSystemGeometrically(SystemMatrix, ray::Vector{GeometricBeam{T}} where T)
    traced = ABCDMatrixOptics.trace_all(ray, SystemMatrix);  # collection of traces of all divergent rays
    coordinates = ABCDMatrixOptics.getGeoRayCoordinates(traced);
    return coordinates
end


"""
    getStabilityParameter()

returns the stability parameter (A+D)/2 for a given optical System Roundtrip Matrix
"""

function getStabilityParameter(System_RT_Matrix::Matrix)
    return (System_RT_Matrix[1,1]+System_RT_Matrix[2,2])/2
end

function getStabilityParameter(System_RT_Vect::Vector{Any})
    M = transfer_matrix(removeMA(System_RT_Vect))
    return getStabilityParameter(M)
end


"""
    getqParameter()

returns the q-parameter of the resonator
"""
function getqParameter(System_RT_Matrix::Matrix)
    A = System_RT_Matrix[1,1]
    B = System_RT_Matrix[1,2]
    C = System_RT_Matrix[2,1]
    D = System_RT_Matrix[2,2]

    q_plus = ((A-D) + sqrt(Complex((A+D)^2-4)))/(2*C)
    q_minus = ((A-D) - sqrt(Complex((A+D)^2-4)))/(2*C)

    if imag(q_plus) > 0
        return q_plus
    else
        return q_minus
    end
end

function getqParameter(System_RT_Vect::Vector{Any})
    M = transfer_matrix(removeMA(System_RT_Vect))
    return getqParameter(M)
end
