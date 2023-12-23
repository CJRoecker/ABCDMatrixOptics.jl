
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

# Genereller Workflow: mit Rays Raw-Beams erzeugen. Einen Raw-Beam durch system tracen => in Geometric_quiver schmeißen, fertig!

# funktioniert nicht da plots.jl nicht teil des package ist => Recipe
# ersatzworkflow: Funktion kopieren und in REPL einfügen
function Geometric_quiver(v::Vector{GeometricBeam{T}} where T; createnew=true)
    i=2
    if createnew
        p1 = plot()
    end
    for i in eachindex(v)
        if i==1 
            continue 
        elseif i==2
            if createnew
                p1 = quiver([v[i-1].zpos], [v[i-1].w], quiver=([v[i].zpos-v[i-1].zpos], [v[i].w-v[i-1].w]), show=true, c=:black)
            else    
                p1 = quiver!([v[i-1].zpos], [v[i-1].w], quiver=([v[i].zpos-v[i-1].zpos], [v[i].w-v[i-1].w]), show=true, c=:black)
            end
        else
            quiver!([v[i-1].zpos], [v[i-1].w], quiver=([v[i].zpos-v[i-1].zpos], [v[i].w-v[i-1].w]), c=:black)
        end
    end
    return p1
end

function create_rays(;n = 5, xmin = -10, xmax = 10, kmin = -0.1, kmax = 0.1)
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

function trace_all(rays::Vector{GeometricBeam{T}} where T, M)
    collection = Vector{Vector{GeometricBeam{Float64}}}(undef,length(rays))
    for i in eachindex(collection)
        collection[i] = trace(M, rays[i])
    end
    return collection
end

# same wie für oben: muss in REPL kopiert werden 
function quiver_all(collection)
    for i in eachindex(collection)
        if i == 1
            display(Geometric_quiver(collection[i]))
        else
            display(Geometric_quiver(collection[i], createnew = false))
        end    
    end
end

#exampe usage
# p,d = ABCDMatrixOptics.create_rays(n=10)
# cld = ABCDMatrixOptics.trace_all(d,M)  # collection of traces of all divergent rays
# clp = ABCDMatrixOptics.trace_all(p,M)  # collection of traces of all parallel rays
# quiver_all(cld)
# quiver_all(clp)