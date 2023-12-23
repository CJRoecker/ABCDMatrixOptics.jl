function removeMA(vma::Vector{<:Any})
    v = []
    for element in vma
        push!(v, sortElement(element))    
    end
    return convert(Vector{Element},v)
end

function sortElement(el::MAElement)
    return el.e 
end

function sortElement(el::Element)
    return el
end