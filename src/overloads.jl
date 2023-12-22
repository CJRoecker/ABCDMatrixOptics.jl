import Base.+

+(x::GeometricBeam, y::GeometricBeam) =
begin
    w = x.w + y.w
    k = x.k + y.k
    res = GeometricBeam(w = w, k = k)
    return res  
end

