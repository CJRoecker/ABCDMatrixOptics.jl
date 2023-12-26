using RecipesBase

@recipe function f(rc::ABCDMatrixOptics.GeoRayCoords{T} where T)
    # set a default value for an attribute with `-->`
    xlabel --> "Distance z"
    yguide --> "Distance from reference axis"
    markershape --> :cross

    #set the stage for the quiver plot: quiver(x,y,quiver=(u,v))
    seriestype := :quiver
    quiver --> (rc.dz,rc.dw)

    # get the seriescolor passed by the user
    c = get(plotattributes, :seriescolor, :auto)

    # return data
    rc.z, rc.w
end