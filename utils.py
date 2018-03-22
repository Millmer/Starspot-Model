# Modified tissot function to take custom globe and draw multiple spots at once
# Also, adds spots to given axis and returns matplotlib paths
def tissot(rads_km=None, lons=None, lats=None, n_samples=80, globe=None, ax=None, draw=False, **kwargs):
    import numpy as np
    import cartopy.feature as feature
    import shapely.geometry as sgeom
    from cartopy import geodesic
    import cartopy.mpl.patch as cpatch

    geod = geodesic.Geodesic(radius=globe.semimajor_axis, flattening=globe.flattening)
    geoms = []

    rads_km = np.asarray(rads_km)
    lons = np.asarray(lons)
    lats = np.asarray(lats)

    for i in range(len(lons)):
        circle = geod.circle(lons[i], lats[i], rads_km[i],
                             n_samples=n_samples)
        geoms.append(sgeom.Polygon(circle))

    polys = cpatch.geos_to_path(geoms)

    if draw:
        if ax==None: ax = plt.gca()
        f = feature.ShapelyFeature(geoms, ccrs.Geodetic(globe=globe),
                                       **kwargs)
        ax.add_feature(f)

    return polys

def splitPoly(points, edge, **kwargs):
    import numpy as np
    import matplotlib.patches as patches
    
    if (np.all(points[:,0] >= 0) and np.all(points[:,0] <= edge)) or (np.all(points[:,0] <= 0) and np.all(points[:,0] >= -edge)):
        # All points between 0 and 180 or between -180 and 0
        poly = patches.Polygon(points,closed=True, **kwargs)
        return [poly]
    else:
        mask = np.ma.masked_inside(points[:,0],0,edge).mask
        mask2 = np.ma.masked_inside(points[:,0],-edge,0).mask
        vs1 = points[mask]
        vs2 = points[mask2]

        startx_pts1 = vs1[0,0]
        startx_pts2 = vs2[0,0]

        if startx_pts1 - startx_pts2 < edge:
            poly = patches.Polygon(points,closed=True, **kwargs)
            return [poly]
        else:
            poly1 = patches.Polygon(vs1,closed=True, **kwargs)
            poly2 = patches.Polygon(vs2,closed=True, **kwargs)
            return [poly1, poly2]

# Check if a value is inside an interval
def checkInBetween(val, interval):
    if interval[0] > interval[1]:
        return interval[1] <= val <= interval[0]
    else:
        return interval[0] <= val <= interval[1]
