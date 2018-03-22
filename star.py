import random
import numpy as np
from utils import splitPoly
import matplotlib.patches as patches
import matplotlib.path as path
from matplotlib.transforms import Bbox
import cartopy.crs as ccrs
from spot import Spot

class Star:
    # Stellar Radius in RSun, inclincation in degrees
    # Limb darkening grid resolution (pixel*pixel grid)
    # Rotation period in days
    def __init__(self, params):
        self.radius = params.rad_star
        self.inc = params.sinc
        self.res = params.res
        self.period = params.prot
        self.u = params.u
        self.spots = None
        self.initial_band = params.high_band
        self.low_band = params.low_band
        self.cycle = params.stellar_cycle
        self.active_region = list(self.initial_band)
        self.active_region_vel = [-(params.high_band[0]-params.low_band[0])/self.cycle, -(params.high_band[1] - params.low_band[1])/self.cycle]
        self.params = params # Needed for new spot generation

        # Create globe structure and set up initial projections
        self.globe = ccrs.Globe(semimajor_axis=self.radius, semiminor_axis=self.radius, ellipse='sphere', flattening=1e-9)
        self.rotated_proj = ccrs.RotatedPole(pole_longitude=180, pole_latitude=90-self.inc, central_rotated_longitude=0, globe=self.globe)
        self.geodetic_proj = ccrs.Geodetic(globe=self.globe)
        self.orth_proj = ccrs.Orthographic(globe=self.globe, central_latitude = self.inc, central_longitude=0)

        # Visible surface
        edge = 90
        self.lon1, self.lat1, self.lon2, self.lat2 = -edge, -edge, edge, edge
        
        # Circular grid for limb darkening formula scaled to unity
        x = np.linspace(-1,1,self.res)
        x, y = np.meshgrid(x,x)
        self.grid = np.sqrt(x**2 + y**2)
        self.greater_mask = np.ma.masked_greater(self.grid,1).mask
        self.grid[self.greater_mask] = np.nan
        self.totalGridSquares = self.res**2 - self.greater_mask.sum()
        self.grid_x, self.grid_y = (x*self.radius, y*self.radius) # Re-scale grid back to given star radius
        
        # Unspotted Flux
        self.unspottedFlux = self.limbDarken()
        self.totalUnspottedFlux = self.totalFlux(self.unspottedFlux)
        
        # Spotted Flux
        self.spottedFlux = None
        self.totalSpottedFlux = None

    # Apply quadratic limb darkening to model
    def limbDarken(self):
        mu = np.sqrt(1-self.grid**2)
        mu_1 = 1-mu
        u1 = self.u[0]
        u2 = self.u[1]
        unspottedFlux = 1-u1*mu_1-u2*(mu_1**2)
        
        return unspottedFlux
    
    # Add spots
    def addSpots(self, spots):
        self.spots = spots
        self.spottedFlux = self.mapSpots()
        self.totalSpottedFlux = self.totalFlux(self.spottedFlux)

    # Life Cycle management
    def update(self, cur_phase, t):
        # Update projections
        cur_long = 360*((cur_phase)%1)
        self.updateProjections(cur_long)
        
        # If spots, update them
        if not self.spots == None:
            self.updateSpots(t)
            self.spottedFlux = self.mapSpots()
            self.totalSpottedFlux = self.totalFlux(self.spottedFlux)

    def updateProjections(self, cur_long):
        # Calculte Projections based on current rotation
        self.rotated_proj = ccrs.RotatedPole(pole_longitude=cur_long-180, pole_latitude=90-self.inc, central_rotated_longitude=0, globe=self.globe)
        self.orth_proj = ccrs.Orthographic(globe=self.globe, central_latitude = self.inc, central_longitude=cur_long)
    
    def updateSpots(self, t, dt=0):
        # If no spots then ignore
        if not self.spots == None:
            # Update active latitudes first
            if dt > 0: self.updateActiveRegion(dt)
            
            # Update spots and remove if dead
            doCull = []
            for spot in self.spots:
                if dt > 0: spot.update(self, t, dt)
                if spot.dead: doCull.append(spot)
            
            # Remove dead spots and replace
            if len(doCull) > 0:
                spotsToAddBack = len(doCull)
                for obj in doCull:
                    self.spots.remove(obj)
                for i in range(spotsToAddBack):
                    self.spots.append(Spot.gen_spot(self.params, self, t))
    
    def updateActiveRegion(self, dt):
        self.active_region[0] += dt*self.active_region_vel[0]
        self.active_region[1] += dt*self.active_region_vel[1]
        # Reset when lower than lower band limit
        if self.active_region[0] < self.low_band[0] or self.active_region[1] < self.low_band[1]:
            self.active_region = list(self.initial_band)
    
    # Spot masking and mapping
    def maskPixels(self, path):
        XY = np.dstack((self.grid_x, self.grid_y))
        XY_flat = XY.reshape((-1, 2))
        mask_flat = path.contains_points(XY_flat)
        mask = mask_flat.reshape(self.grid_x.shape)
        return mask
    
    def mapSpots(self):
        # Create new flux array
        spottedFlux = self.unspottedFlux*np.ones(self.unspottedFlux.shape)

        # Map Spots
        for i, spot in enumerate(self.spots):
            # Get polygon
            spotPoly = spot.poly
            
            # Transform spot coords from Geodetic coord system to rotated projection
            spot_vs = self.rotated_proj.transform_points(self.geodetic_proj, spotPoly.vertices[:,0], spotPoly.vertices[:,1])[:,0:2]

            # Split poly to avoid issues at boundary
            polys = splitPoly(spot_vs, 180)
            
            for poly in polys:
                # Get vertices of spot/tissot polygon
                spot_vs = poly.get_xy()

                # Mask in rotated projection (use mpl.Path.clip_to_bbox function)
                spot_path = patches.Path(spot_vs).clip_to_bbox(Bbox([[self.lon1,self.lat1],[self.lon2,self.lat2]]))

                # If spot in visible area calculate flux change
                if len(spot_path.vertices):
                    # Transform masked path to orth projection as this is coordinate space LD grid is in
                    spot_vs = self.orth_proj.transform_points(self.rotated_proj, spot_path.vertices[:,0], spot_path.vertices[:,1])[:,0:2]
                    spot_path = patches.Path(spot_vs)

                    # Find pixels contained in mask and multiply by spot brightnesss
                    mask = self.maskPixels(spot_path)
                    spottedFlux[mask] = spottedFlux[mask]*spot.brightness
        
        return spottedFlux

    # Manage transit
    def transit(self, planet, time, dt):
        I = []
        D = []
        Time = []
        planetPoly = patches.CirclePolygon((0,0),1,100)
        while (planet.isTransiting(time)):
            # Carry on now integrating planet across surface but don't rotate star
            planetFlux = self.unspottedFlux*np.ones(self.unspottedFlux.shape) if self.spottedFlux is None else self.spottedFlux*np.ones(self.spottedFlux.shape)

            # Find position of planet and scale to star's radius
            X, Y = planet.skyPosAtTime(time)
            planet_vx = self.radius*(planetPoly.get_path().vertices[:,0]*planet.rad + X)
            planet_vy = self.radius*(planetPoly.get_path().vertices[:,1]*planet.rad + Y)
            planet_path = path.Path(np.column_stack((planet_vx,planet_vy)))
            
            # Find pixles contained within planet's disk and set to 0
            mask = self.maskPixels(planet_path)
            planetFlux[mask] = 0

            totalTransitFlux = self.totalFlux(planetFlux)
            I.append(totalTransitFlux)
            if self.spots is None:
                D.append(self.totalUnspottedFlux - totalTransitFlux)
            else:
                D.append(self.totalSpottedFlux - totalTransitFlux)
            
            Time.append(time)
            time += dt
            
        return I, D, Time, time
    
    # Helper func to sum over grid of flux values
    def totalFlux(self, flx):
        totalFlux = flx[~self.greater_mask].sum()/self.totalGridSquares
        return totalFlux