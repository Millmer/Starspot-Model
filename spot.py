import random
import copy
import numpy as np
from utils import tissot
import settings

class Spot:
    # lon, lat in terms of stellar lon lat
    # Radius in terms of stellar radius.
    # Brightness in relation to max stellar intensity = 1
    # Latitude velocity in deg/day
    # Lifetime in days
    def __init__(self, lon, lat, r, b, lv, life, globe, time):
        self.lon = lon
        self.lat = lat
        self.rad = r*globe.semimajor_axis
        self.brightness = b
        self.lat_vel = lv
        self.life = life
        self.dead = False
        self.globe = globe
        self.poly = tissot([r],[lon],[lat],100,globe)[0]
        settings.GLOBAL_ALL_SPOTS.append([copy.deepcopy(self), time, time + self.life])
        self.global_index = len(settings.GLOBAL_ALL_SPOTS)-1
    
    def update(self, star, time, dt):
        # Update life
        self.life -= dt
        
        # Check if dead
        if self.life <= 0 or np.fabs(self.lat) < star.active_region[0]:
            settings.GLOBAL_ALL_SPOTS[self.global_index][2] = time # Update death time
            self.dead = True
        else:
            # Update lat position and recalculate polygon
            self.lat += dt*self.lat_vel
            self.poly = tissot([self.rad],[self.lon],[self.lat],100,self.globe)[0]
    
    def __str__(self):
        desc = 'r = {0:.2f}, (lon, lat) = ({1:.2f}, {2:.2f})\nb = {3:.2f}, lv = {4:.2f}, l = {5:.2f}'.format(self.rad, self.lon, self.lat, self.brightness, self.lat_vel, self.life)
        return desc
    
    @classmethod
    def gen_spot(cls, params, star, time):
        b = random.uniform(params.b[0], params.b[1]) # In relation to max stellar intensity
        r = random.uniform(params.spot_rads[0], params.spot_rads[1]) # In terms of stellar radius
        lat = random.uniform(star.active_region[0], star.active_region[1]) # Following sun patterns
        # Decide hemisphere of spot generation
        lat = -1*lat if random.uniform(0, 1) > params.spot_asymmetry else lat
        lon = random.uniform(0, 359) # Randomly anywhere for starting
        lv = star.active_region_vel[0] # Velocity based on given stellar cycle
        # Ensure migration to centre
        if np.sign(lat) == np.sign(lv):
            lv *= -1
        life = random.uniform(params.spot_lives[0], params.spot_lives[1]) # Life
        return Spot(lon, lat, r, b, lv, life, star.globe, time)