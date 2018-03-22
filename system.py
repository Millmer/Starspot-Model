import numpy as np
import cartopy.crs as ccrs
from peakdetect import peakdetect
from utils import checkInBetween
import settings

class System:
    def __init__(self, star, planet, params):
        self.star = star # Star object class
        self.planet = planet # Planet object class
        self.rot_dt = params.rotation_dt # Time step outside transit
        self.tra_dt = params.transit_dt # Time step inside transit
        self.adj_dist = params.adjustDistance # Distance from star at which to adjust timestep (recommend 1.5Rstar)
    
    # Helper func to manage if no spots
    def calculateFlux(self, time):
        star_phase = time/self.star.period
        self.star.update(star_phase, time)
        i_avg = self.star.totalUnspottedFlux if self.star.totalSpottedFlux is None else self.star.totalSpottedFlux
        return i_avg

    # Params:
    #   Duration - Length of sim [days]
    #   transitOnly - Only calc flux for transits [bool], Default True
    #   progUpdate - print sim progres every x%, Default 25
    # Returns:
    #   I - list of intensity values
    #   D - list of differences between outside transit intensities and inside
    #   Time - list of time steps
    #   TTimes - list of tuples of transit times (start, end)
    #   SAR - list of star's active region over time
    def simulate(self, duration, transitOnly=True, progUpdate=25):
        I = []
        D = []
        Time = []
        TTimes = []
        SAR = []
        time = 0
        old_prog = None
        while (time < duration):
            # Print status
            prog = int(np.floor(time/duration*100))
            if prog % progUpdate == 0 and not (prog == old_prog):
                print('Simulation {0}% complete'.format(prog))
                old_prog = prog
            
            SAR.append((list(self.star.active_region), time))
            
            if self.planet.withinDistanceFromStar(time, self.adj_dist):
                I.append(self.calculateFlux(time))
                D.append(np.nan)
                Time.append(time)

                if self.planet.isTransiting(time):
                    t1 = time
                    i,d,t,time = self.star.transit(self.planet, time, self.tra_dt)
                    I.extend(i)
                    D.extend(d)
                    Time.extend(t)
                    t2 = time
                    TTimes.append((t1,t2))

                    I.append(self.calculateFlux(time))
                    D.append(np.nan)
                    Time.append(time)
                    
                time += self.tra_dt
                self.star.updateSpots(time, self.tra_dt)
            else:
                i_avg = np.nan if transitOnly else self.calculateFlux(time)
                I.append(i_avg)
                D.append(np.nan)
                Time.append(time)

                time += self.rot_dt
                self.star.updateSpots(time, self.rot_dt)

        # Normalise
        I = [i/self.star.totalUnspottedFlux for i in I]

        return I, D, Time, TTimes, SAR

    def butterfly(self, residuals, noise, Time, TimeError, TTimes):
        peaks = peakdetect(residuals, Time, 1, noise)[0]
        if len(peaks) == 0: 
            raise RuntimeError("No peaks found")

        foundSpots = []
        for peak in peaks:
            t = peak[0]
            h = peak[1]
            for interval in TTimes:
                if checkInBetween(t, interval):
                    X, Y = self.planet.skyPosAtTime(t)
                    X_min, Y_min = self.planet.skyPosAtTime(t - TimeError)
                    X_max, Y_max = self.planet.skyPosAtTime(t + TimeError)
                    fs = MultipurposeSpot((X,Y), ((X_min, Y_min), (X_max, Y_max)), h, t, self.star)
                    foundSpots.append(fs)
                    continue

        realSpots = []
        for t in Time:
            for obj in settings.GLOBAL_ALL_SPOTS:
                # Grab spot and creation/death times
                realSpot = obj[0]
                creation = obj[1]
                death = obj[2]
                # Check spot didn't die
                if checkInBetween(t, [creation, death]):
                    # Correct for rotation star centre long
                    timeCorrectedLong = 360*((t/self.star.period)%1)
                    # Visible region bounds
                    bound_left = 270
                    bound_right = 90
                    # Correct spot and shift based on moving star long
                    spotCorrectedLong = (realSpot.lon - timeCorrectedLong)%360
                    timeCorrectedLat = realSpot.lat + (t-creation)*realSpot.lat_vel
                    # Check if spot was visible at that time
                    if spotCorrectedLong > bound_left or spotCorrectedLong < bound_right:
                        rs = MultipurposeSpot(None, None, realSpot.brightness, t, None, (None, timeCorrectedLat))
                        realSpots.append(rs)

        matchedSpots = []
        for fs in foundSpots:
            distance = 1000
            idx = np.nan
            for i, rs in enumerate(realSpots):
                val = np.sqrt(np.abs(fs.timeFound-rs.timeFound)**2 + np.abs(fs.coords[1]-rs.coords[1])**2)
                if val < distance:
                    distance = val
                    idx = i
            matchedSpots.append(realSpots[idx])

        return foundSpots, realSpots, matchedSpots

class MultipurposeSpot():
    def __init__(self, pos, pe, b, tf, star, crds=None):
        self.coords = crds
        self.latWidth = None

        # If not given coords, calculate them
        if crds == None:
            self.coords = self.findProjectedPosition(pos, star)
            min_lat = self.findProjectedPosition(pe[0], star)
            max_lat = self.findProjectedPosition(pe[1], star)
            if min_lat is np.nan:
                min_lat = self.coords
            if max_lat is np.nan:
                max_lat = self.coords
            self.latWidth = np.abs(max_lat[1] - min_lat[1])


        self.brightness = b
        self.timeFound = tf
    
    def findProjectedPosition(self, pos, star):
        # Check inside disk
        if np.sqrt(pos[0]**2 + pos[1]**2) <= 1:
            # Scale calculated X,Y to actual star radius
            pos *= star.radius
            # Use CartoPy to transform from Ortho coords
            unitPole = ccrs.RotatedPole(pole_longitude=180, pole_latitude=90,
                                    central_rotated_longitude=0, globe=star.globe)
            return unitPole.transform_points(star.orth_proj, np.asarray([pos[0]]), np.asarray([pos[1]]))[:,0:2][0]
        else:
            return np.nan