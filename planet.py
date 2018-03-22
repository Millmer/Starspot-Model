import numpy as np

class Planet:
    # Beta is spin-orbit alignment in degrees (lambda is reserved in python)
    # Impact parameter in star_radii (shift in latitude on stellar disk)
    # Orbit period in days
    def __init__(self, params):
        # Planet transit chord
        self.rad = params.rp # Planet radius [Rstar]
        self.beta = params.beta # Spin-orbit alignement [deg]
        self.per = params.per # Orbital period [days]
        self.w = params.w # Argument of periapsis [deg]
        self.a = params.a # Semi-major axis [Rstar]
        self.inc = params.oinc # Orbirtal inclincation [deg]
        self.t0 = params.t0 # Time of inferior conjunction [days]

    def t_periastron(self):
        f = np.pi/2. - np.deg2rad(self.w) # True Anomaly at mid transit
        phase = f/2/np.pi
        return self.t0 - self.per*phase

    def inFrontOfStar(self, f, omega, inc):
        return np.sin(f + omega) * np.sin(inc) > 0 # Z > 0

    def inFrontOfStarAtTime(self, time):
        tp = self.t_periastron()
        omega = np.deg2rad(self.w)
        i = np.deg2rad(self.inc)
        f = ((time - tp)/self.per - (int)((time-tp)/self.per))*2*np.pi
        return self.inFrontOfStar(f, omega, i)

    def r_skyAtTime(self, time):
        # Update orbital params
        tp = self.t_periastron()
        omega = np.deg2rad(self.w)
        i = np.deg2rad(self.inc)
        f = ((time - tp)/self.per - (int)((time-tp)/self.per))*2*np.pi
        if self.inFrontOfStar(f, omega, i):
            d = self.a*np.sqrt(1 - np.sin(omega + f)*np.sin(omega + f)*np.sin(i)*np.sin(i));
            return d
        else:
            return 100 # Return massive number to ensure planet doesn't 'transit' when behind star

    def withinDistanceFromStar(self, time, r):
        d = self.r_skyAtTime(time)
        x_in = max(d - self.rad, 0)
        return not (x_in >= r)

    def isTransiting(self, time):
        return self.withinDistanceFromStar(time, 1)
    
    def skyPosAtTime(self, time):
        tp = self.t_periastron()
        omega = np.deg2rad(self.w)
        i = np.deg2rad(self.inc)
        f = ((time - tp)/self.per - (int)((time-tp)/self.per))*2*np.pi
        X = -self.a*np.cos(f + omega)
        Y = -self.a*np.sin(f + omega) * np.cos(i)
        # Rotate due to spin-orbit misalignment
        theta = -self.beta
        Xfin = X*np.cos(np.deg2rad(theta)) - Y*np.sin(np.deg2rad(theta))
        Yfin = X*np.sin(np.deg2rad(theta)) + Y*np.cos(np.deg2rad(theta))
        return (Xfin,Yfin)