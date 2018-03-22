class Params(object):
    def __init__(self):
        # Stellar Params
        self.prot = None # Star's Rotational Period [days]
        self.sinc = None # Star's inclination [degrees]
        self.u = None # Quad limb darkening coeffs [u1, u2]
        self.rad_star = None # Radius [RSol]
        self.res = None # Resolution of LD grid
        self.high_band = None # Highest active band [deg]
        self.low_band = None # Lowest active band [deg]
        self.stellar_cycle = None # Used to calculate drift in active latitudes [days]
        self.spot_asymmetry = None # [0, 1] if > 0.5 then in southern hemisphere
        
        # Planet Params
        self.t0 = None # Time of inferior conjunction [days]
        self.per = None # Orbital period [days]
        self.rp = None # Radius of planet [RStar]
        self.a = None # Semi-major axis [RStar]
        self.oinc = None # Orbital inclination [deg]
        self.w = None # Argument of periapsis [deg]
        self.beta = None # Sky spin-orbit alignment [deg] (lambda is reserved in python)
        
        # Spot Params (Min - Max values)
        self.b = None # Realtive intensity to stellar average
        self.spot_rads = None # Radii of spots [RStar]
        self.spot_lives = None # How long spots live on surface [days]
        self.nspots = None # Number of spots on surface at once
        
        # Integrator/System Params
        self.rotation_dt = None # Time step for general rotation curve [days]
        self.transit_samples = None # Time step during transit [days]
        self.adjustDistance = None # Distance from star to change time step [RStar]