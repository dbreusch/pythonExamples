from netCDF4 import Dataset, num2date
import numpy as np
from geopy.distance import great_circle
from scipy import interpolate
from mpl_toolkits.basemap import Basemap
import math
from datetime import datetime, timedelta
import pandas as pd

class Model:
    """A class for model objects"""
    
    def __init__(self, metaFN = None, dataFN = None):
        self.metaFN = metaFN
        self.dataFN = dataFN
        self.gData = None
        self.mGrid = None
        self.modelMeta = None
        self.modelData = None
        self.time = None
#         print self.metaFN, self.dataFN
        
    def setDomain( self, corners ):
        self.minLat = corners[0]
        self.maxLat = corners[1]
        self.minLon = corners[2]
        self.maxLon = corners[3]
        self.midLat = (self.minLat + self.maxLat) / 2
        self.midLon = (self.minLon + self.maxLon) / 2
        
    def setGrid( self ):
        self.mGrid = Basemap(projection='stere', \
                     width=2400000, height=3300000, \
                     lat_ts=self.midLat, lat_0=self.midLat, lon_0=self.midLon, \
                     resolution=None )
        
    def loadMeta(self):
        if self.metaFN is None:
            print "Warning: Metadata filename is not defined!"
            print "No metadata variables have been created."
        else:
            self.modelMeta = Dataset( self.metaFN, "r")

            self.lat = self.modelMeta.variables['lat'][:]
            self.lon = self.modelMeta.variables['lon'][:]
            self.lon = np.where( self.lon < 0, self.lon+360, self.lon )
            self.elev = self.modelMeta.variables['orog'][:]

            ndim = len(self.lat.shape)
            if ndim > 1:  # 2-d lat/lon
                print "2d dataset"
                self.lat2d = self.lat
                self.lat = self.lat2d.flatten()
                self.nlat = self.lat2d.shape[0]
                
                self.lon2d = self.lon
                self.lon = self.lon2d.flatten()
                self.nlon = self.lon2d.shape[1]
                
#                 self.elev2d = self.elev
#                 self.elev = self.elev2d.flatten()
                
                self.npts = len( self.lat )
                
# 6/5/18: taking out dx,dy because I don't use it for erai/cesm and it doesn't
#         mean anything for wrf (when calculated like this).
#             self.dy = self.lat[1] - self.lat[0]
#             if self.dy < 0:
#                 self.dy = self.dy * -1
#             self.dx = self.lon[1] - self.lon[0]
#             if self.dx < 0:
#                 self.dx = self.dx * -1
            
            else:
                print "1d dataset"
                self.nlat = len(self.lat)
                self.nlon = len(self.lon)
                self.npts = self.nlat * self.nlon
                z = np.zeros(self.elev.shape)
                self.lat2d = np.transpose(np.transpose(z) + self.lat)
                self.lon2d = z + self.lon
            
#             print ndim, self.nlat, self.nlon, self.npts

            self.ixLat = np.array(range( self.nlat ),dtype=np.uint16)
            self.ixLon = np.array(range( self.nlon ),dtype=np.uint16)

            z = np.zeros(self.elev.shape,dtype=np.uint16)
            self.ixLon2d = z + self.ixLon
            self.ixLat2d = np.transpose(np.transpose(z) + self.ixLat)

            self.metaDtype = np.dtype( [ ('ixLat',self.ixLat2d.dtype), ('lat',self.lat2d.dtype), \
                            ('ixLon',self.ixLon2d.dtype), ('lon',self.lon2d.dtype), \
                            ('distance',np.float) ] )

            self.gData = np.zeros( self.npts, dtype = self.metaDtype )
            self.gData['ixLat'] = self.ixLat2d.reshape( self.npts )
            self.gData['ixLon'] = self.ixLon2d.reshape( self.npts )
            self.gData['lat'] = self.lat2d.reshape( self.npts )
            self.gData['lon'] = self.lon2d.reshape( self.npts )
            
            
    def convertTime( self, timeCF ):
        """
           *** Superseded by call to netcdf4.num2str() ***
           Convert CF-compliant time to Python datetime.
           NOTE: for unknown reasons, timedelta only accepts numpy.float64!
        """
        tUnitStr = timeCF.units
        tUnits = tUnitStr.split()[0]
        tBase = ' '.join(tUnitStr.split()[-2:])
        tStart = datetime.strptime( tBase,"%Y-%m-%d %H:%M:%S")
        if tUnits == "hours":
            try:
                time = np.array( [tStart+timedelta(hours=x) for x in timeCF ])
            except TypeError:
                temp1 = np.array(timeCF).astype('float64')
                time = np.array( [tStart+timedelta(hours=x) for x in temp1 ])
            return time
        else:
            if tUnits == "days":
                try:
                    time = np.array( [tStart+timedelta(days=x) for x in timeCF ])
                except TypeError:
                    temp1 = np.array(timeCF).astype('float64')
                    time = np.array( [tStart+timedelta(days=x) for x in temp1 ])
                return time
            else:
                print 'Time units "'+tUnits+'" not recognized!'
                return None
    
    def loadData( self, varName, dFrame = None, mbr = 1 ):
        """
           Load a particular variable
           If it's time, return **both** a DateTime array and the original
           CF-compliant variable.
        """
        
        if self.modelData is None:
            if self.dataFN is None:
                print "Data filename not defined!"
                return None
            else:
                self.modelData = Dataset(self.dataFN, "r")
        X = self.modelData.variables[varName]
        if varName == "time":  # save converted time locally
#             self.time = self.convertTime( X )
            self.time = num2date(X[:], X.units, X.calendar)
            if dFrame is None:
                vtype = type(self.time[0])
                # noleap calendars need special handling to convert from "netcdftime" objects
                # back to normal "datetime" objects (otherwise matplotlib fails)
                if "netcdftime" in str(type(self.time[0])):
                    d = [ datetime( x.year, x.month, x.day, x.hour, x.minute) for x in self.time ]
                    self.time = np.array( d )
                return ( self.time, X )
            else:
                s = pd.Series( self.time, name="Time" )
                return s
        else:
            if dFrame is None:
                return X
            else:
                Y = np.array( X )
                Y = np.where( Y > 1.e36, np.nan, Y )
                if varName == "tas":
                    if np.min(Y) > 100:
                        Y = Y - 273.15
                ndims = Y.shape
                ncols = ndims[-1]
                if len(ndims) > 2:  
                    if len(ndims) == 3: # 3-d var, i.e., ensemble members
                        nens, ntimes, nstns = Y.shape
                        if mbr < 0:
                            if mbr == -1:
                                # calculate ensemble average
                                Y = np.average( Y, axis=0)
                                M = np.concatenate( [self.time[:,None], Y], axis=1 )
                            else:
                                if mbr == -2:
                                    # return the whole ensemble reshaped
                                    Y = np.reshape( Y, (nens*ntimes, nstns), order='C')
                                    try:
                                        timeRep = np.tile( self.time, nens )
                                        M = np.concatenate( [timeRep[:,None], Y], axis=1 )
                                    except:
                                        M = Y
                                else:
                                    M = Y
                        else:
                            mbrIX = mbr - 1
                            if mbrIX < 0:
                                mbrIX = 0
                            if mbrIX > (ncols-1):
                                mbrIX = ncols - 1
                            Y = Y[mbrIX,:,:]
                            try:
                                M = np.concatenate( [self.time[:,None], Y], axis=1 )
                            except:
                                M = Y
                    else:
                        M = Y  # 4-d or higher
                else:
                    if len(ndims) > 1:  # 2-d
                        try:
                            M = np.concatenate( [self.time[:,None], Y], axis=1 )
                        except:
                            M = Y
                    else:
                        M = Y  # 1-d
                df = pd.DataFrame( M )
                return df        

    def getLat2d(self, y=None, x=None):
        if any( ( y is None, x is None ) ): 
            return self.lat2d
        else:
            return self.lat2d[ y, x]
    
    def getLon2d(self, y=None, x=None):
        if any( ( y is None, x is None ) ): 
            return self.lon2d
        else:
            return self.lon2d[ y, x]
    
    def getIXLat2d(self):
        return self.ixLat2d
    
    def getIXLon2d(self):
        return self.ixLon2d

    def getElev(self, y=None, x=None):
        if any( ( y is None, x is None ) ): 
            return self.elev
        else:
            return self.elev[ y, x]

    def getNPts(self):
        return( self.npts )

#     def getDY(self):
#         return( self.dy )

#     def getDX(self):
#         return( self.dx )

    def info(self):
        print "nlat, nlon, npts:",self.nlat, self.nlon, self.npts
        print "lat.shape:",self.lat.shape
        print "lon.shape:",self.lon.shape
        print "elev.shape:",self.elev.shape
        
    def getMeta(self):
        return self.gData

    def closestPoints( self, siteLatLon ):
        if self.gData is None:
            self.loadMeta()
        distList = [ ( great_circle( siteLatLon, ( x['lat'], x['lon'] ) ).km ) for x in self.gData ]
        self.gData['distance'] = distList
        nearest = np.sort( self.gData, order = 'distance')
        
        """ four nearest surrounding grid points
        this is tricky!  above certain latitudes the gridpoints are so close along the
        latitude that four closest may all be on same latitude (which crashes the
        interpolator.  in other cases, can end up with a point not on the box around the
        station (fourth point is closer than the last point on the box).  so I ended up
        with a dictionary-based approach that only adds a point if its latitude and
        longitude haven't already appeared twice.
        Note that for WRF, all of this complexity is irrelevant since it's not on a
        rectangular lat/lon grid like erai/cesm, so first four closest are the best!
        """
        Cl =  nearest.tolist()
        C = [Cl[0]]  # absolute closest is first element on the list...
        latDict = {}
        lonDict = {}
        latDict[ str( C[0][1] ) ] = 1
        lonDict[ str( C[0][3] ) ] = 1
        ix = 1       # ...so start from second element when searching for other points!
        while len(C) < 4:
            latx = str( Cl[ix][1] )
            lonx = str( Cl[ix][3] )
            if latx not in latDict:     # found new latitude
                if lonx not in lonDict:
                    C.append(Cl[ix])
                    latDict[ latx ] = 1
                    lonDict[ lonx ] = 1
                else:
                    if lonDict[ lonx ] < 2:
                        C.append(Cl[ix])
                        latDict[ latx ] = 1
                        lonDict[ lonx ] += 1
            else:
                if latDict[ latx ] < 2:
                    if lonx not in lonDict:
                        C.append(Cl[ix])
                        latDict[ latx ] += 1
                        lonDict[ lonx ] = 1
                    else:
                        if lonDict[ lonx ] < 2:
                            latDict[ latx ] += 1
                            lonDict[ lonx ] += 1
                            C.append(Cl[ix])
            ix += 1
            if ix > 30:
                break

        """ convert list back to a structured array """
        closest = np.zeros( 4, dtype=nearest.dtype )
        closest['ixLat'] = [ L[0] for L in C ]
        closest['lat'] = [ L[1] for L in C ]
        closest['ixLon'] = [ L[2] for L in C ]
        closest['lon'] = [ L[3] for L in C ]
        closest['distance'] = [ L[4] for L in C ]
        
        return closest
    
    def interpolate( self, V, closest, siteLatLon ):
        """
        interpolate data at four closest grid points to site location.
        to avoid lat/lon-related artifacts, this is done on lat/lon converted to grid x,y
        """
        
        """ point to interpolate to """
        siteLat = siteLatLon[0]
        siteLon = siteLatLon[1]

        """ convert lat/lon to map projection grid coordinates """
        if self.mGrid is None:
            self.setGrid()   # create map projection (using Basemap) if it doesn't exist yet
        cX, cY = self.mGrid( closest['lon'], closest['lat'] )
        siteXY = self.mGrid( siteLon, siteLat )

        """ do the interpolation """
        if len( V ) == 4:
            Vnew = None
            try:
                Vnew = interpolate.griddata((cX, cY), V, siteXY, method='linear')
            except ValueError, e:
                print "Problem with interpolation"
                print e.args
        else:
            nRec = len( V )
            Vnew = np.empty( nRec, dtype = 'float' )
            Vnew[:] = None
            for r in range( nRec ):
                try:
                    Vnew[r] = interpolate.griddata((cX, cY), V[r], siteXY, method='linear')
                except ValueError, e:
                    print "Problem with interpolation"
                    print e.args

        return Vnew