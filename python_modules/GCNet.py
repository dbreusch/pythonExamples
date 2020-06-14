from netCDF4 import Dataset, num2date
from datetime import datetime, timedelta
import numpy as np
import pandas as pd
    
class GCNet:
    """A class for working with GCNet AWS sites"""
    
    def __init__(self, metaFN = None, dataFN = None ):
        self.metaFN = metaFN
        self.dataFN = dataFN
        self.awsMeta = None
        self.awsData = None
        self.time = None

    def loadMeta( self, varName = None ):
        """
           Load standard metadata for all GCNet sites
           or
           load a specific metadata variable.
        """
        if self.metaFN is None:
            print "Warning: Metadata filename is not defined!"
            print "No metadata variables have been created."
        else:
            if varName is None:
                self.awsMeta = Dataset(self.metaFN, "r")
                self.awsLat  = self.awsMeta.variables['site_lat'][:]
                self.awsLon  = self.awsMeta.variables['site_lon'][:] + 360
                self.awsElev = self.awsMeta.variables['site_elev'][:]
                self.awsName = self.awsMeta.variables['site_name'][:]
                self.nSites  = len( self.awsLat )
            else:
                try:
                    X = self.awsMeta.variables[varName][:]
                except:
                    print 'Variable '+varName+' not found.'
                    return None
                return X

    def getMetaDS( self ):
        """ 
           Return the netCDF4 Dataset object for the metadata file.
           Should only be needed in unusual situations!
        """
        
        if self.metaFN is None:
            print "Warning: Metadata filename is not defined!"
            print "Unable to return metadata Dataset variable."
        else:
            return( self.awsMeta )
        
    def convertTime( self, timeCF ):
        """
           *** Superseded by call to netcdf4.num2date() ***
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
    
    def qcStats( self, X ):
        qcPct = []
        nSites = X.shape[1]
        M = X.values
        for ix in (range(nSites)):
            y = M[:,ix].astype("float")
            yAll = len(y)
            yVal = np.count_nonzero(~np.isnan(y))
            yPct = yVal*1.0/yAll*100.
            qcPct.append(yPct)
        return np.array( qcPct )
        
    def awsQC( self, X, varName ):
        if varName == "AirT1":
            # remove anything greater than 20 deg
            Xqc = np.where( X > 20., np.nan, X )
            # if it looks large and has a NaN neighbor, remove it
            Xrow = Xqc.shape[0]
            Xcol = Xqc.shape[1]
            for m in range( Xcol ):
#                 print m
                for n in range( 1, Xrow-1 ):
                    if Xqc[n,m] > 10:
                        if np.isnan( Xqc[n-1,m] ):
                            print "LHS: Resetting value "+str(Xqc[n,m])
                            Xqc[n,m] = np.nan
#                             print "     "+str(Xqc[n,m])
                        if np.isnan( Xqc[n+1,m] ):
                            print "RHS: Resetting value "+str(Xqc[n,m])
                            Xqc[n,m] = np.nan
#                             print "     "+str(Xqc[n,m])
        else:
            Xqc = None
        return Xqc
            
    def loadData( self, varName, dFrame = None ):
        """
           Load a particular variable
           If it's time, return **both** a DateTime array and the original
           CF-compliant variable.
        """
        
        if self.awsData is None:
            if self.dataFN is None:
                print "Data filename not defined!"
                return None
            else:
                self.awsData = Dataset(self.dataFN, "r")
        X = self.awsData.variables[varName]
        if varName == "time":  # save converted time locally
#             self.time = self.convertTime( X )
            print X.units
            print X.calendar
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
                Yqc = self.awsQC( Y, varName )
                Yqc = np.where( Yqc > 1.e36, np.nan, Yqc )

                try:
                    M = np.concatenate( [self.time[:,None], Yqc], axis=1 )
                    awsNames = ["Time"]
                except:
                    M = Y
                    awsNames = []
                for awsIX in range(self.nSites):
                    awsNames.append( self.getName( awsIX ) )
                df = pd.DataFrame( M, columns=awsNames )
                return df

    def getDataDS( self ):
        """ 
           Return the netCDF4 Dataset object for the data file.
           Should only be needed in unusual situations!
        """
        if self.DataFN is None:
            print "Warning: Datadata filename is not defined!"
            print "Unable to return data Dataset variable."
        else:
            return( self.awsData )
        
    def info( self, awsIX ):         
        """Print metadata for specific AWS"""
        name = ''.join(self.awsName[awsIX,:])
        awsStr = str(awsIX+1)
        print "Site "+awsStr+": "+name 
        print "Lat, Lon:", self.awsLat[awsIX], self.awsLon[awsIX]
        print "Elevation:", self.awsElev[awsIX] 

    def getAWS( self, awsIX ):
        """Return metadata for specific AWS"""
        name = self.getName(awsIX)
        return ( self.awsLat[awsIX], self.awsLon[awsIX], self.awsElev[awsIX], name )
    
    def getNSites( self ):
        return self.nSites
    
    def getTime( self, timeCF = None ):
        """Return time for all AWS"""
        if timeCF is None:
            if self.time is None:
                print "Must call loadData() for time before retrieving it!"
                return none
            else:
                return ( self.time )
        else:
            self.time = self.convertTime( timeCF )
            return( self.time )

    def getLat( self, awsIX ):
        """Return latitude for specific AWS"""
        return ( self.awsLat[awsIX] )
    
    def getLon( self, awsIX ):
        """Return longitude for specific AWS"""
        return ( self.awsLon[awsIX] )

    def getElev( self, awsIX ):
        """Return elevation for specific AWS"""
        return ( self.awsElev[awsIX] )

    def getName( self, awsIX ):
        """Return name for specific AWS"""
        name = ''.join(self.awsName[awsIX,:]).replace('_',' ')
        return ( name )
    
    def getNames( self ):
        """Get all site names"""
        awsNames = []
        for awsIX in range(self.nSites):
            awsNames.append( self.getName( awsIX ) )
        return awsNames