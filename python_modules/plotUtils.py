import numpy as np
import pandas as pd
import seaborn as sns
from collections import OrderedDict

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D
from matplotlib.text import Text

# from netCDF4 import Dataset, num2date
# from geopy.distance import great_circle
# from scipy import interpolate
# from mpl_toolkits.basemap import Basemap
# import math
# from datetime import datetime, timedelta

class PlotUtils:
    """
       A class for plotting.
    """
    
    def __init__( self, plotOptions ):
        """
           Class constructor.
           Supports MANY default values.
           
           saveFigure is a GLOBAL option to save/not save figures to output files
        """
        self.fig = None
        self.ax = None
        self.dfData = None
        self.dfData1 = None
        self.dfData2 = None

        # set up default options dictionary
        #  var: variable to plot
        #  saveFigure: make a plot of the figure (boolean)
        #  type: type of plot (all, jja, jul, ovr)
        #  figWidth, figHeight: figure size
        #  yTickInterval: y-axis tick interval
        #  xLabel, yLabel: axis labels
        #  labelFontSize: axis label size
        #  tickFontSize: axis tick mark size
        #  inner: inner plot for makePlots (quartile, None)
        #  inner1, inner2: inner plots for makePlots2 (quartile, None); don't really expect to use inner1!
        #  cutLimit, bandWidth, scale, style, palette, context: violinplot options
        #  colorLHS1, colorLHS2, colorRHS1, colorRHS2: colors for each side of plot and each violinplot
        # legLoc: main legend location
                                    
        self.paramDict = dict( [ \
                                   ('var', "tas"), \
                                   ('saveFigure', False), \
                                   ('type', "all"), \
                                   ('figWidth', 16), \
                                   ('figHeight', 10), \
                                   ('yTickInterval', 5), \
                                   ('xLabel', "GCNet AWS"), \
                                   ('yLabel', "tas ($^\circ$C)"), \
                                   ('inner', "quartile"), \
                                   ('inner1', None), \
                                   ('inner2', "quartile"), \
                                   ('cutLimit', 0.0), \
                                   ('bandWidth', 0.1), \
                                   ('scale', "count"), \
                                   ('colorLHS1', "b"), \
                                   ('colorLHS2', "g"), \
                                   ('colorRHS1', "y"), \
                                   ('colorRHS2', "r"), \
                                   ('style', "whitegrid"), \
                                   ('palette', "pastel"), \
                                   ('context', "talk"), \
                                   ('labelFontSize', 22), \
                                   ('tickFontSize', 18), \
                                   ('legLoc', 'lower right'), \
                                   ('legLoc2', 'upper left') \
                               ] )
#         self.paramDict = dict( zip( self.paramKW, self.paramDef) )
        self.paramKW = self.paramDict.keys()
        
        # look for keywords in passed options and override defaults as found
        for k in self.paramKW:
            if k in plotOptions.keys():
                self.paramDict[k] = plotOptions[k]
                
        # common y-axis values, by plot "type"
        self.yAxisValues = dict( [ \
                                    ('all', ( -70, 15 ) ), \
                                    ('jja', ( -35, 15 ) ), \
                                    ('jun', ( -30, 20 ) ), \
                                    ('jul', ( -30, 20 ) ), \
                                    ('aug', ( -30, 20 ) ), \
                                    ('ovr', ( -27, 20 ) ) \
                               ] )
            
    """
    Utility function to set option value(s) with a dictionary
    """
    def setOption( self, newOpt ):
        for k in self.paramKW:
            if k in newOpt.keys():
                self.paramDict[k] = newOpt[k]

    """
    Utility function to get an option value
    """
    def getOption( self, optKW ):
        if optKW in self.paramKW:
            return self.paramDict[optKW]
        else:
            print "Option keyword "+optKW+" not found"
            return None
        
    """
    Utility function for xlabel
    """
    def setXlabel( self, xlabel="GCNet AWS" ):
        self.xLabel = xlabel

    """
    Utility function for ylabel
    """
    def setYlabel( self, ylabel="tas ($^\circ$C)" ):
        self.yLabel = ylabel

    """
    Build a dataframe to support plotting
    """
    def buildDF( self, X, Y, names, modelName1, modelName2 ):
        nSite = X.shape[1]-1
        for n in range(nSite):
            df1 = pd.DataFrame(columns=[self.paramDict["var"],'Source','Site','month'])
            df1[self.paramDict["var"]] = X[names[n]].astype("float")
            df1['Source'] = modelName1
            df1['Site'] = names[n]
            df1['month'] = X['month']

            df2 = pd.DataFrame(columns=[self.paramDict["var"],'Source','Site','month'])
            df2[self.paramDict["var"]] = Y[names[n]].astype("float")
            df2['Source'] = modelName2
            df2['Site'] = names[n]
            df2['month'] = Y['month']

            if n == 0:
                dfData = pd.concat( [df1, df2] )
            else:
                dfData = pd.concat( [dfData, df1, df2] )

        dfData.dropna(inplace=True)
        return dfData

    """
    Standard seaborn split violin plots, one dataset on each side of figure
    """
    def makePlots( self, X, Y, titleStr, modelName1, modelName2 ):

        # replace spaces in modelName (allows for flexible naming)
        mName1 = modelName1.replace(' ','_')
        mName2 = modelName2.replace(' ','_')

        # first build a dataframe for the plotting data
        awsList = X.columns.values[:-1]
        self.dfData = self.buildDF( X, Y, awsList, mName1, mName2 )

        # now do the plots
        sns.set(style=self.paramDict["style"], palette=self.paramDict["palette"], \
                context = self.paramDict["context"], color_codes=True)
        self.fig, self.ax = plt.subplots(figsize=(self.paramDict["figWidth"], self.paramDict["figHeight"]))

        sns.violinplot(x="Site", y=self.paramDict["var"], hue="Source", data=self.dfData, \
                       split=True, cut=self.paramDict["cutLimit"], \
                       bw=self.paramDict["bandWidth"], scale=self.paramDict["scale"], scale_hue=True, \
                       palette={mName1: self.paramDict["colorLHS1"], mName2: self.paramDict["colorRHS1"]}, \
                       linewidth=3, legend=False, inner=self.paramDict["inner"])

        # set y axis range to common values
        try:
            yLimits = self.yAxisValues[ self.paramDict["type"].lower() ]
            plt.ylim( yLimits[0], yLimits[1] )
            print "y-axis limits: "+str(yLimits[0])+", "+str(yLimits[1])
        except:
            print "Plot type "+self.paramDict["type"]+" not in yAxisValues table"
            print "Skipping setting of y axis limits"
        
        # are there inner quartile lines?
        if self.paramDict["inner"] == "quartile":
            # change quartile line characteristics
            k=0
            for o in self.fig.findobj(Line2D):
                if o.get_linestyle() == '--':
                    o.set_color('red')
                    o.set_linewidth(2)
                    if k == 1:
                        o.set_linestyle('-')
                    else:
                        o.set_linestyle('dotted')
                    k = (k + 1) % 3

            # custom legend for quartile lines
            lineArtist = Line2D((0,1),(0,0), color='red', linestyle='-', linewidth=2)
            dotArtist  = Line2D((0,1),(0,0), color='red', linestyle='dotted', linewidth=2)
            legend1 = plt.legend( [lineArtist, dotArtist], ['Mean', 'Quartile'], \
                                   loc=self.paramDict["legLoc2"], frameon=True, \
                                   fancybox=True, edgecolor='black', fontsize=18)
            self.ax.add_artist(legend1)

        # tweak the plots and add a title
        sns.despine(offset=10)
        plt.xticks(rotation=45, fontsize=self.paramDict["tickFontSize"] )
        plt.yticks(fontsize=self.paramDict["tickFontSize"])
        self.ax.set_xlabel(self.paramDict["xLabel"], size = self.paramDict["labelFontSize"], alpha=0.7)
        self.ax.set_ylabel(self.paramDict["yLabel"], size = self.paramDict["labelFontSize"], alpha=0.7);
        self.ax.yaxis.set_major_locator(ticker.MultipleLocator(self.paramDict["yTickInterval"]))
        self.fig.suptitle(titleStr, fontsize=18, fontweight='bold');

        # add legend
        legend2 = plt.legend(loc=self.paramDict["legLoc"],ncol=2, frameon=True,  \
                             fancybox=True, edgecolor='black', fontsize=20)
        self.ax.add_artist(legend2)

        # fix legend labels
        for o in self.fig.findobj(Text):
            txt = o.get_text()
            if '_' in txt:
                txt2 = txt.replace('_',' ')
                o.set_text( txt2 )

        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        

    """
    Complex seaborn split violin plots, one or two datasets on each side of figure
    """
    def makePlots2( self, X1, X2, Y1, Y2, titleStr, modelName1, modelName2, modelSuff ):
        '''
        X1: LHS model data, full dataset
        X2: LHS model data, ensemble average (if None, then set to X1, e.g., for AWS)
        Y1: RHS model data, full dataset
        Y2: RHS model data, ensemble average
        '''

        # replace spaces in modelName (allows for flexible naming)
        mName1a = modelName1.replace(' ','_')+"_"+modelSuff[0]
        mName1b = modelName1.replace(' ','_')+"_"+modelSuff[1]
        mName2a = modelName2.replace(' ','_')+"_"+modelSuff[0]
        mName2b = modelName2.replace(' ','_')+"_"+modelSuff[1]

        color1a = self.paramDict["colorLHS1"]
        color1b = self.paramDict["colorLHS2"]
        nCol = 2
        LHScount = 2
        RHScount = 2
        if X2 is None:
            LHScount = 1
            X2 = X1 # allow for X being the same in both plots, e.g., for AWS comparisons
            mName1a = modelName1.replace(' ','_')
            mName1b = mName1a
            color1b = color1a
            nCol = 3
        color2a = self.paramDict["colorRHS1"]
        color2b = self.paramDict["colorRHS2"]
        if Y2 is None:
            RHScount = 1
            Y2 = Y1 # allow for Y being the same in both plots, e.g., for CESM LW comparisons
            mName2a = modelName2.replace(' ','_')
            mName2b = mName2a
            color2b = color2a
            nCol = 3 

        # first build a dataframes for the plotting data
        awsList = X1.columns.values[:-1]
        self.dfData1 = self.buildDF( X1, Y1, awsList, mName1a, mName2a )
        self.dfData2 = self.buildDF( X2, Y2, awsList, mName1b, mName2b )

        # now do the plots
        sns.set(style=self.paramDict["style"], palette=self.paramDict["palette"], color_codes=True)
        sns.set_context(self.paramDict["context"])
        self.fig, self.ax = plt.subplots(figsize=(self.paramDict["figWidth"], self.paramDict["figHeight"]))

        # first dataset (usually larger plots)
        sns.violinplot(x="Site", y=self.paramDict["var"], hue="Source", data=self.dfData1, \
                       split=True, cut=self.paramDict["cutLimit"], ax = self.ax, \
                       bw=self.paramDict["bandWidth"], scale=self.paramDict["scale"], scale_hue=True, \
                       palette={mName1a: color1a, mName2a: color2a}, \
                       linewidth=3, legend=False, inner=self.paramDict["inner1"])

        # second dataset (usually smaller plots)
        sns.violinplot(x="Site", y=self.paramDict["var"], hue="Source", data=self.dfData2, \
                       split=True, cut=self.paramDict["cutLimit"], ax = self.ax, \
                       bw=self.paramDict["bandWidth"], scale=self.paramDict["scale"], scale_hue=True, \
                       palette={mName1b: color1b, mName2b: color2b}, \
                       linewidth=3, legend=False, inner=self.paramDict["inner2"])

        # tweak the plots and add a title
        sns.despine(offset=10)
        plt.xticks(rotation=45, fontsize=self.paramDict["tickFontSize"])
        plt.yticks(fontsize=self.paramDict["tickFontSize"])
        self.ax.set_xlabel(self.paramDict["xLabel"], size = self.paramDict["labelFontSize"], alpha=0.7)
        self.ax.set_ylabel(self.paramDict["yLabel"], size = self.paramDict["labelFontSize"], alpha=0.7);
        self.ax.yaxis.set_major_locator(ticker.MultipleLocator(self.paramDict["yTickInterval"]))
        self.fig.suptitle(titleStr, fontsize=18, fontweight='bold');

        # set y axis range to common values
        try:
            yLimits = self.yAxisValues[ self.paramDict["type"].lower() ]
            plt.ylim( yLimits[0], yLimits[1] )
            print "y-axis limits: "+str(yLimits[0])+", "+str(yLimits[1])
        except:
            print "Plot type "+self.paramDict["type"]+" not in yAxisValues table"
            print "Skipping setting of y axis limits"
            
         # are there inner quartile lines?
        if self.paramDict["inner2"] == "quartile":
            # change quartile line characteristics
            k=0
            for o in self.fig.findobj(Line2D):
                if o.get_linestyle() == '--':
                    o.set_color('red')
                    o.set_linewidth(2)
                    if k == 1:
                        o.set_linestyle('-')
                    else:
                        o.set_linestyle('dotted')
                    k = (k + 1) % 3

            # custom legend for quartile lines
            lineArtist = Line2D((0,1),(0,0), color='red', linestyle='-', linewidth=2)
            dotArtist  = Line2D((0,1),(0,0), color='red', linestyle='dotted', linewidth=2)
            legend1 = plt.legend( [lineArtist, dotArtist], ['Mean', 'Quartile'], \
                                   loc=self.paramDict["legLoc2"], frameon=True, \
                                   fancybox=True, edgecolor='black', fontsize=18)
            self.ax.add_artist(legend1)

        # add legend
        handles, labels = plt.gca().get_legend_handles_labels()
        if nCol == 2: # need to reorder
            order = [0,2,1,3]
            legend2 = plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order], \
                                  loc=self.paramDict["legLoc"], ncol=nCol, \
                                  frameon=True, fancybox=True, edgecolor='black', fontsize=18)
        else:  # need to remove duplicates (for LHS)
            if LHScount == 1:
                by_label = OrderedDict(zip(labels, handles))
                legend2 = plt.legend(by_label.values(), by_label.keys(), \
                                     loc=self.paramDict["legLoc"], ncol=nCol, \
                                     frameon=True, fancybox=True, edgecolor='black', fontsize=18)
            else:
                if RHScount == 1:
                    ll = zip(labels, handles).sort()
                    order = [ 0, 2, 1 ]
                    legend2 = plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order], \
                                         loc=self.paramDict["legLoc"], ncol=nCol, \
                                         frameon=True, fancybox=True, edgecolor='black', fontsize=18)
        self.ax.add_artist(legend2)
        
        # fix legend labels
        for o in self.fig.findobj(Text):
            txt = o.get_text()
            if '_' in txt:
                txt2 = txt.replace('_',' ')
                o.set_text( txt2 )

        plt.tight_layout(rect=[0, 0.03, 1, 0.95])

    """
    Save figure to a file
    """
    def saveFigure( self, modelName, varName, srcFile, yrRange, enss, monStr = None, fmt = "eps"):
        if self.paramDict["saveFigure"]:
            # build up a filename string from variable components
            if srcFile is not None:
                srcFile = srcFile.lower()
            else:
                srcFile = "closeinterp"

            fileStr = modelName.lower()
            if enss is not None:
                fileStr = fileStr+"_"+enss
            fileStr = fileStr+"_"+str(yrRange[0])+"-"+str(yrRange[1])+"_"+monStr.lower()+"_"+srcFile
            fileStr = fileStr.replace(' ','_')
            ofn = fileStr+"."+fmt
            print "Saving figure to "+ofn
            plt.savefig(ofn, format=fmt)
        else:
            pass
    
        """
    Save figure to a file (provide output file name)
    """
    def saveFigure2( self, fileStr, fmt = "eps"):
        if self.paramDict["saveFigure"]:
            ofn = fileStr+"."+fmt
            print "Saving figure to "+ofn
            plt.savefig(ofn, format=fmt)
        else:
            pass