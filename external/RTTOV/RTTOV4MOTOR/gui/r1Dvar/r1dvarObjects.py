'''
Created on Jan 19, 2015

@author: pascale
'''
import rttov
import numpy as np


class R1dvarKmatrix(rttov.kmatrix.Kmatrix):
    '''
    Class derived from Kmatrix
    Has ist own methods for creatif item need by the 1d var
    '''

    def __init__(self):
        '''
        '''
        rttov.kmatrix.Kmatrix.__init__(self)
        self.nlevels = None

    def toKmatrix1dvar(self):
        """ compute the Kmatrix for safnwp 1dvar """
        """ the Kmatrix is a matrix nchannel x vector1d """
        self.nlevels = self.profile["P"].shape[0]

        if isinstance(self.chanprof['CHANNELS'], np.int32):
            nchannels = 1
        else:
            nchannels = len(self.chanprof['CHANNELS'])

        kmatrix1dvar = np.zeros((nchannels, self.nlevels + 29 + 3), float)

        for channel in range(nchannels):
            kprofile = self.getchanprof(channel)
            vector = kprofile.to1DvarVect()

            kmatrix1dvar[channel, :] = vector

        return kmatrix1dvar


class R1dvarRadiance(rttov.radiance.Radiance):

    def __init__(self):
        '''
        '''
        rttov.radiance.Radiance.__init__(self)

    def toObsfile(self, filename):
        """ create an Obsfile.dat wich can be used bu the sanwp
            1dvar software as input """
        """ see https://nwpsaf.eu/deliverables/nwpsaf_1dvar/
            nwpsaf-mo-ud-032_NWPSAF_1DVar_Manual.html#FilesIn """
        """ can't be used for atovs because rttov can deal with only
             one instrument """

        if self.misc is None:
            raise RuntimeError(
                "cannot create an 1dvar obsfile.dat from Radiance Object")
        f = open(filename, "w")
        f.write("This is a simulated observation dataset for " +
                self.misc["INSTRUMENT"] + "\n")
        f.write("Generated from h5 file created by RTTOV version 11.2" + "\n")
        for l in range(0, 8):
            f.write("\n")
        f.write("Number of Observations in File:     1" + "\n")
        f.write("No. of Chans per Observation:" +
                ("%d" % self.misc["NCHANNELS"]).rjust(8) + "\n")
        f.write("Number of instruments making up observations : 1" + "\n")
        f.write(
            "*** In the following Series, Platform and Instrument are "
            "defined  ***" + "\n")
        f.write(
            "*** according to the relevant RT Model definitions "
            "(if required): ***" + "\n")
        f.write(
            "Sat. Series   Platform   Instrument First_Channel   "
            "Last_Channel  Sat ID" + "\n")
        f.write(
            "10            2          16           1            "
            "8461           4" + "\n")
        f.write("Channels:" + "\n")
        chan = 0
        while chan < self.misc["NCHANNELS"]:
            for k in range(0, 16):
                chan = chan + 1
                if chan > self.misc["NCHANNELS"]:
                    continue
                f.write(("%d" % (chan)).rjust(5))
            f.write("\n")
        f.write(
            "-------------------------------------------------------"
            "---------------\n")
        f.write(
            "Obs ID:              1 Obs Type:          "
            "3 Satellite ID:     4" + "\n")
        f.write("Latitude:  -90.000 Longitude:    0.000 "
                "Elevation:    0.0" + "\n")
        f.write(
            "Surface Type:   1 Sat Zen Angle:    0.000 Solar Zen."
            " Ang.:    0.000" + "\n")
        f.write("Brightness Temperatures:" + "\n")
        chan = 0
        while chan < self.misc["NCHANNELS"]:
            for k in range(0, 6):
                chan = chan + 1
                if chan > self.misc["NCHANNELS"]:
                    continue
                f.write(("%6.3f" % (self["BT"][chan - 1])).rjust(13))
            f.write("\n")

        f.close()
