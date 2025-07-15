# -*- coding: utf-8 -*-
import numpy


class Layeritem ():

    def __init__(self, layeritem, pression):
        self.layeritem = layeritem
        self.pression = pression

    def computeLayerLine(self, layeritem=None, pression=None, layers=None):
        """ compute a specific line to represent layeritems """
        if layeritem is not None:
            self.layeritem = layeritem
        else:
            layeritem = self.layeritem
        if pression is not None:
            self.pression = pression
        else:

            pression = self.pression

        if layers is None:
            layers = self.computeLayers(self.pression)
        nlevels = pression.shape[0]
        nlayers = layers.shape[0]
        y = numpy.zeros(nlevels * 2)
        x = numpy.zeros(nlevels * 2)

        x[0] = 0
        y[0] = pression[0]
        y[2] = pression[0]
        x[1] = layeritem[0]

        for i in range(0, nlayers):
            y[2 * i + 1] = pression[i]
            x[2 * i + 1] = layeritem[i]
            y[2 * i + 2] = pression[i + 1]
            x[2 * i + 2] = layeritem[i]
        x[2 * nlevels - 1] = 0
        y[2 * nlevels - 1] = pression[nlevels - 1]
        self.xlayeritem = x
        self.ylayeritem = y
        return (x, y)

    def getLayeritem(self, x=None):
        """ compute layeritem sur layers """
        if x is None:
            x = self.xlayeritem
        else:
            self.xlayeritem = x
        nlayers = x.shape[0] // 2 - 1
        layeritem = numpy.zeros(nlayers)
        for i in range(0, nlayers - 1):
            layeritem[i] = x[i * 2 + 1]
        return layeritem

    def update(self, xlayeritem, ylayeritem):
        self.xlayeritem = xlayeritem
        self.ylayeritem = ylayeritem

    def computeLayers(self, pression=None):
        """ Compute the mean value of pression in a layer """
        if pression is not None:
            self.pression = pression
        foo = numpy.empty(self.pression.shape[0] - 1)
        for i in range(foo.shape[0]):
            foo[i] = (self.pression[i + 1] + self.pression[i]) / 2
        return foo


if __name__ == "__main__":
    import matplotlib
    import rmodel
    import matplotlib.pyplot as plt
    matplotlib.use('WXAgg')

    print("version matplotlib :", matplotlib.__version__)
    p = rmodel.project.Project()

    p.openProfile(p.config.ENV["RTTOV_GUI_PROFILE_DIR"] +
                  "/cldaer101lev_allgas.H5", 1)
    mydata = p.myProfile['STCO']
    print("mydata STCO shape", mydata.shape)
    pression = p.myProfile['P']
    myLayeritem = Layeritem(mydata, pression)
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_yscale('log')
    ax1.set_ylim((pression[-1], 100))
    ax1.set_ylim((pression[-1], 100))
    ax1.set_yticks((100, 200, 300, 500, 1000))
    label = ('100', '200', '300', '500', '1000')
    layers = myLayeritem.computeLayers()
    ax1.set_yticklabels(label)
    for i in range(0, layers.shape[0] - 1):
        ax1.axhline(y=pression[i])
    ax1.plot(mydata, layers, 'g')
    (x, y) = myLayeritem.computeLayerLine()
    ax1.plot(x, y, 'r')
    newx = myLayeritem.getLayeritem(x)
    ax1.plot(newx, layers, 'r')
    print("x", x.shape)
    print(x)
    print("mydata", mydata.shape)
    print(mydata)
    print("newx", newx.shape)
    print(newx)
    for i in range(0, newx.shape[0]):
        if newx[i] != mydata[i]:
            print("pb ici ", i, newx[i])

    print(y)

    plt.show()
