import numpy as np
import pyqtgraph as pg
import pyqtgraph.opengl as gl
from PyQt5 import QtGui, QtCore

class Sibyl2DViewer(gl.GLViewWidget):
    App = None
    can_idx = []
    can_clr = []
    def __init__(self, plot, app=None):
        if self.App is None:
            if app is not None:
                self.App = app
            else:
                self.App = QtGui.QApplication([])
        super(Sibyl2DViewer,self).__init__()
        self.Plot = plot
        self.addItem(self.Plot)
        self._downpos = []
        self._init_camera()

    # Fix for "high dpi":
    def width(self):
        trueWidth = super(Sibyl2DViewer,self).width()
        pixRatio = super(Sibyl2DViewer,self).devicePixelRatio()
        return trueWidth * pixRatio

    def height(self):
        trueHeight = super(Sibyl2DViewer,self).height()
        pixRatio = super(Sibyl2DViewer,self).devicePixelRatio()
        return trueHeight * pixRatio

    def _init_camera(self):
        # init camera -- this is pseudo 2D
        self.opts['center'] = QtGui.QVector3D(0, 0, 0)
        self.setCameraPosition(pos=np.array([0, 0, 0]))
        self.opts['elevation'] = 90
        self.opts['distance'] = 35000
        self.opts['azimuth'] = -90 
        self.Plot.phi = 0
        self.Plot.update()

    # All keyboard and mouse event handlers
    def keyPressEvent(self, ev):
        super(Sibyl2DViewer,self).keyPressEvent(ev)
        if ev.key() == ord(' '):
            self._init_camera()

    def mousePressEvent(self, ev):
        ''' Store the position of the mouse press for later use '''
        super(Sibyl2DViewer, self).mousePressEvent(ev)
        self._downpos = self.mousePos

    def mouseReleaseEvent(self, ev):
        ''' Allow for single click to move and right click for context menu '''
        super(Sibyl2DViewer, self).mouseReleaseEvent(ev)
        if self._downpos == ev.pos():
            x = ev.pos().x()
            y = ev.pos().y()
            if ev.button() == 2:
                self.mPosition()
            elif ev.button() == 1:
                x = x - self.width() / 2
                y = y - self.height() / 2
        self._prev_zoom_pos = None
        self._prev_pan_pos = None

    def mouseMoveEvent(self, ev):
        pos = ev.pos().x(), ev.pos().y()
        if not hasattr(self, '_prev_pan_pos') or not self._prev_pan_pos:
            self._prev_pan_pos = pos
            self._prev_phi = self.Plot.phi
            return
        dx = pos[0] - self._prev_pan_pos[0]
        norm_width = self.width() / 2*np.pi
        self.Plot.phi = self._prev_phi - (dx/self.width()*np.pi*2)
        self.Plot.update()

    def mPosition(self):
        # See: 
        # Step 0: 2d viewport coordinates
        # these are x,y from -1 to 1
        mx = self._downpos.x()
        my = self._downpos.y()
        self.Candidates = []
        view_w = self.width()
        view_h = self.height()
        x = 2.0 * mx / view_w - 1.0
        y = 1.0 - (2.0 * my / view_h)
        PMi = self.projectionMatrix().inverted()[0]
        VMi = self.viewMatrix().inverted()[0]
        ray_clip = QtGui.QVector4D(x, y, -1.0, 1.0)
        ray_eye = PMi * ray_clip
        ray_eye.setW(0)
        # Convert to world coordinates
        ray_world = VMi * ray_eye
        ray_world = QtGui.QVector2D(ray_world.x(), ray_world.y(), ray_world.z())
        ray_world.normalize()
        origin = np.matrix(self.cameraPosition())
        ray_world = np.matrix([ray_world.x(), ray_world.y(), ray_world.z()])
        r = 125 # millimeters
        # Restore old colors before selecting new colors
        for idx,clr in zip(self.can_idx, self.can_clr):
            self.Plot.color[idx] = clr

        self.can_idx, self.can_clr = [], []
        for i, C in enumerate(self.Plot.pos):
            OC = origin - C
            b = np.inner(ray_world, OC)
            b = b.item(0)
            c = np.inner(OC, OC)
            c = c.item(0) - r**2
            bsqr = np.square(b)
            if (bsqr - c) >= 0:
                self.Candidates.append(self.Plot.pos[i])
                self.can_idx.append(i)
        for i in self.can_idx:
            self.can_clr.append(copy.deepcopy(self.Plot.color[i]))
            self.Plot.color[i] = np.array([1, 1, 1, 1])
        self.Plot.update()

class SibylWatchmanFlat(gl.GLScatterPlotItem):
    ''' Special version of the GLScatterPlotItem for flat map '''
    position = np.array([0,0,0,0])
    realcolors = np.array([0])
    weights = np.array([0])
    phi = 0 #Angle to rotate about 0 to 2pi or something
    boundaryRadius = 1.0
    boundaryHalfz = 1.0
    def __init__(self):
        super(SibylWatchmanFlat,self).__init__(pxMode=False, size=300)

    def update(self):
        super(SibylWatchmanFlat,self).update()
        self.project()

    def setPosition(self, position):
        self.position = position

    def setColor(self, colors):
        self.realcolors = colors

    def setWeights(self, weights):
        self.weights = weights

    def setBoundaries(self, radius, halfz):
        self.boundaryRadius = radius
        self.boundaryHalfz = halfz

    def project(self):
        # Cylindrical coordinates
        posX = (self.position.T)[0]
        posY = (self.position.T)[1]
        posZ = (self.position.T)[2]
        posRho = (posX**2 + posY**2)**0.5
        angle = np.arctan2(posX, posY) + self.phi
        posX = posRho * np.cos(angle)
        posY = posRho * np.sin(angle)
        angle = np.arctan2(posX, posY)
        #angle = angle%(2*np.pi)-np.pi
        theta = angle * posRho
        ## Move this to a detector class, but for now, arbitrarily define
        ## various regions.
        ## displace top and bottom
        zCut = self.boundaryHalfz-100 #6300
        rCut = self.boundaryRadius #6400
        vetoH = zCut + 200 #6500
        flatSelect = (posZ<zCut)&(posZ>-zCut)&(posRho<rCut)
        flatPOS = np.array([theta[flatSelect], posZ[flatSelect], np.zeros(len(theta[flatSelect]))]).T
        flatCLR = self.realcolors[flatSelect]
        #brushes = self.realcolors[flatSelect].tolist()
        #brushes = [tuple(br) for br in brushes]
        #print(brushes[0])
        # Flat region
        #self.plotFlatMap.plot(theta[flatSelect], posZ[flatSelect], pen=None, 
        #        symbol='o')
                #symbolBrushes=[pg.mkBrush(br) for br in brushes] )
        # Top Circle
        topSelect = (posZ>=zCut)&(posZ<vetoH)
        topShift = zCut*2.1 #1.5
        topX, topY = posX[topSelect], -posY[topSelect] + topShift
        topPOS = np.array([topX, topY, np.zeros(len(topX))]).T
        topCLR = self.realcolors[topSelect]
        #self.plotFlatMap.setData(x=topX, y=topY, z=self.CHARGE[topSelect])
        #self.plotFlatMap.update()
        #self.plotFlatMap.plot(topX, topY, pen=None, symbol='o')
        # Bot Circle
        botSelect = (posZ<=-zCut)&(posZ>-vetoH)
        botShift = -zCut*2.1 #1.5
        botX, botY = posX[botSelect], posY[botSelect] + botShift
        botPOS = np.array([botX, botY, np.zeros(len(botX))]).T
        botCLR = self.realcolors[botSelect]
        #self.plotFlatMap.plot(botX, botY, pen=None, symbol='o')

        ## Lets add VETO pmts in a scaled version.
        scaleFactor = 0.25 # This should fit, lets go upper left
        v_flatSelect = (posZ<zCut)&(posZ>-zCut)&(posRho>rCut)
        xShift = -2.1*zCut
        yShift = 2*zCut
        v_flatPOS = np.array(
                [theta[v_flatSelect]*scaleFactor + xShift, 
                posZ[v_flatSelect]*scaleFactor + yShift, 
                np.zeros(len(theta[v_flatSelect]))]).T
        v_flatCLR = self.realcolors[v_flatSelect]
        # Top Veto
        v_topSelect = (posZ>=vetoH)
        v_topShift = zCut*2.1 * scaleFactor + yShift
        v_topX = posX[v_topSelect]*scaleFactor+xShift
        v_topY = -posY[v_topSelect]*scaleFactor + v_topShift
        v_topPOS = np.array([v_topX, v_topY, np.zeros(len(v_topX))]).T
        v_topCLR = self.realcolors[v_topSelect]
        # Bot Veto
        v_botSelect = (posZ<=-vetoH)
        v_botShift = -zCut*2.1 * scaleFactor + yShift
        v_botX = posX[v_botSelect]*scaleFactor+xShift
        v_botY = posY[v_botSelect]*scaleFactor + v_botShift
        v_botPOS = np.array([v_botX, v_botY, np.zeros(len(v_botX))]).T
        v_botCLR = self.realcolors[v_botSelect]


        totalPOS = np.concatenate([topPOS, botPOS, flatPOS, 
                                   v_flatPOS, v_topPOS, v_botPOS])
        totalCLR = np.concatenate([topCLR, botCLR, flatCLR, 
                                   v_flatCLR, v_topCLR, v_botCLR])

        self.pos = totalPOS
        self.color = totalCLR
