import numpy as np
import pyqtgraph as pg
import pyqtgraph.opengl as gl
from PyQt5 import QtGui, QtCore

class Sibyl3DViewer(gl.GLViewWidget):
    App = None
    can_idx = []
    can_clr = []
    def __init__(self, plot, paths, app=None):
        if self.App is None:
            if app is not None:
                self.App = app
            else:
                self.App = QtGui.QApplication([])
        super(Sibyl3DViewer,self).__init__()
        self.paths = paths
        self.Poss = []
        self.Plot = plot
        self.addItem(self.Plot)
        self.addItem(self.paths)
        self._downpos = []
        self._init_camera()
        self.update()

    # Fix for "high dpi":
    def width(self):
        trueWidth = super(Sibyl3DViewer,self).width()
        pixRatio = super(Sibyl3DViewer,self).devicePixelRatio()
        return trueWidth * pixRatio

    def height(self):
        trueHeight = super(Sibyl3DViewer,self).height()
        pixRatio = super(Sibyl3DViewer,self).devicePixelRatio()
        return trueHeight * pixRatio

    def _init_camera(self):
        # init camera
        self.opts['center'] = QtGui.QVector3D(0, 0, 0)
        self.setCameraPosition(pos=np.array([0, 0, 0]))
        self.opts['elevation'] = 25 
        self.opts['distance'] = 30000
        self.opts['azimuth'] = 0 

    # All keyboard and mouse event handlers
    def keyPressEvent(self, ev):
        super(Sibyl3DViewer,self).keyPressEvent(ev)
        if ev.key() == ord(' '):
            self._init_camera()
        if ev.key() == ord('O'):
            self.orbit(45, 45)

    def mousePressEvent(self, ev):
        ''' Store the position of the mouse press for later use '''
        super(Sibyl3DViewer, self).mousePressEvent(ev)
        self._downpos = self.mousePos

    def mouseReleaseEvent(self, ev):
        ''' Allow for single click to move and right click for context menu '''
        super(Sibyl3DViewer, self).mouseReleaseEvent(ev)
        if self._downpos == ev.pos():
            x = ev.pos().x()
            y = ev.pos().y()
            if ev.button() == 2:
                self.mPosition()
            elif ev.button() == 1:
                x = x - self.width() / 2
                y = y-self.height() / 2
        self._prev_zoom_pos = None
        self._prev_pan_pos = None

    def mouseMoveEvent(self, ev):
        ''' Allow shift to move and ctrl to pan '''
        shift = ev.modifiers() & QtCore.Qt.ShiftModifier
        ctrl = ev.modifiers() & QtCore.Qt.ControlModifier
        if shift:
            y = ev.pos().y()
            if not hasattr(self, '_prev_zoom_pos') or not self._prev_zoom_pos:
                self._prev_zoom_pos = y
                return
            dy = y - self._prev_zoom_pos
            def delta():
                return -dy * 5
            ev.delta = delta
            self._prev_zoom_pos = y
            self.wheelEvent(ev)
        elif ctrl:
            pos = ev.pos().x(), ev.pos().y()
            if not hasattr(self, '_prev_pan_pos') or not self._prev_pan_pos:
                self._prev_pan_pos = pos
                return
            dx = pos[0] - self._prev_pan_pos[0]
            dy = pos[1] - self._prev_pan_pos[1]
            self.pan(dx, dy, 0, relative=True)
            self._prev_pan_pos = pos
        else:
            super(Sibyl3DViewer, self).mouseMoveEvent(ev)

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
        ray_world = QtGui.QVector3D(ray_world.x(), ray_world.y(), ray_world.z())
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
