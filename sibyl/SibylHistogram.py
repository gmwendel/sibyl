from PyQt5 import QtGui
from PyQt5.QtGui import *
from PyQt5.QtCore import *
import pyqtgraph as pg
import pyqtgraph.opengl as gl
import numpy as np
import matplotlib.cm as cm

class SibylHistogram(gl.GLViewWidget):
    '''
    SibylHistogram is a GLViewWidget (3D viewport), which uses the GLMeshItem
    to draw the vertices and faces of the histogram in different colors.
    '''
    App = None
    data  = 0   # Input data to be histogrammed
    xmin  = 0
    xmax  = 5.5
    nbins = 20 # Default binning
    txt   = '( )'
    zoom_y = 1
    zoom_x = 1
    autoMode = True
    logy = False
    def __init__(self, app=None, parent=None):
        if self.App is None:
            if app is not None:
                self.App = app
            else:
                self.App = QtGui.QApplication([])
        super(SibylHistogram,self).__init__()
        self._parent=parent
        self._init_camera()
        self._axes()
        self._menu()
        self.setMouseTracking(True)
        #self.setAutoFillBackground(False)
        self.histMesh = gl.GLMeshItem(smooth=False)
        self.addItem(self.histMesh)
        self.update()

    # Fix for "high dpi":
    def width(self, fixed=False):
        trueWidth = super(SibylHistogram,self).width()
        if fixed:
            return trueWidth
        pixRatio = super(SibylHistogram,self).devicePixelRatio()
        return trueWidth * pixRatio

    def height(self, fixed=False):
        trueHeight = super(SibylHistogram,self).height()
        if fixed:
            return trueHeight
        pixRatio = super(SibylHistogram,self).devicePixelRatio()
        return trueHeight * pixRatio

    def _axes(self):
        self.ax = gl.GLGridItem()
        self.ax.setSize(5.5,10,1)
        self.ax.setSpacing(5.5/5, 1, 1)
        self.ax.translate(5.5/2,10/2.,0)
        #self.ax.setSize(x=10)
        #self.ax = gl.GLBoxItem(glOptions='opaque')
        #self.ax.setSize(5.5, 0, 1)
        self.addItem(self.ax)

    def _menu(self):
        self.menu = QMenu()
        # Reset (like auto)
        resetAction = self.menu.addAction("Reset")
        resetAction.triggered.connect(self._autoset)
        # Set auto-mode
        autoModeAction = self.menu.addAction("Auto-coordinates")
        autoModeAction.setCheckable(True)
        if self.autoMode:
            autoModeAction.setChecked(True)
        autoModeAction.triggered.connect(self._toggle_auto)
        # Log-scale-y
        logyAction = self.menu.addAction("Log-y")
        logyAction.setCheckable(True)
        if self.logy:
            logyAction.setChecked(True)
        logyAction.triggered.connect(self._toggle_logy)
        # Color picker
        colorAction = self.menu.addAction("Colormap")
        colorAction.triggered.connect(self._selectColor)

    def _toggle_auto(self):
        self.autoMode = not self.autoMode

    def _toggle_logy(self):
        self.logy = not self.logy
        self.resetHist()
        self.resetMesh()

    def _autoset(self):
        # Set best x and y (x to edge, y at 90%)
        x, y = self.screenToWorld( (0, self.height()*0.1) )
        end_x = np.max(self.data)
        self.zoom_x = self.xmax/end_x if end_x > 0 else 1
        self.zoom_x = np.max([0, self.zoom_x])
        self.resetHist()
        self.zoom_y = y/np.max(self.hy)
        self.resetMesh()

    def _selectColor(self):
        self.selector = SibylColorSelector(self.App, self)

    def update(self):
        super(SibylHistogram,self).update()

    def setData(self, data):
        self.data = data
        if self.autoMode:
            self._autoset()
        else:
            self.resetHist()
            self.resetMesh()

    def resetHist(self):
        bins = np.linspace(self.xmin, self.xmax/self.zoom_x,self.nbins)
        hy, hx = np.histogram(
                    #self.data,
                    np.clip(self.data, -999, bins[-1]), 
                    bins=bins,
                    density=False)
        self.hx = hx
        self.hy = hy
        if self.logy:
            self.hy = np.log10(self.hy)
            self.hy[self.hy==-np.inf] = 0
            #self.hy += np.min(self.hy)

    def resetMesh(self):
        vtx, colrs, faces = self._verts(self.hx*self.zoom_x, self.hy*self.zoom_y)
        self.histMesh.setMeshData(vertexes=vtx, faces=faces, faceColors=colrs)
        self.histMesh.update()
        self._parent.parameters['histXMax'] = self.xmax/self.zoom_x
        if self._parent is not None:
            self._parent.drawColors()

    def _verts(self, x, y):
        # points, N bins -> 3*(N-2)+7
        # for x, [0,0,1,1,1, ..., N, N]
        nbins = len(y)
        x_verts     = np.repeat(x, 3)[1:-1]
        temp_y      = np.repeat(y, 3)
        temp_y[::3] = 0
        y_verts     = np.concatenate([temp_y, [0]])
        z_verts     = np.zeros(x_verts.shape)
        pos = np.array([x_verts, y_verts, z_verts]).T
        # Color stuff
        #clrs = np.repeat(cm.jet(x[:-1]), 2, axis=0)
        #clrs = cm.cividis(np.linspace(0,1,nbins))
        clrs = self._parent.parameters['colorMap'](np.linspace(0, 1, nbins))
        # Faces
        k,j = np.arange(nbins*3), np.arange(nbins*3)
        face1 = k.reshape(nbins,3)
        k[1::3] += 2
        face2 = j.reshape(nbins,3)
    
        face = np.concatenate([face1, face2])
        clrs = np.concatenate([clrs, clrs])
        return pos, clrs, face

    def _init_camera(self):
        # init camera -- this is pseudo 2D
        self.opts['center'] = QVector3D(0, 0, 0)
        self.setCameraPosition(pos=np.array([0, 0, 0]))
        self.opts['elevation'] = 90
        self.opts['distance']  = 350
        self.opts['fov']       = 1
        self.opts['azimuth']   = -90

    def paintGL(self, *args, **kwargs):
        #self.paintEvent(*args, **kwargs)
        super(SibylHistogram,self).paintGL(*args, **kwargs)
        self.drawLegend()
    #    pass

    def paintEvent(self, event, *args, **kwargs):
        #super(SibylHistogram,self).paintGL(*args, **kwargs)
        super(SibylHistogram,self).paintEvent(event, *args, **kwargs)
        #self.drawLegend()


    def drawLegend(self):
        margin = 11
        padding = 6
        pixRatio = super(SibylHistogram,self).devicePixelRatio() #needed to scale font size and wstart
        self.qglColor(Qt.white)
        cmask = (self._parent.parameters['colorMask']).capitalize()
        fsize = int(self.width()/(14.0*pixRatio))
        font = QFont("Times", fsize, QFont.Bold)
        wstart = self.width() - fsize*8*pixRatio
        self.renderText(wstart/pixRatio, 0.1*self.height(fixed=True),
                cmask, font=font)
        self.renderText(wstart/pixRatio, 0.2*self.height(fixed=True), 
                self.txt, font=font)

    '''
    Interactive events: mouse, keys, wheel, resize
    '''

    def contextMenuEvent(self, event):
        self.menu.exec_(self.mapToGlobal(event.pos()))

    def resizeEvent(self, event):
        self.plotToBottom()
        if self.autoMode:
            self._autoset()

    def wheelEvent(self, event):
        delta = event.angleDelta().x()
        if delta == 0:
            delta = event.angleDelta().y()
        # Change binning with shift-scroll
        if event.modifiers() & Qt.ShiftModifier:
            self.nbins += delta/20.0
            self.nbins = int(np.max([4, self.nbins]))
            self.resetHist()
        # Pan x-axis with fixed nbins (change xmax)
        elif event.modifiers() & Qt.ControlModifier:
            self.zoom_x *= 1.001**delta
            self.resetHist()
        # Zoom along y-axis
        else:
            self.zoom_y *= 1.001**delta
        self.resetMesh()

    def mouseDoubleClickEvent(self, event):
        if event.button() == 1:
            self._autoset()

    def mousePressEvent(self, event):
        super(SibylHistogram,self).mousePressEvent(event)
        if event.button() == 1:
            self._mouseDown = True

    def mouseReleaseEvent(self, event):
        super(SibylHistogram,self).mouseReleaseEvent(event)
        if event.button() == 1:
            self._mouseDown = False
            try:
                del(self._prev_pan_pos)
            except AttributeError:
                pass

    def mouseMoveEvent(self, event):
        mousePos = event.pos()
        x, y = mousePos.x(), mousePos.y()
        mx, my = self.screenToWorld( (x, y) )
        my /= self.zoom_y
        mx /= self.zoom_x
        self.txt = '(%0.2f, %0.2f)'%(mx, my)
        self.update()
        if hasattr(self, '_mouseDown') and self._mouseDown:
            if not hasattr(self, '_prev_pan_pos') or not self._prev_pan_pos:
                self._prev_pan_pos = x, y
                return
            dx = x - self._prev_pan_pos[0]
            dy = y - self._prev_pan_pos[1]
            self._prev_pan_pos = x, y
            self.zoom_x *= 1.003**dx
            self.zoom_y *= 0.997**dy
            self.resetHist()
            self.resetMesh()

    def leaveEvent(self, event):
        super(SibylHistogram,self).leaveEvent(event)
        self.txt = '()'
        self.update()

    def plotToBottom(self):
        self.opts['center'] = QVector3D(0, 0, 0)
        bx, by = self.screenToWorld( (0, self.height() ) )
        self.opts['center'] = QVector3D(-bx*0.9, -by*0.90, 0)

    def worldToScreen(self, world):
        # world is a 4D vector, (x,y,z,w)
        pass

    def screenToWorld(self, screen):
        mx, my = screen
        view_w = self.width()
        view_h = self.height()
        x = 2.0 * mx / view_w - 1.0
        y = 1.0 - (2.0 * my / view_h)
        PMi = self.projectionMatrix().inverted()[0]
        VMi = self.viewMatrix().inverted()[0]
        ray_clip = QVector4D(x, y, -1.0, 1.0)
        ray_eye = PMi * ray_clip
        ray_eye.setW(0)
        # Convert to world coordinates
        ray_world = VMi * ray_eye
        ray_world = QVector3D(ray_world.x(), ray_world.y(), ray_world.z())
        ray_world.normalize()
        origin = np.matrix(self.cameraPosition())
        origin = self.cameraPosition()
        # intercept with XY-plane
        normal = QVector3D(0, 0, 1)
        pzero = QVector3D(0, 0, 0)
        d = self.dot( (pzero - origin), normal ) / self.dot(ray_world, normal)
        intercept = d*ray_world + origin
        return intercept.x(), intercept.y()

    def dot(self, v1, v2):
        return QVector3D.dotProduct(v1, v2)

class SibylColorSelector(QScrollArea):
    App = None
    cmaps = [('Perceptually Uniform Sequential', [
                'viridis', 'plasma', 'inferno', 'magma', 'cividis']),
             ('Sequential', [
                'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
                'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
                'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn']),
             ('Sequential (2)', [
                'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone', 'pink',
                'spring', 'summer', 'autumn', 'winter', 'cool', 'Wistia',
                'hot', 'afmhot', 'gist_heat', 'copper']),
             ('Diverging', [
                'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu',
                'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic']),
             ('Cyclic', ['twilight', 'twilight_shifted', 'hsv']),
             ('Qualitative', [
                'Pastel1', 'Pastel2', 'Paired', 'Accent',
                'Dark2', 'Set1', 'Set2', 'Set3',
                'tab10', 'tab20', 'tab20b', 'tab20c']),
             ('Miscellaneous', [
                'flag', 'prism', 'ocean', 'gist_earth', 'terrain', 'gist_stern',
                'gnuplot', 'gnuplot2', 'CMRmap', 'cubehelix', 'brg',
                'gist_rainbow', 'rainbow', 'jet', 'nipy_spectral', 'gist_ncar'])]
    def __init__(self, app=None, parent=None):
        super(SibylColorSelector,self).__init__()
        if self.App is None:
            if app is not None:
                self.App = app
            else:
                self.App = QtGui.QApplication([])
        self.parent = parent
        self.setGeometry()
        self.addButtons()
        self.show()

    def setGeometry(self):
        self.setMinimumSize(300, 800)
        self.setMaximumSize(300, 800)
        k = QApplication.desktop()
        w = k.screenGeometry().width()
        h = k.screenGeometry().height()
        self.move( (w-self.width())*(0.85), (h-self.height())/2)
        self.setWindowTitle("Color Selector")

    def addButtons(self):
        scrollWidget = QWidget()
        scrollLayout = QFormLayout()
        scrollWidget.setLayout(scrollLayout)
        self.setWidget(scrollWidget)
        self.setWidgetResizable(True)
        # Header
        title = QLabel("Color Palette Selection")
        title.setAlignment(Qt.AlignCenter)
        title.setFont(QFont("Arial", 14, QFont.Bold))
        scrollLayout.addRow(title)

        # Button
        for name,maplist in self.cmaps:
            try:
                tlabel = QLabel(name)
                tlabel.setAlignment(Qt.AlignCenter)
                tlabel.setFont(QFont("Arial", 12))
                scrollLayout.addRow(tlabel)
                for c in maplist:
                    mp = cm.get_cmap(c)
                    clr_button = QPushButton()
                    clr_button.setToolTip(c)
                    pos = np.linspace(0,1,100)
                    clr = mp(pos, bytes=True)
                    color = [QColor(*(x)).name() for x in clr]
                    stops = [',stop:%0.2f %s'%(x,y) for x,y in zip(pos,color)]
                    styleString = 'QPushButton {background-color:' + \
                            'qlineargradient(x1:0,y1:0,x2:1,y2:0' + \
                            '%s);}' % ''.join(stops)
                    clr_button.setStyleSheet(styleString)
                    clr_button.setFixedHeight(20)
                    clr_button.clicked.connect(self.setColor(mp) )
                    #clr_button.setPalette(pal)
                    #clr_button.update()
                    #clr_button.setAutoFillBackground(True)
                    ## Test gradient
                    scrollLayout.addRow(clr_button)
            except ValueError:
                pass

    def setColor(self, color):
        def sCol():
            self.parent._parent.parameters['colorMap'] = color
            self.parent.resetMesh()
        return sCol

    def keyPressEvent(self, ev):
        if ( ev.key() == ord('W') ) and ( ev.modifiers() & Qt.ControlModifier ):
            self.close()
