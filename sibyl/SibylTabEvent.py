from PyQt5.QtGui import *
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
from .SibylTab import SibylTab
from .Sibyl2DViewer import Sibyl2DViewer, SibylWatchmanFlat
from .Sibyl3DViewer import Sibyl3DViewer
from .SibylHistogram import SibylHistogram
import pyqtgraph.opengl as gl
import os.path as path

class SibylTabEvent(SibylTab):
    '''
    Main event display tab.
    '''
    def __init__(self, parent=None):
        super(SibylTabEvent,self).__init__(parent)
        self.buildWidget()
        #self._parent.newEvent()
        #self.drawEvent()

    def buildWidget(self):
        formGroupBox = QGroupBox()
        formLayout = QFormLayout()
        self.layout = QGridLayout()
        self.hboxchooser = QHBoxLayout()
        formLayout.addRow(self.hboxchooser)

        ## Choose charge mode
        chargeTime_button = QPushButton("Q/T Toggle")
        chargeTime_button.clicked.connect(self.toggleColorMode)
        self.hboxchooser.addWidget(chargeTime_button)

        ## Fit current event
        fit_button = QPushButton()
        fit_button.clicked.connect(self.performFit)
        logoFile = path.join( path.dirname( path.abspath(__file__) ),
                'assets/bonsai.png')
        fit_button.setIcon(QIcon(logoFile))
        fit_button.setToolTip('Bonsai!')
        self.hboxchooser.addWidget(fit_button)

        ## Put form together
        formGroupBox.setLayout(formLayout)
        self.layout.addWidget(formGroupBox,5,12,4,4)
        self.figCanvas = SibylHistogram(self.App, self)
        self.layout.addWidget(self.figCanvas,5,8,4,4)

        # If we want to pop-out widgets instead of putting them
        # into a layout, simply use widget.show() i.e.
        # self.figCanvas.show()

        ## Now the 3D graph
        self.plot3DView = gl.GLScatterPlotItem(pxMode=False)
        self.plotTracks = gl.GLLinePlotItem(mode='lines', color=(1,0,0,0.1))
        self.plotFlatMap = SibylWatchmanFlat()
        self.glWin = Sibyl3DViewer(self.plot3DView, self.plotTracks, self.App)
        self.flatMapWindow = Sibyl2DViewer(self.plotFlatMap, self.App)
        self.layout.addWidget(self.flatMapWindow,0,0,9,8)
        self.layout.addWidget(self.glWin,0,8,5,8)

        # Fix the column and row stretch
        for r in range(9):
            self.layout.setRowStretch(r, 1)
        for c in range(16):
            self.layout.setColumnStretch(c, 1)

        self.setLayout(self.layout)

    def performFit(self):
        print('Bonsai!')

    def toggleColorMode(self):
        if self.parameters['colorMask'] == 'charge':
            self.parameters['colorMask'] = 'time'
        else:
            self.parameters['colorMask'] = 'charge'
        self.drawEvent()

    def drawEvent(self):
        self.colorize()
        ## 3D Viewport
        if self.plot3DView is None:
            self.plot3DView = gl.GLScatterPlotItem(pos=self.parameters["posArray"],
                    color=self.colorArray, size=self.parameters["plWeights"], pxMode=False)
        else:
            self.plot3DView.pos = self.parameters["posArray"]
            self.plot3DView.color = self.colorArray
            self.plot3DView.size = self.parameters["plWeights"]
            self.plot3DView.update()
        ## 2D Viewport
        self.plotFlatMap.setPosition(self.parameters["posArray"])
        self.plotFlatMap.setColor(self.colorArray)
        self.plotFlatMap.setWeights(self.parameters["plWeights"])
        self.plotFlatMap.update()

        if self.parameters['colorMask'] == 'charge':
            variable = self.parameters["charge"]
        else:
            variable = self.parameters["time"]
        self.figCanvas.setData(variable)

        ## Tracking
        if self.parameters["trackingEnabled"]:
            self.plotTracks.pos = self.parameters["trackPosition"]
            self.plotTracks.color = self.parameters["trackColors"]
            self.plotTracks.update()

    def drawColors(self):
        '''
        Only update colors, not positions or sizes
        '''
        self.colorize()
        self.plot3DView.color = self.colorArray
        self.plot3DView.update()
        self.plotFlatMap.setColor(self.colorArray)
        self.plotFlatMap.update()

    def colorize(self):
        # choose color mask, charge for now
        if self.parameters['colorMask'] == 'charge':
            variable = self.parameters["charge"]
        else:
            variable = self.parameters["time"]
        true_min = 0.0
        min_var = true_min
        max_var = self.parameters["histXMax"]
        if self.parameters["onlyHits"]:
            null_color = [1,1,1,0.0]
        else:
            null_color = [1,1,1,0.25]
        var_range = max_var - min_var
        varNorm = (variable-min_var)*1/var_range
        self.colorArray = self.parameters['colorMap'](varNorm)
        self.colorArray[variable<true_min] = null_color
        ## Invisible hack
        if self.parameters["invisible"]:
            self.colorArray *= 0
