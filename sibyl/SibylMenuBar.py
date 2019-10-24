from PyQt5.QtWidgets import QMenuBar
#from PyQt5 import QtGui
from PyQt5.QtGui import *
from PyQt5.QtCore import *

class SibylMenuBar:
    '''
    Menu bar
    '''
    def __init__(self, parent):
        self.parent = parent
        self.menubar = self.parent.menuBar()
        self.setupMenu()

    def setupMenu(self):
        self.fileMenu()

    def fileMenu(self):
        fMenu = self.menubar.addMenu("File")
        # Open
        openAction = fMenu.addAction("Open")
        # Options
        optionsAction = fMenu.addAction("Options")
        optionsAction.triggered.connect(self.openOptions)
        # Close
        closeAction = fMenu.addAction("Close")
        closeAction.triggered.connect(self.parent.close)
        closeAction.setShortcut("CTRL+W")

    def openOptions(self):
        self.optionsMenu = SibylOptions(self.parent, self)

class SibylOptions(QDialog):
    App=None
    def __init__(self, app=None, parent=None):
        super(SibylOptions,self).__init__()
        if self.App is None:
            if app is not None:
                self.App = app
            else:
                self.App = QtGui.QApplication([])
        self.parent = parent
        self.setGeometry()
        self.show()

    def setGeometry(self):
        self.setMinimumSize(1200, 800)
        self.setMaximumSize(1200, 800)
        k = QApplication.desktop()
        w = k.screenGeometry().width()
        h = k.screenGeometry().height()
        self.move( (w-self.width())*0.1, (h-self.height())/2 )
        self.setWindowTitle("Options")

    def keyPressEvent(self, ev):
        if ( ev.key() == ord('W') ) and ( ev.modifiers() & Qt.ControlModifier ):
            self.close()
