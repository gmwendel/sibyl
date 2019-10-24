from PyQt5.QtGui import *
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *

class SibylTab(QWidget):
    '''
    ABS for tabs added to the main window.
    '''
    App = None
    def __init__(self, parent=None):
        super(SibylTab,self).__init__()
        self._parent = parent
        if self.App is None:
            self.App = parent.App
        self.parameters = self._parent.parameters

    def drawEvent(self):
        pass
