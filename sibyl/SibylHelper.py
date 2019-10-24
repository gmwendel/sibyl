'''
SyblyHelper is a collection of stand-alone functions that
are not part of a particular class, that can be used everywhere
'''

import matplotlib.cm as cm
from PyQt5.QtGui import *
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *

class SibylColorChoice(QtWidget):
    def __init__(self, parent=None):
        self._parent=parent
        self.buildWindow()


