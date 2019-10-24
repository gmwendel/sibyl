from PyQt5 import QtGui, QtCore
import os.path as path

class SibylTrayIcon(QtGui.QSystemTrayIcon):
    def __init__(self, parent=None):
        iconfile = path.join( path.dirname( path.abspath(__file__) ), 
                'assets/logo.png' )
        icon = QtGui.QIcon(iconfile)
        QtGui.QSystemTrayIcon.__init__(self, icon, parent)
        self.menu = QtGui.QMenu(parent)
        exitAction = self.menu.addAction("Exit")
        self.setContextMenu(self.menu)
