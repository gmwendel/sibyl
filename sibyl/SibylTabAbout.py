from PyQt5.QtGui import *
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
from .SibylTab import SibylTab
import os.path as path

class SibylTabAbout(SibylTab):
    '''
    Main event display tab.
    '''
    def __init__(self, parent=None):
        super(SibylTabAbout,self).__init__(parent)
        self.buildWidget()

    def buildWidget(self):
        filename = path.join( path.dirname( path.abspath(__file__) ), "assets", "readme.md" )
        with open(filename) as rr:
            text = rr.read()
        browser = QTextBrowser()
        try:
            import markdown as md
            text = "<head><link rel='stylesheet' type='text/css' href='style.css'></head>"+\
                    md.markdown(text)
        except ModuleNotFoundError as e:
            print(e)
            text = "<h1>INSTALL python-markdown to render this page properly</h1>" \
                    + "<p>" + text + "</p>"
        browser.setHtml(text)
        mainLayout = QVBoxLayout()
        mainLayout.addWidget(browser)
        self.setLayout(mainLayout)
