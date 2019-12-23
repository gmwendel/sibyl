import ctypes as ct
import os
import numpy as np
from glob import glob

c_int = ct.c_int
c_double = ct.c_double
c_void_p = ct.c_void_p
c_char_p = ct.c_char_p
c_double_pointer = ct.POINTER(c_double)
c_int_pointer = ct.POINTER(c_int)
libname = glob(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "fastrat.*.so")
)[0]


class snake:
    def __init__(self):
        self.lib = ct.cdll.LoadLibrary(libname)
        self._init_functions()

    def _init_functions(self):
        self.l_square = self.lib.square
        self.l_square.restype = c_double_pointer
        self.l_square.argstypes = [c_double_pointer, c_int]

        self._freeSquare = self.lib.freeSquare
        self._freeSquare.restype = c_void_p
        self._freeSquare.argstypes = []

        self._openFile = self.lib.openFile
        self._openFile.restype = c_void_p
        self._openFile.argstypes = [c_char_p]

        self._getEntries = self.lib.getEntries
        self._getEntries.restype = c_int
        self._getEntries.argstypes = []

        self._getBoundaryRadius = self.lib.getBoundaryRadius
        self._getBoundaryRadius.restype = c_double
        self._getBoundaryRadius.argstypes = []

        self._getBoundaryHalfz = self.lib.getBoundaryHalfz
        self._getBoundaryHalfz.restype = c_double
        self._getBoundaryHalfz.argstypes = []

        self._getEvent = self.lib.getEvent
        self._getEvent.restype = c_void_p
        self._getEvent.argstypes = [c_int]

        self._getNHit = self.lib.getNHit
        self._getNHit.restype = c_int
        self._getNHit.argstypes = []
        self._getPMTCount = self.lib.getPMTCount
        self._getPMTCount.restype = c_int
        self._getPMTCount.argstypes = []
        self._getPMTX = self.lib.getPMTX
        self._getPMTX.restype = c_double_pointer
        self._getPMTX.argstypes = []
        self._getPMTY = self.lib.getPMTY
        self._getPMTY.restype = c_double_pointer
        self._getPMTY.argstypes = []
        self._getPMTZ = self.lib.getPMTZ
        self._getPMTZ.restype = c_double_pointer
        self._getPMTZ.argstypes = []
        self._getCharge = self.lib.getCharge
        self._getCharge.restype = c_double_pointer
        self._getCharge.argstypes = []
        self._getTime = self.lib.getTime
        self._getTime.restype = c_double_pointer
        self._getTime.argstypes = []

        self._getTracking = self.lib.getTracking
        self._getTracking.restype = c_void_p
        self._getTracking.argstypes = []
        self._getTrackCount = self.lib.getTrackCount
        self._getTrackCount.restype = c_int
        self._getTrackCount.argstypes = []
        self._getTrackX = self.lib.getTrackX
        self._getTrackX.restype = c_double_pointer
        self._getTrackX.argstypes = []
        self._getTrackY = self.lib.getTrackY
        self._getTrackY.restype = c_double_pointer
        self._getTrackY.argstypes = []
        self._getTrackZ = self.lib.getTrackZ
        self._getTrackZ.restype = c_double_pointer
        self._getTrackZ.argstypes = []
        self._getTrackNames = self.lib.getTrackNames
        self._getTrackNames.restype = c_int_pointer
        self._getTrackNames.argstypes = []

    def openFile(self, fname):
        print(fname)
        self._openFile(fname.encode())

    def getEvent(self, idx):
        self._getEvent(idx)

    def getEntries(self):
        return self._getEntries()

    def getBoundaryRadius(self):
        return self._getBoundaryRadius()

    def getBoundaryHalfz(self):
        return self._getBoundaryHalfz()

    def getXYZ(self):
        pmtc = self._getPMTCount()
        x_ret = self._getPMTX()
        x_arr = np.fromiter(x_ret, dtype=np.float64, count=pmtc)
        y_ret = self._getPMTY()
        y_arr = np.fromiter(y_ret, dtype=np.float64, count=pmtc)
        z_ret = self._getPMTZ()
        z_arr = np.fromiter(z_ret, dtype=np.float64, count=pmtc)
        return x_arr, y_arr, z_arr

    def getHitInfo(self):
        pmtc = self._getPMTCount()
        c_ret = self._getCharge()
        c_arr = np.fromiter(c_ret, dtype=np.float64, count=pmtc)
        t_ret = self._getTime()
        t_arr = np.fromiter(t_ret, dtype=np.float64, count=pmtc)
        return c_arr, t_arr

    def square(self, arr):
        arr = np.array(arr, dtype=np.float64)
        asend = arr.ctypes.data_as(c_double_pointer)
        areturn = self.l_square(asend, len(arr))
        newarr = np.fromiter(areturn, dtype=np.float64, count=len(arr))
        self._freeSquare()
        return newarr

    ## Tracking information
    def getTracking(self):
        self._getTracking()

    def getTrackSteps(self):
        tsteps = self._getTrackCount()
        x_ret = self._getTrackX()
        x_arr = np.fromiter(x_ret, dtype=np.float64, count=tsteps)
        y_ret = self._getTrackY()
        y_arr = np.fromiter(y_ret, dtype=np.float64, count=tsteps)
        z_ret = self._getTrackZ()
        z_arr = np.fromiter(z_ret, dtype=np.float64, count=tsteps)
        n_ret = self._getTrackNames()
        n_arr = np.fromiter(n_ret, dtype=np.int64, count=tsteps)
        return x_arr, y_arr, z_arr, n_arr


if __name__ == "__main__":
    a = np.random.rand(10)
    cobra = snake()
    cobra.square(a)
    cobra.openFile("cobalt_watchman.root")
    cobra.getEvent(0)
