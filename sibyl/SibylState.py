import pickle

class SibylState:
    '''
    Save and load the sibyl parameter state. Control default
    config files and exclude particular parameters
    '''
    def __init__(self, state=None):
        if state is None:
            self.state = {}
        else:
            self.state = state
        self.createBlacklist()

    def createBlacklist(self):
        self.bl = [
                'trackingEnabled', # This could be removed from this list
                'filename',
                'entries',
                'charge',
                'time',
                'posArray',
                'plWeights',
                'trackPosition',
                'trackColors',
                ]

    def saveState(self, filename='.default'):
        cleanState = {k:v for k,v in self.state.items() if k not in self.bl}
        with open(filename, 'wb') as outfile:
            pickle.dump(cleanState, outfile)

    def loadState(self, realState, filename='.default'):
        try:
            with open(filename, 'rb') as infile:
                copyState = pickle.load(infile)
            realState.update(copyState)
            self.state = realState
        except FileNotFoundError:
            pass


