class Nexus():
    def __init__(self, nexusID = None, order = -1):
        self.nexusID = nexusID
        self.order = order

class DendriticReachHeadNexus(Nexus): # as DRN
    '''Extend this to have 0-upstreams and 1-downstream reaches'''
    def __init__(upstreams = None, downstream = None):
        self.upstreams = upstreams
        self.downstream = downstream

class DendriticReachJuntionNexus(Nexus): # as DRN
    '''Extend this to have N-upstreams and 1-downstream reaches'''
    def __init__():
        pass
    # dict upstreams Reaches
    # dict downstream Reach
    def passFlow(isThisThePredictorOrCorrectorStep = 'Predictor'):
        '''
        method to handle flow passing
        with a flag to distinguish behavior in the
        predictor vs corrector steps'''
        pass

class DendriticReachTerminalNexus(Nexus): # as DRN
    '''Extend this to have N-upstreams and 0-downstream reaches'''
    def __init__():
        pass
    # dict upstreams Reaches
    # dict downstream Reach

