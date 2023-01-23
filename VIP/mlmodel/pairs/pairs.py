'''
Description goes here
'''


class Pairs:

    def __init__(self, virus, host) -> None:
        
        self.virus = virus
        self.host = host

        # features
        self.GCdiff = None
        self.k3dist = None
        self.k6dist = None
        
        self.spacers = False
        self.blastn = False
    
    

