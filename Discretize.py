# Trieda na diskretizaciu spojitych distribucii

class Discretize:
    def __init__(self, cdf):
        self.cdf = cdf
    
    def pmf(self, x, *params):
        return (self.cdf(x+0.5, *params) - self.cdf(x-0.5, *params))/(1-self.cdf(0.5, *params))