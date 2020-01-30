# Udaje pre jednu distribuciu

# distribucia = scipy.stats.distribucia
# p0 -sa pouziva pri optimize.curve_fit
# name - na zobrazenie mena pri vykresleni grafu
# color - farba pouzita pri vykresleni

class DistributionTuple:
    def __init__(self, distribution, p0, name, color):
        self.distribution = distribution
        self.p0 = p0
        self.name = name
        self.color = color