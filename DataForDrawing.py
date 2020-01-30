# Udaje, ktore naplni a vracia DrawDistributions funkcia draw(), ktore potom pouzije funkcia drawErrors()

# data - pole objektov DistributionData
# num_of_shifts - pocet posunuti
# distributions - pole objektov DistributionTuple

class DataForDrawing:
    def __init__(self, data, num_of_shifts, distributions):
        self.data = data
        self.num_of_shifts = num_of_shifts
        self.distributions = distributions
