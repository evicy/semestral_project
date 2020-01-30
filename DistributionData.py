# Trieda pre udaje pri vykreslovani

# distribution_num - ocislovanie distribucie
# shift - pocet shiftnuti
# popt - Optimal values for the parameters 
# li_error - error pri danej norme

class DistributionData:
    def __init__(self, distribution_num, shift, popt, l1_error, l2_error):
        self.distribution_num = distribution_num
        self.shift = shift
        self.popt = popt
        self.l1_error = l1_error
        self.l2_error = l2_error