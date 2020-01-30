#Trieda Tailer nieco ako shifter

import numpy as np

class Tailer:
    #konstruktor, dostane distribuciu a numpy array normalizovane y hodnoty 
    def __init__(self, distr, y, q):
        self.distr = distr
        self.y = y
        self.q = q
        
    # vrati novy y, kde prvych q je povodnych a ostatne su normalizovane a prerobene podla distr
    # *params nam umozni mat viac parametrov
    def pmf(self, x, *params):
        y_temp = self.distr.pmf(x, *params)
        y_result = y_temp * (1-sum(self.y[0:self.q]))/(1-sum(y_temp[0:self.q]))
        for i in range(self.q):
            y_result[i] = self.y[i]
        return np.array(y_result)       

