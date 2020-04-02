import scipy.stats as st
from numpy import linalg as LA
import Tailer as tailer
import Discretize as disc
import numpy as np
import DistributionData as dist_data
import matplotlib.pyplot as plt
from scipy import optimize
import DataForDrawing as DataForDrawing
import DistributionTuple as dTuple


def draw(histo, max_x = 30, num_of_shifts = 4):
    # Funkcia vykresli vsetky distribucie * num_of_shifts pre data z histo (jellyfish histo file)
    # vypise 3 najlepsie distribucie s najmensim errorom (pre L1 a L2 normy)
    
    x = np.array(sorted(histo.keys()))
    y_raw = [histo[i] for i in x]
    y_sum = sum(y_raw)
    y = np.array([h/y_sum for h in y_raw])   
    # y obsahuje normalizovane hodnoty

    # vytvorime si DistributionTuples, ktore si ulozime do pola
    distributions = []
    
    distributions.append(dTuple.DistributionTuple(st.poisson, [50], 'Poisson', 'red')) #[2]
    distributions.append(dTuple.DistributionTuple(st.binom, [2000, 0.01], 'Binomial', 'gold')) #[1000, 0.01]
    distributions.append(dTuple.DistributionTuple(st.nbinom, [100, 0.7], 'Negative binomial', 'fuchsia')) #[2000, 0.01]
    #distributions.append(dTuple.DistributionTuple(st.geom, [0.0001], 'Geometric', 'gold'))
    #distributions.append(dTuple.DistributionTuple(disc.Discretize(st.pareto.cdf), [0.0001], 'Pareto', 'forestgreen'))
    #distributions.append(dTuple.DistributionTuple(disc.Discretize(st.gamma.cdf), [0.5], 'Gamma', 'royalblue'))
    #distributions.append(dTuple.DistributionTuple(disc.Discretize(st.cauchy.cdf), [0.0001], 'Cauchy', 'fuchsia'))
    distributions.append(dTuple.DistributionTuple(disc.Discretize(st.norm.cdf), [20, 1], 'Gauss', 'black'))
    #distributions.append(dTuple.DistributionTuple(st.betabinom, [len(x), 1, 1], 'BetaBinom', 'black'))

    #pole objektov DistributionData
    data = []

    # vykreslovanie vsetky distribucie * num_of_shifts
    for i in range(num_of_shifts):
        plt.figure(figsize=(10,5))
        plt.bar(x, y, color = 'lightblue')
        plt.yscale("log")
        
        for j in range(len(distributions)):
            tail = tailer.Tailer(distributions[j].distribution, y, i)
            popt, pcov = optimize.curve_fit(tail.pmf, xdata=x, ydata=y, p0=distributions[j].p0)
            xspace = np.linspace(1, x[len(x)-1], num=x[len(x)-1])
            if(j % 2 == 1):
                plt.plot(xspace, tail.pmf(xspace,*popt), linestyle='dashed', color=distributions[j].color, linewidth=2, label=distributions[j].name)
            else:
                plt.plot(xspace, tail.pmf(xspace,*popt), color=distributions[j].color, linewidth=2, label=distributions[j].name)
            a = tail.pmf(xspace,*popt)
            l2_error = LA.norm(a[:len(y)] - y[:len(y)])
            l1_error = LA.norm(a[:len(y)] - y[:len(y)], ord=1)
            data.append(dist_data.DistributionData(j, i, popt,  l1_error, l2_error)) 
        
        plt.ylim(min(y), 2)
        plt.xlim(0, max_x)
        plt.title('Histogram o početnosti čítaní', fontsize=20)
        plt.xlabel('Čítania na intervaloch', fontsize=15)
        plt.ylabel('Pocet čítaní', fontsize=15)
        plt.legend(loc='best', prop={'size': 15})
        plt.show()

    #for i in range(len(data)):
    #    print(data[i].l1_error)

    # usporiadame podla L1 normy
    #data.sort(key=lambda x: x.l1_error)
    
    # vypis 3 najlepsich distribucii podla L1 normy
    #print("Top 3 for L1 norm:")
    #print(distributions[data[0].distribution_num].name, ", shifted: ", data[0].shift)
    #print(distributions[data[1].distribution_num].name, ", shifted: ", data[1].shift)
    #print(distributions[data[2].distribution_num].name, ", shifted: ", data[2].shift)
    #print()

    # usporiadame podla L2 normy
    #data.sort(key=lambda x: x.l2_error)

    # vypis 3 najlepsich distribucii podla L2 normy
    #print("Top 3 for L2 norm:")
    #print(distributions[data[0].distribution_num].name, ", shifted: ", data[0].shift)
    #print(distributions[data[1].distribution_num].name, ", shifted: ", data[1].shift)
    #print(distributions[data[2].distribution_num].name, ", shifted: ", data[2].shift)

    # vratime data, ktoru potom pouzije funkcia DrawErrors
    data_for_drawing = DataForDrawing.DataForDrawing(data, num_of_shifts, distributions)
    return data_for_drawing
