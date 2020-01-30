import scipy.stats as st
import Discretize as disc
import DistributionData as dist_data
import matplotlib.pyplot as plt
import DataForDrawing as DataForDrawing
import DistributionTuple as dTuple


def drawErrors(data_for_drawing):
    # Funkcia, ktora vykresli pre kazdu distribuciu errory
    # ako parameter dostane objekt triedy DataForDrawing, ktory vytvoril drawDistribution
    
    data = data_for_drawing.data
    num_of_shifts = data_for_drawing.num_of_shifts
    distributions = data_for_drawing.distributions

    data.sort(key=lambda x: x.shift)
    data.sort(key=lambda x: x.distribution_num)
    # python sort is in place and stable

    temp_x = list(range(num_of_shifts))

    plt.figure(figsize=(15,8))
    plt.yscale("log")

    # vykreslime errory pri L1 norme pre shiftnutia
    for dist in range(len(distributions)):
        l1_errors = []
        for i in range(dist*num_of_shifts, dist*num_of_shifts + num_of_shifts):
            l1_errors.append(data[i].l1_error)
        plt.plot(temp_x, l1_errors, label = distributions[dist].name, color=distributions[dist].color, linewidth=3.0)

    plt.title(r'L1 Errors for shifts', fontsize=20)
    plt.xlabel('Shifts', fontsize=15)
    plt.ylabel('Errors', fontsize=15)
    plt.legend(loc='best', prop={'size': 15})
    plt.show()

    ######################

    plt.figure(figsize=(15,8))
    plt.yscale("log")
    
    # vykreslime errory pri L2 norme pre shiftnutia
    for dist in range(len(distributions)):
        l2_errors = []
        for i in range(dist*num_of_shifts, dist*num_of_shifts + num_of_shifts):
            l2_errors.append(data[i].l2_error)
        plt.plot(temp_x, l2_errors, label = distributions[dist].name, color=distributions[dist].color, linewidth=3.0)
                                                
    plt.title(r'L2 Errors for shifts', fontsize=20)
    plt.xlabel('Shifts', fontsize=15)
    plt.ylabel('Errors', fontsize=15)
    plt.legend(loc='best', prop={'size': 15})
    plt.show()