import numpy as np
import csv
import matplotlib.pyplot as plt
import math
import pandas as pd  
import seaborn as sns
from DrawDistributions import draw
from scipy.stats import entropy


######## BEDGRAPH FUNCTIONS ########

def read_bedgraph(path):
    """Reads the bedgraph from the path

    Returns:
        array containing on each index 4 entries:
        chromName  chromStartPosition  chromEndPosition  dataValue
    """
    with open(path) as bed: 
        bed_reader = csv.reader(bed, delimiter='\t')
        bed_columns = list(zip(*bed_reader))
    return bed_columns
    
    
def get_bedG_ranges(bed_columns):
    """Returns numpy array containing the 2nd and 3rd column 
    (chromStartPosition, chromEndPosition) from the bedgraph's columns """
    return np.array([(int(bed_columns[1][i]), int(bed_columns[2][i])) for i in range(len(bed_columns[1]))])


def get_bedG_ranges_size(bed_columns):
    """Returns numpy array containing the size of the ranges of the bedgraph's columns:
    3rd column - 2nd column =  (chromEndPosition - chromStartPosition)"""
    return np.array([(int(bed_columns[2][i]) - int(bed_columns[1][i])) for i in range(len(bed_columns[1]))])


def get_bedG_scores(bed_columns):
    """Returns numpy array containing the 4th column 
    (dataValue) from the bedgraph's columns """
    return np.array([int(bed_columns[3][i]) for i in range(len(bed_columns[3]))])


# vrati najvacsiu chromEndPosition
def get_bedG_maxPosition(bed_columns):
    """Returns the ending position of the bedgraph (biggest chromEndPosition)"""
    return bed_columns[2][len(bed_columns[0])-1]


#function inserts bedG_range_sizes[i] times value bedG_scores[i] 
def filler(bedG_range_sizes, bedG_scores):
    """Function returns numpy array where
    bedG_scores[i] is inserted bedG_range_sizes[i] times"""
    return np.repeat(bedG_scores, bedG_range_sizes)	
    
    
    

######## DISTRIBUTION DRAWING FUNCTIONS ########

def drawHistogramsForReads(ranges, score, title='Title', x_axis_name='x axis', y_axis_name='y axis'):
    """Draws distributions 
    To change the distributions go to: DrawDistributions.py - function draw """
    if len(ranges) != len(score):
        raise ValueError('Wrong files, the length of arrays should be equal')

    # Zadaj limitu pre os x:
    max_x=100

    #spravim histogram a zobrazim ho pomocou draw
    histo = {}
    for i in range(len(ranges)) :
        if int(score[i]) in histo:
            histo[int(score[i])] += (int(ranges[i][1])-int(ranges[i][0]))
        else:
            histo[int(score[i])] = (int(ranges[i][1])-int(ranges[i][0]))
    # treba vymazat 0, inak to pada
    if 0 in histo:
        del histo[0]
    data = draw(histo, max_x, title, x_axis_name, y_axis_name)
    #print(len(histo))    
    print('Maximum number of reads:', max(histo))



######## HEATMAP ########

def make_heatmap(x, y, title='Heatmap', colorbar_title = 'values', x_axis_name='x', y_axis_name='y', are_integers=True):
    """Draws a 2D heatmap
    are_integers = bins have to be changed for data with non-integer values on y axis
    """
    #predpokladam, ze tie mena chromozom su v rovnakom poradi a velkost suboru je rovanaka

    #bins = nasekat rovnomerne

    plt.figure(figsize=(20,10))  
    if are_integers:
        plt.hist2d(x, y,bins=(range(max(x)), range(max(y))), cmap='plasma')
    else:
        plt.hist2d(x, y,bins=(range(max(x)), np.linspace(min(y),max(y))), cmap='plasma')


    cb = plt.colorbar()
    cb.set_label(colorbar_title, fontsize=15)

    plt.title(title, fontsize=20)
    plt.xlabel(x_axis_name, fontsize=15)
    plt.ylabel(y_axis_name, fontsize=15)

    plt.show()




######## FASTA FUNCTIONS ########

def read_fasta_file(filename):
    """Converts fasta file to list. 
    The returned list contains only characters of the sequence"""
    data = []
    with open(filename) as fasta:
        for line in fasta: 
            if not line.startswith(">"):
                for c in line.strip():
                    data.append(c)
    return data


def create_values_from_fasta(fasta_data, kmer_length, function):
    """Function creates a list of values. 
    These values are created by applying a function on a window (sublist). 
    We save the value to the middle position (middle of k-mer length)"""
    if (kmer_length % 2 != 1):
        raise ValueError('The length of k-mer should be odd!')
    
    values = []

    half = int(kmer_length/2) 
    
    for i in range(len(fasta_data) - kmer_length +1):  # +1lebo ideme aj za posledne pismenko
        window = fasta_data[i:i+kmer_length]

        #ak su to kraje, tak napln half-krat values
        if i == 0 or i == (len(fasta_data) - kmer_length):
            # !!! pomale: values += [function(window)] * half
            values.extend(function(window) for i in range(half))
            
        #ak su to vnutorne casti
        if i >= 0 and i <= (len(fasta_data) - kmer_length):
            values.append(function(window))  
    return values
    

def count_GC(string):
    """Simple function to count G and C in list of characters"""
    count = 0
    for c in string:
        if(c == 'G' or c == 'C'):
            count += 1
    return count



def count_entropy(string):
    """Function that counts Shannon Entropy = 
    estimation of the average amount of information stored in a random variable"""
    length = len(string)
    if (length == 0):
        raise ValueError('Empty window!')
    probabilities = []
    probabilities.append(string.count('C')/length)
    probabilities.append(string.count('G')/length)
    probabilities.append(string.count('A')/length)
    probabilities.append(string.count('T')/length)

    return entropy(probabilities, base=2)


    
######## BOXPLOT ########

def make_boxplot(x, y, title='Boxplot',  x_axis_name='x', y_axis_name='y'):
    """Creates a boxplot"""
    plt.figure(figsize=(50,30))
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=30)
    plt.title(title, fontsize=80)
    plt.xlabel(x_axis_name, fontsize=50)
    plt.ylabel(y_axis_name, fontsize=50)
    sns.boxplot(x=x, y=y)
    plt.show
