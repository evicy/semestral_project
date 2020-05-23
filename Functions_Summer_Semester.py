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
        numpy array where elements look like this:
        [ chromStartPosition  chromEndPosition  dataValue]
    """
    return np.loadtxt(path, usecols=range(1,4), dtype=int)

    
    
def get_bedG_ranges(bed_columns):
    """Returns numpy array containing (chromStartPosition, chromEndPosition) 
    from the bedgraph's columns 
    Looks like this:[ [chromStartPosition  chromEndPosition] ...]
    """
    return bed_columns[: , 0:2]


def get_bedG_ranges_size(bed_columns):
    """Returns numpy array containing the size of the ranges of the bedgraph's columns
    Looks like this: [ chromEndPosition-chromStartPosition, ...] 
    """
    ranges = get_bedG_ranges(bed_columns)
    result = np.zeros(len(ranges), dtype=int)
    for i in range(len(ranges)):
        result[i] = int(ranges[i][1]-ranges[i][0])    
    return result


def get_bedG_scores(bed_columns):
    """Returns numpy array containing the values(number of reads) 
    Looks like this: [value, ...] """
    return bed_columns[:,2]


# vrati najvacsiu chromEndPosition
def get_bedG_maxPosition(bed_columns):
    """Returns the ending position of the bedgraph (biggest chromEndPosition)"""
    return bed_columns[len(bed_columns)-1][1]



#function inserts bedG_range_sizes[i] times value bedG_scores[i] 
def filler(bedG_range_sizes, bedG_scores):
    """Function returns numpy array where
    bedG_scores[i] is inserted bedG_range_sizes[i] times"""
    return np.repeat(bedG_scores, bedG_range_sizes)
    
    
############################### NOT USED ##################################################### 
# something like the filler, but the avarage of step values will be stored in the result array
# the resulting array will be step times smaller
def compressed_filler(bedG_range_sizes, bedG_scores):
    step = 40
    read_from_prev = 0
    can_read = step
   
    result = np.zeros(math.ceil(np.sum(bedG_range_sizes)/step), dtype=int)
    print(math.ceil(np.sum(bedG_range_sizes)/step))
    res_index = 0
   
    for i in range(len(bedG_range_sizes)):
        #print("i = ", i)
        while read_from_prev > 0:
            #print("while: read_from_prev = ", read_from_prev)
            read = min(step, read_from_prev)
            result[res_index] += (read * bedG_scores[i-1])
            read_from_prev -= read
            if read_from_prev > 0:
                res_index +=1
            if read_from_prev == 0:
                if read != step:
                    can_read -= read
                   
        if (can_read - bedG_range_sizes[i]) >= 0:
            #print("if: can_read  = ", can_read, " bedG_range_sizes[i] = ", bedG_range_sizes[i])
            result[res_index] += bedG_range_sizes[i] * bedG_scores[i]
            can_read -= bedG_range_sizes[i]
            if can_read == 0:
                can_read = step
                res_index += 1

        else:
            #print("else: can_read  = ", can_read, " bedG_range_sizes[i] = ", bedG_range_sizes[i])
            #print("read_from_prev = ", read_from_prev)
            result[res_index] += can_read * bedG_scores[i]
            read_from_prev = bedG_range_sizes[i] - can_read
            res_index += 1
            can_read = step
        #print("\n")
       
    for i in range(len(result)):
        if i == len(result)-1 and (len(bedG_range_sizes) % step != 0) :
            print(result[i])
            result[i] = result[i]/(step-(len(bedG_range_sizes) % step))
            print("result", result[i])
            print((step-(len(bedG_range_sizes) % step)))
        else:
            result[i] = result[i]/step
    return result
############################### NOT USED ##################################################### 
    

    

######## DISTRIBUTION DRAWING FUNCTIONS ########

def drawHistogramsForReads(ranges, score, title='Title', x_axis_name='x axis', y_axis_name='y axis', max_x=100):
    """Draws distributions 
    To change the distributions go to: DrawDistributions.py - function draw """
    if len(ranges) != len(score):
        raise ValueError('Wrong files, the length of arrays should be equal')

    # Zadaj limitu pre os x:
    max_x=max_x

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

# if memory error occures then the problem may be because of the number of bins
# specify a small value, for example bins=(300,300)
def make_heatmap(x, y, title='Heatmap', colorbar_title = 'values', bins=-1, x_axis_name='x', y_axis_name='y', x_axis_are_integers=True, y_axis_are_integers=True):
    """Draws a 2D heatmap
    are_integers = bins have to be changed for data with non-integer values on y axis
    """
    plt.figure(figsize=(20,10)) 
    
    #ak bin neni nastaveny
    if bins == -1:
        # moze nastat memory error

        if x_axis_are_integers:
            x_bin = range(max(x))

        if x_axis_are_integers==False:
            x_bin = np.linspace(min(x),max(x))
            
        if y_axis_are_integers:
            y_bin = range(max(y))
            
        if y_axis_are_integers==False:
            y_bin = np.linspace(min(y),max(y))
    
        bins = (x_bin, y_bin)

    
    plt.hist2d(x, y,bins=bins, cmap='plasma')


    cb = plt.colorbar()
    cb.set_label(colorbar_title, fontsize=15)

    plt.title(title, fontsize=20)
    plt.xlabel(x_axis_name, fontsize=15)
    plt.ylabel(y_axis_name, fontsize=15)

    plt.show()



######## FASTA FUNCTIONS ########

def read_fasta_file(filename):
    data = []
    with open(filename) as fasta:
        chrom = []
        for line in fasta: 
            if line.startswith(">"):
                if len(chrom) != 0:
                    data.append(chrom)
                chrom = []
            else:
                for c in line.strip():
                    chrom.append(c)
        data.append(chrom)
    #print(len(data))
    return data

def create_values_from_fasta(fasta_data, kmer_length, function):
    if (kmer_length % 2 != 1):
        raise ValueError('The length of k-mer should be odd!')
    
    values = []

    half = int(kmer_length/2) 
    
    for chrom in fasta_data:
        for i in range(len(chrom) - kmer_length +1):  # +1lebo ideme aj za posledne pismenko
            window = chrom[i:i+kmer_length]

            #ak su to kraje, tak napln half-krat values
            if i == 0 or i == (len(chrom) - kmer_length):
                values.extend(function(window) for i in range(half))

            #ak su to vnutorne casti
            if i >= 0 and i <= (len(chrom) - kmer_length):
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

    return round(entropy(probabilities, base=2),2)

    
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
