import numpy as np
import datetime as dt

number_of_points = 100000000
a = np.random.randn(number_of_points) # generate random data
sliding_window = 2000
sliding_window_list = np.arange(0,number_of_points,sliding_window)

running_mean_ist = []
bins = number_of_points/2000

def sliding_mean_hist(a, sliding_window):
    print('my way')
    bins = number_of_points/2000
    for i in range(sliding_window):
        print(i)
        hist, edges = np.histogram(a[i:], bins=int(bins))


def calculate_sliding_mean(data, window_size):
    print('calculate_sliding_mean')
    """
    Calculate the sliding mean of an array.

    Parameters:
    data (numpy.ndarray): The input array.
    window_size (int): The size of the sliding window.

    Returns:
    numpy.ndarray: An array of sliding means.
    """
    if window_size <= 0:
        raise ValueError("Window size must be a positive integer.")
    
    kernel = np.ones(window_size) / window_size
    sliding_means = np.convolve(data, kernel, mode='valid')
    
    return sliding_means

calculate_sliding_mean(a, sliding_window)