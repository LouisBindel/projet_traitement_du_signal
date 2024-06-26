"""
Description
"""
import os
import random
import numpy as np
import matplotlib.pyplot as plt

from scipy.io.wavfile import read
from algorithm import *

# ----------------------------------------------
# Run the script
# ----------------------------------------------
if __name__ == '__main__':

    # 1: Load the database
    with open('songs.pickle', 'rb') as handle:
        database = pickle.load(handle)

    # 2: Create an instance of the class Encoder
    # Insert code here
    nperseg=128
    noverlap=32
    min_distance=50
    time_window=1.
    freq_window=1500
    encoder = Encoding(nperseg, noverlap, min_distance,time_window,freq_window)

    

    # 3: Randomly get an extract from one of the songs of the database
    songs = [item for item in os.listdir('./samples') if item[:-4] != '.wav']
    song = random.choice(songs)
    print('Selected song: ' + song[:-4])
    filename = './samples/' + song

    fs, s = read(filename)
    tmin = int(50*fs) # We select an extract starting at 50s ...
    duration = int(10*fs) # ... which lasts 10s

    # 4: Use the encoder to extract a fingerprint of the sample
    encoder.process(fs, s[tmin:tmin + duration])
    hashes = encoder.hashes

    # 5: Using the class Matching, compare the fingerprint to all the 
    # fingerprints in the database
    for song in os.listdir('./samples') :
        encoder2 = Encoding(nperseg, noverlap, min_distance,time_window,freq_window)
        filename = './samples/' + song
        fs, s = read(filename)
        encoder2.process(fs, s)
        hashes2 = encoder.hashes
        match = Matching(hashes,hashes2)
        match.display_scatterplot()
        match.display_histogram()






