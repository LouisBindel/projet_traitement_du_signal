"""
Algorithm implementation
"""
import pickle
import numpy as np
import matplotlib.pyplot as plt

from scipy.io.wavfile import read
from scipy.signal import spectrogram
from skimage.feature import peak_local_max

# ----------------------------------------------------------------------------
# Create a fingerprint for an audio file based on a set of hashes
# ----------------------------------------------------------------------------


class Encoding:

   """
    Class implementing the procedure for creating a fingerprint 
    for the audio files

    The fingerprint is created through the following steps
    - compute the spectrogram of the audio signal
    - extract local maxima of the spectrogram
    - create hashes using these maxima

   """

   def __init__(self, fs, s):
      #Initialize Fingerprinting object with sampling frequency (fs) and sampled signal (s)
      self.fs = fs
      self.s = s
      self.spectrogram = None
      self.hashes = []
      


   def process(self):
      """

        To Do
        -----

        This function takes as input a sampled signal s and the sampling
        frequency fs and returns the fingerprint (the hashcodes) of the signal.
        The fingerprint is created through the following steps
        - spectrogram computation
        - local maxima extraction
        - hashes creation

        Implement all these operations in this function. Keep as attributes of
        the class the spectrogram, the range of frequencies, the anchors, the 
        list of hashes, etc.

        Each hash can conveniently be represented by a Python dictionary 
        containing the time associated to its anchor (key: "t") and a numpy 
        array with the difference in time between the anchor and the target, 
        the frequency of the anchor and the frequency of the target 
        (key: "hash")


        Parameters
        ----------

        fs: int
           sampling frequency [Hz]
        s: numpy array
           sampled signal
      """

      f, t, Sxx = spectrogram(self.s, self.fs)
      self.spectrogram = (f, t, Sxx)
      local_maxima = peak_local_max(Sxx,min_distance=50)
      hashes = []
      # Generate hashes
      delta_t=1
      delta_f=1500
      for point in local_maxima:
         t_anchor, f_anchor = point[1], point[0]
         for other_point in local_maxima:
            if other_point[1] > f_anchor and other_point-f_anchor<=delta_t and abs(other_point[0]-f_anchor)<delta_f :
               t_target, f_target = other_point[1], other_point[0]
               time_diff = t_target - t_anchor
               freq_diff = f_target - f_anchor
               hashes.append({"t": t_anchor, "hash": np.array([time_diff, f_anchor, f_target])})
      self.hashes = hashes



   def display_spectrogram(self):
      f, t, Sxx = self.spectrogram
      plt.pcolormesh(t, f, Sxx, shading='gouraud')
      plt.ylabel('Frequency [Hz]')
      plt.xlabel('Time [sec]')
      plt.show()



# ----------------------------------------------------------------------------
# Compares two set of hashes in order to determine if two audio files match
# ----------------------------------------------------------------------------

class Matching:

    """
    Compare the hashes from two audio files to determine if these
    files match

    Attributes
    ----------

    hashes1: list of dictionaries
       hashes extracted as fingerprints for the first audiofile. Each hash 
       is represented by a dictionary containing the time associated to
       its anchor (key: "t") and a numpy array with the difference in time
       between the anchor and the target, the frequency of the anchor and
       the frequency of the target (key: "hash")

    hashes2: list of dictionaries
       hashes extracted as fingerprint for the second audiofile. Each hash 
       is represented by a dictionary containing the time associated to
       its anchor (key: "t") and a numpy array with the difference in time
       between the anchor and the target, the frequency of the anchor and
       the frequency of the target (key: "hash")

    matching: numpy array
       absolute times of the hashes that match together

    offset: numpy array
       time offsets between the matches
    """

    def __init__(self, hashes1, hashes2):

        """
        Compare the hashes from two audio files to determine if these
        files match

        Parameters
        ----------

        hashes1: list of dictionaries
           hashes extracted as fingerprint for the first audiofile. Each hash 
           is represented by a dictionary containing the time associated to
           its anchor (key: "t") and a numpy array with the difference in time
           between the anchor and the target, the frequency of the anchor and
           the frequency of the target

        hashes2: list of dictionaries
           hashes extracted as fingerprint for the second audiofile. Each hash 
           is represented by a dictionary containing the time associated to
           its anchor (key: "t") and a numpy array with the difference in time
           between the anchor and the target, the frequency of the anchor and
           the frequency of the target
          
        """


        self.hashes1 = hashes1
        self.hashes2 = hashes2

        times = np.array([item['t'] for item in self.hashes1])
        hashcodes = np.array([item['hash'] for item in self.hashes1])

        # Establish matches
        self.matching = []
        for hc in self.hashes2:
             t = hc['t']
             h = hc['hash'][np.newaxis, :]
             dist = np.sum(np.abs(hashcodes - h), axis=1)
             mask = (dist < 1e-6)
             if (mask != 0).any():
                 self.matching.append(np.array([times[mask][0], t]))
        self.matching = np.array(self.matching)

        # TODO: complete the implementation of the class by
        # 1. creating an array "offset" containing the time offsets of the 
        #    hashcodes that match
        # 2. implementing a criterion to decide whether or not both extracts
        #    match
        self.offset = np.array(self.offset)
        match_threshold = 5  
        self.match = len(self.matching) >= match_threshold
             
    def display_scatterplot(self):

        """
        Display through a scatterplot the times associated to the hashes
        that match
        """
    
        plt.scatter(self.matching[:, 0], self.matching[:, 1])
        plt.show()


    def display_histogram(self):

        """
        Display the offset histogram
        """
    
        plt.hist(self.offsets, bins=100, density=True)
        plt.xlabel('Offset (s)')
        plt.show()


