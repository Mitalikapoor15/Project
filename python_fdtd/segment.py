import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft, fftfreq
import math as m

class Segment:
    def __init__(self, SIZE):
        self.ez = np.zeros([SIZE], dtype=complex)
        self.hy = np.zeros([SIZE], dtype=complex)
        self.sx = 49
        self.size = SIZE
        
    def hy_update(self):
        imp0 = 377
        self.hy[self.size-1] = self.hy[self.size-2] # simple ABC for hy
        for j in range(self.size-1):
            self.hy[j] = self.hy[j] + (self.ez[j+1] - self.ez[j]) / imp0;
            
    def ez_update(self):
        imp0 = 377
        for j in range(1,self.size):
            self.ez[j] = self.ez[j] + (self.hy[j] - self.hy[j-1]) * imp0;
        self.ez[0] = self.ez[1]  # simple ABC for ez

    def hy_sources(self, s, sx):
        imp0 = 377
        self.sx = sx
        self.hy[self.sx] -= s /imp0 
        
    def ez_sources(self, s,sx):   
        self.sx = sx
        self.ez[self.sx] += s  