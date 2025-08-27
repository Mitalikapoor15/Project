import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib import colormaps
from scipy.fft import fft, ifft, fftfreq
import math as m

def gaussian(qTime, maxTime):
    imp0 = 377
    start = int(0.20 * maxTime) # to start after 20% of the main signal 
    return np.exp(-(qTime - start) * (qTime - start)/ 10000.) /imp0

def cw(qTime):
    imp0 = 377
    start = 30
    return np.exp(-(qTime - start) * (qTime - start) / 100000000.) /imp0

def cw2(qTime, del_t, f0):
    
    f = f0 #resonant frequency (only frequency in this case)
    t = qTime * del_t
    A = 1.0 #amplitude

    phi0 = 0 #phase if required
    # s = A*np.sin(2*np.pi*f*t*dt + phi0)
    s_complex = A*np.sin(2.0*np.pi*f0*t + phi0)
    return s_complex

#Finding the resonant frequency of the CROW hoping for it to be around the resonant frequency of the single ring resonator (on-site frequency)
#We want to do this so that the system can be driven at the zero mode frequency and topological states can be accessed.

def zero_mode_freq(E, dt, f_ref, search_bw = 5e12): 
    #f_ref is the resonant frequency of the single ring around which the res freq of the CROW should land.
    
    N = len(E)
    t = np.arange(N) * dt
    
    start = int(0.30 * N) #using the last 70 percent of the signal so that steady state is reached.
    Eg = np.asarray(E[start:], dtype = np.complex128)
    Ng = len(Eg)
    w = np.hanning(Ng)
    
    #demodulating the baseband
    t_g = t[start:]
    Ebb = Eg * np.exp(-1j*2*np.pi*f_ref*t_g)
    
    #Taking the fft in detuned coordinates
    F = np.fft.fftshift(np.fft.fft(Ebb * w))
    freqs = np.fft.fftshift(np.fft.fftfreq(Ng, d=dt))
    detuning = freqs - f_ref
    mag = np.abs(F)
    
    
    #Searching for the frequency closest to f_ref i.e. close to zero detuning
    win = np.where(np.abs(detuning) <= search_bw)[0]
    k = win[np.argmax(max(win))]  #index of max in the window
    #quadratic interpolation around (k-1, k, k+1)
    if 0 < k < len(mag)-1:
        y1, y2, y3 = mag[k-1], mag[k], mag[k+1]
        denom = (y1 - 2*y2 + y3)
        delta = 0.5*(y1-y3)/denom if denom != 0 else 0.0
    else:
        delta = 0.0
    
    df = freqs[1] - freqs[0]
    detuning_peak = detuning[k] + delta*df
    f_zero = f_ref + detuning_peak
    df_fft = 1.0 /(Ng *dt)
    
    return f_zero, df_fft, detuning, mag
    

def cosMod(qTime, complex_signal, f0, sigma, del_t):
    dt=del_t
    # sigma= 10e-15 # Width
    phase=0.0
    # f0 = 300
    # Time array
    t = qTime * dt
    t0 = 4 * sigma  # Center time
    
    # Gaussian envelope
    g = np.exp(-0.5 * ((t - t0) / sigma)**2)
    
    # Signal
    if complex_signal:
        signal = g * np.exp(1j * (2 * np.pi * f0 * t + phase))
    else:
        signal = g * np.cos(2 * np.pi * f0 * t + phase)
    
    return signal

def create_colored_arc(center, radius, values, angle_range, cmap='plasma'):
    """
    Creates a colored arc (full or partial ring).
    angle_range: tuple (start_angle, end_angle) in radians.
    """
    N = len(values)
    theta = np.linspace(angle_range[0], angle_range[1], N + 1)
    x = center[0] + radius * np.cos(theta)
    y = center[1] + radius * np.sin(theta)

    # Make line segments
    segments = [[[x[i], y[i]], [x[i+1], y[i+1]]] for i in range(N)]
    segments = np.array(segments)

    # Normalize and map to colormap
    norm_values = (values - np.min(values)) / (np.ptp(values) + 1e-8)
    colors = colormaps[cmap](norm_values)

    return LineCollection(segments, colors=colors, linewidths=14), x, y

def plot_field_ring(ez_tab_tp, N_rings):
    num_rings = N_rings - 2  # Number of full rings in the middle (will be 2 less than the total number of rings as they will act as i/p and o/p)
    N_seg = ez_tab_tp.shape[0]
    print(N_seg)
    
    radius = 1.0
    spacing = 2.3  # spacing between ring centers
    cmap = 'Reds'

    fig, ax = plt.subplots(figsize=(12, 4))
    all_x = []
    all_y = []

    # Left half-ring (start)
    center = (0, 0)
    values = np.real(np.hstack((ez_tab_tp[0,:], ez_tab_tp[1,:])) )

    arc, x, y = create_colored_arc(center, radius, values, (-np.pi/2, np.pi/2), cmap)
    ax.add_collection(arc)
    all_x.extend(x)
    all_y.extend(y)

    # Full rings
    for i in range(num_rings):
        center = ((i + 1) * spacing, 0)
        if i%2 == 0:
        # values = np.real(np.hstack((ez_tab_tp[2*(i+1),:], ez_tab_tp[2*(i+1)+1,:])))
        # arc, x, y = create_colored_arc(center, radius, values, (0, 2*np.pi), cmap)
        # arc, x, y = create_colored_arc(center, radius, values, (2*np.pi,0), cmap)
            arc, x, y = create_colored_arc(center, radius, values, (0*np.pi, 2*np.pi), cmap)
        else:
            arc, x, y = create_colored_arc(center, radius, values, (-1*np.pi, 1*np.pi), cmap)

        ax.add_collection(arc)
        all_x.extend(x)
        all_y.extend(y)

    # Right half-ring (end)
    center = ((num_rings + 1) * spacing, 0)
    if num_rings%2 == 0: #even no of rings in the middle
        values = np.real( np.hstack((ez_tab_tp[N_seg-2,:], ez_tab_tp[N_seg-1,:])) )
        arc, x, y = create_colored_arc(center, radius, values,(1*np.pi/2, 3*np.pi/2) , cmap)
    else:
        values = np.real( np.hstack((ez_tab_tp[N_seg-1,:], ez_tab_tp[N_seg-2,:])) )
        arc, x, y = create_colored_arc(center, radius, values,(3*np.pi/2, 1*np.pi/2) , cmap)
        
    ax.add_collection(arc)
    all_x.extend(x)
    all_y.extend(y)
    
    #need to fix the last (right) port for even number of rings in the middle.

    # Set axis limits based on all arc coordinates
    margin = 0.5
    ax.set_xlim(min(all_x) - margin, max(all_x) + margin)
    ax.set_ylim(min(all_y) - margin, max(all_y) + margin)
    ax.set_aspect('equal')
    ax.axis('off')

    plt.title(f'{num_rings} Coupled Optical Rings with Half Arcs at Ends')
    plt.tight_layout()
    plt.show()

def Sources(N_rings):
    seg_no = N_rings*2
    s = np.zeros((seg_no,2), dtype=int)
    s1 = np.zeros((N_rings,2), dtype=int) #forward transmission
    s2 = np.zeros((N_rings,2), dtype=int) #backward propagation
    s1[0][1] = -1
    s1[0][0] = -1
    for i in range (N_rings):
        if i > 0:
            s1[i][0] = (2*i + 1)  #odd indices belong to the s2 segments which are responsible for back propagation
            s1[i][1] = 2*i - 2 #because if it is 0, we will have -2 index which is not what we want.
        
        s2[i][0] = 2*i  #for the transmission of odd segments  
        s2[i][1] = 2*i  + 3
        if i==(N_rings-1):
            s2[i][1] =  2*i - 1
        
    for i in range(N_rings):
        s[2*i][0] = s1[i][0]
        s[2*i][1] = s1[i][1]
        s[2*i + 1][:] = s2[i][:]

    s[2*N_rings-1][0] = -2 
    s[2*N_rings-1][1] = -2 
    return s

def Couplings(N_rings, tau):
    c = np.zeros((N_rings*2,2), dtype=complex)
    t = tau
    k = 1j* m.sqrt(1-t**2)
    c[:][:] = (t,k)
    return c
    
def SSH_Couplings(N_rings, tau_alt):
    c = np.zeros((N_rings*2,2), dtype=complex)
    
    t = tau_alt
    k = np.zeros(len(t), dtype=complex)
    for i in range(len(t)):
        k[i] = 1j* m.sqrt(1-t[i]**2)
    for i in range(N_rings*2):
        if i%2 == 0:
            c[i][:] = (t[0],k[0])
        else:
            c[i][:] = (t[1],k[1])
    return c