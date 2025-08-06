import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib import colormaps
from scipy.fft import fft, ifft, fftfreq
import math as m

def gaussian(qTime):
    imp0 = 377
    return np.exp(-(qTime - 30.) * (qTime - 30.) / 100.) /imp0

def cw(Time):
    imp0 = 377
    return np.exp(-(Time - 30.) * (Time - 30.) / 100000000.) /imp0

def cosMod(qTime, maxTime):
    dt = 10e-17  # Time step (s) //added
    t = qTime * dt
    A = 1
    f0 = 300e12
    t0 = 4 * dt * maxTime / 50  # Center time
    tau = 2 * dt * maxTime / 50 # Width
    source = A * np.cos(2 * np.pi * f0 * (t - t0)) * np.exp(-((t - t0)**2) / (tau**2))
    source = A * np.exp( 1j* 2 * np.pi * f0 * (t - t0)) * np.exp(-((t - t0)**2) / (tau**2))
    return source

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

def plot_field_ring(ez_tab_tp):
    num_rings = 1  # Number of full rings in the middle
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
        values = np.real(np.hstack((ez_tab_tp[2*(i+1),:], ez_tab_tp[2*(i+1)+1,:])))
        arc, x, y = create_colored_arc(center, radius, values, (0, 2*np.pi), cmap)
        ax.add_collection(arc)
        all_x.extend(x)
        all_y.extend(y)

    # Right half-ring (end)
    center = ((num_rings + 1) * spacing, 0)
    values = np.real( np.hstack((ez_tab_tp[N_seg-2,:], ez_tab_tp[N_seg-1,:])) )
    arc, x, y = create_colored_arc(center, radius, values,(np.pi/2, 3*np.pi/2) , cmap)
    ax.add_collection(arc)
    all_x.extend(x)
    all_y.extend(y)

    # Set axis limits based on all arc coordinates
    margin = 0.5
    ax.set_xlim(min(all_x) - margin, max(all_x) + margin)
    ax.set_ylim(min(all_y) - margin, max(all_y) + margin)
    ax.set_aspect('equal')
    ax.axis('off')

    plt.title(f'{num_rings} Coupled Optical Rings with Half Arcs at Ends')
    plt.tight_layout()
    plt.show()
