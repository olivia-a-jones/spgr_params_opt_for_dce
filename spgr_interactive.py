import numpy as np
import matplotlib.pyplot as plt
from ipywidgets import interact, interactive, fixed, interact_manual, Layout
import ipywidgets as widgets
from IPython.display import display
import seaborn as sns

# Set bigger font sizes globally
plt.rcParams.update({'font.size': 12,
                    'axes.labelsize': 14,
                    'axes.titlesize': 16,
                    'xtick.labelsize': 12,
                    'ytick.labelsize': 12})

# Constants and conversions remain the same
R1_TO_RELAXIVITY = 3.4  # s^-1 mM^-1
MS_TO_S = 1e-3  # conversion from ms to s

def calculate_signal(fa, TR_ms, T1_ms, M0=1):
    """
    Calculate SPGR signal.
    
    Parameters:
    fa (float): Flip angle in radians
    TR_ms (float): Repetition time in milliseconds
    T1_ms (float): T1 relaxation time in milliseconds
    M0 (float): Equilibrium magnetization (default=1)
    """
    TR_s = TR_ms * MS_TO_S
    T1_s = T1_ms * MS_TO_S
    E1 = np.exp(-TR_s / T1_s)
    return M0 * np.sin(fa) * (1 - E1) / (1 - E1 * np.cos(fa))

def dSdFA(fa, TR_ms, T1_ms, M0=1):
    """Calculate derivative of signal with respect to flip angle."""
    TR_s = TR_ms * MS_TO_S
    T1_s = T1_ms * MS_TO_S
    E1 = np.exp(-TR_s / T1_s)
    num = M0 * (1 - E1) * (np.cos(fa) * (1 - E1 * np.cos(fa)) + E1 * np.sin(fa) * np.sin(fa))
    den = (1 - E1 * np.cos(fa))**2
    return num / den

def dSdR1(fa, TR_ms, T1_ms, M0=1):
    """Calculate derivative of signal with respect to R1."""
    TR_s = TR_ms * MS_TO_S
    T1_s = T1_ms * MS_TO_S
    E1 = np.exp(-TR_s / T1_s)
    R1 = 1 / T1_s
    num = M0 * np.sin(fa) * TR_s * E1 * (1 - np.cos(fa))
    den = (1 - E1 * np.cos(fa))**2
    return num / den

def find_ernst_angle(TR_ms, T1_ms):
    """Find Ernst angle numerically."""
    angles = np.linspace(0.001, np.pi, 1800)
    signals = calculate_signal(angles, TR_ms, T1_ms)
    return angles[np.argmax(signals)]

def find_optimal_R1_angle(TR_ms, T1_ms):
    """Find optimal angle for R1 precision."""
    angles = np.linspace(0.001, np.pi, 1800)
    sensitivity = np.abs(dSdR1(angles, TR_ms, T1_ms))
    return angles[np.argmax(sensitivity)]

def T1_from_R1(R1_s):
    """Convert R1 (s^-1) to T1 (ms)"""
    return 1000 / R1_s  # Convert to ms

def R1_from_concentration(C_mM, R10_s, r1=R1_TO_RELAXIVITY):
    """Calculate R1 from concentration using linear model"""
    return R10_s + r1 * C_mM

def plot_heatmap(ax, data, TR_range, T1_range, title, cbar_label):
    """Create a nicely formatted heatmap using seaborn."""
    # Create the heatmap
    g = sns.heatmap(data, 
                ax=ax,
                cmap='jet',
                xticklabels=TR_range,
                yticklabels=T1_range,
                cbar_kws={'label': cbar_label,
                         'shrink': 0.85,    # Adjusted colorbar height
                         'aspect': 40,        # Make colorbar thinner
                         'pad': 0.02},        # Adjust spacing
                fmt='.1f',                    # Format for annotations
                annot=True,                   # Show annotations
                annot_kws={'size': 12},       # Larger annotation text
                square=True)                  # Keep square aspect ratio
    
    # Customize appearance
    ax.set_xlabel('TR [ms]', fontsize=14)
    ax.set_ylabel('T10 [ms]', fontsize=14)  # Changed to T1₀
    ax.set_title(title, fontsize=16, pad=20)
    
    # Remove cell borders
    ax.collections[0].set_linewidth(0)
    
    # Make axis labels bigger and rotate y-axis labels
    ax.tick_params(labelsize=12)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    
    # Adjust colorbar label size after creation
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=12)
    cbar.set_label(cbar_label, size=14)
    
    return ax

def plot_spgr_comprehensive(TR_ms, T10_ms, C_mM_max=0.005):
    """Create comprehensive SPGR analysis plots with current parameters."""
    fig = plt.figure(figsize=(20, 25))
    gs = plt.GridSpec(4, 2, height_ratios=[2, 2, 1, 1], hspace=0.4, wspace=0.3)
    
    # Plot 1: Signal and dS/dR1 vs flip angle
    ax1a = fig.add_subplot(gs[0, 0])
    ax1b = ax1a.twinx()
    
    angles_deg = np.linspace(1, 90, 90)
    angles_rad = angles_deg * np.pi / 180
    
    signal = calculate_signal(angles_rad, TR_ms, T10_ms)
    dsr1 = dSdR1(angles_rad, TR_ms, T10_ms)
    
    l1 = ax1a.plot(angles_deg, signal, 'b-', label='Signal')
    l2 = ax1b.plot(angles_deg, dsr1, 'r-', label='dS/dR1')
    
    ernst_angle = find_ernst_angle(TR_ms, T10_ms) * 180 / np.pi
    opt_angle = find_optimal_R1_angle(TR_ms, T10_ms) * 180 / np.pi
    ax1a.axvline(ernst_angle, color='g', linestyle='--', alpha=0.5,
                 label=f'Ernst Angle ({ernst_angle:.1f}°)')
    ax1a.axvline(opt_angle, color='k', linestyle='--', alpha=0.5,
                 label=f'Optimal R1 Angle ({opt_angle:.1f}°)')
    
    ax1a.set_xlabel('Flip Angle (degrees)')
    ax1a.set_ylabel('Signal (a.u.)')
    ax1b.set_ylabel('dS/dR1 (a.u. * s)')
    
    lines = l1 + l2 + [plt.Line2D([0], [0], color='g', linestyle='--'),
                       plt.Line2D([0], [0], color='k', linestyle='--')]
    labels = [l.get_label() for l in l1 + l2] + ['Ernst Angle', 'Optimal R1 Angle']
    ax1a.legend(lines, labels)
    ax1a.set_title(f'Signal and R1 Sensitivity (TR={TR_ms}ms, T1₀={T10_ms}ms)')
    
    # Plot 2: Signal vs R1 for different flip angles
    ax2 = fig.add_subplot(gs[0, 1])
    R10_s = 1 / (T10_ms * MS_TO_S)
    conc_range = np.linspace(0, C_mM_max, 100)
    R1_range_s = R1_from_concentration(conc_range, R10_s)
    T1_range_ms = T1_from_R1(R1_range_s)
    
    test_angles = [1, 2, 4, 6, 10, 15, 20]
    for fa_deg in test_angles:
        fa_rad = fa_deg * np.pi / 180
        signal = calculate_signal(fa_rad, TR_ms, T1_range_ms)
        ax2.plot(conc_range * 1000, signal, label=f'{fa_deg}°')  # Convert to μM
    
    ax2.set_xlabel('Concentration (μM)')
    ax2.set_ylabel('Signal (a.u.)')
    ax2.set_title('Signal vs Concentration for Different Flip Angles')
    ax2.legend()
    
    # Plot 3: Optimal angle vs TR for GM/WM
    ax3 = fig.add_subplot(gs[1, :])
    TRs = np.logspace(0, 2, 50)
    T1s = [850, 1250]  # GM and WM at 3T
    
    for t1 in T1s:
        optimal_angles = []
        ernst_angles = []
        for tr in TRs:
            optimal = find_optimal_R1_angle(tr, t1) * 180 / np.pi
            ernst = find_ernst_angle(tr, t1) * 180 / np.pi
            optimal_angles.append(optimal)
            ernst_angles.append(ernst)
        
        label = 'Gray Matter' if t1 > 1000 else 'White Matter'
        ax3.plot(TRs, ernst_angles, '--', label=f'{label} Ernst (T1={t1}ms)')
        ax3.plot(TRs, optimal_angles, '-', label=f'{label} Optimal R1 (T1={t1}ms)')
    
    ax3.set_xscale('log')
    ax3.set_xlabel('TR (ms)')
    ax3.set_ylabel('Flip Angle (degrees)')
    ax3.set_title('Optimal Angles vs TR for Brain Tissues')
    ax3.legend()
    ax3.grid(True)
    
    # Create square axes for heatmaps
    ax4 = fig.add_subplot(gs[2:, 0])
    ax5 = fig.add_subplot(gs[2:, 1])
    
    # Generate data for heatmaps
    # For the heatmaps, update TR and T1 ranges based on slider values
    # For heatmaps, update TR and T1 ranges based on slider values
    TR_range = np.linspace(max(1, TR_ms-10), min(20, TR_ms+10), 10)
    T1_range = np.linspace(max(800, T10_ms-400), min(2000, T10_ms+400), 10)
    
    ernst_angles = np.zeros((len(T1_range), len(TR_range)))
    optimal_angles = np.zeros((len(T1_range), len(TR_range)))
    
    for i, t1 in enumerate(T1_range):
        for j, tr in enumerate(TR_range):
            ernst_angles[i, j] = find_ernst_angle(tr, t1) * 180 / np.pi
            optimal_angles[i, j] = find_optimal_R1_angle(tr, t1) * 180 / np.pi
    
    # Clear previous content
    ax4.clear()
    ax5.clear()
    
    # Create enhanced heatmaps
    plot_heatmap(ax4, ernst_angles, 
                [f'{tr:.1f}' for tr in TR_range],
                [f'{t1:.0f}' for t1 in T1_range],
                'Optimal Flip Angle for Maximum Signal',
                'Ernst Angle (degrees)')
    
    plot_heatmap(ax5, optimal_angles,
                [f'{tr:.1f}' for tr in TR_range],
                [f'{t1:.0f}' for t1 in T1_range],
                'Optimal Flip Angle for Maximum Sensitivity to Change in T1',
                'Optimal R1 Angle (degrees)')
    
    #plt.tight_layout()
    plt.show()

# Create interactive widgets with explicit units
TR_slider = widgets.FloatSlider(
    value=10,
    min=1,
    max=50,
    step=1,
    description='TR (ms):',
    style={'description_width': 'auto'},
    continuous_update=False)

T10_slider = widgets.FloatSlider(
    value=1000,
    min=800,
    max=2000,
    step=50,
    description='T1₀ (ms):',
    style={'description_width': 'auto'},
    continuous_update=False)

C_max_slider = widgets.FloatSlider(
    value=0.005,
    min=0.001,
    max=0.1,
    step=0.001,
    description='Max Conc. (mM):',
    style={'description_width': 'auto'},
    continuous_update=False,
    readout_format='.3f')

# Create interactive plot with layout adjustments
interactive_plot = interactive(
    plot_spgr_comprehensive, 
    TR_ms=TR_slider,
    T10_ms=T10_slider,
    C_mM_max=C_max_slider)

# Add layout adjustments
interactive_plot.children[-1].layout.height = '1000px'

# When you want to display the interactive plot in a notebook, run:
if __name__ == "__main__":
    display(interactive_plot)