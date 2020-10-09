"""
Helper functions for plotting consistency across analysis / notebooks
"""
import numpy as np
import matplotlib as mpl
from . import proc_utils
from . import filter_tools


def set_mpl():
    # Use these settings and create as many subplots on one figure as possible, calling fig.tight_layout()
    # This will mean they are uniformly "tightened" and require less tweaking in illustrator

    # Set all default values.
    mpl.rcdefaults()

    # FONT SETTINGS
    mpl.rcParams['pdf.fonttype'] = 42                   # Makes the fonts editable inside of illustrator
    mpl.rcParams['ps.fonttype'] = 42                    # Makes the fonts editable inside of illustrator
    mpl.rcParams['font.sans-serif'] = 'Arial'           # Arial or Helvetica required by cell press
    mpl.rcParams['font.size'] = 14.0                    # For cell graphical abstracts "Font size: 12–16 points"
    mpl.rcParams['axes.labelsize'] = 14.0               # Can be set distinctly from font.size, otherwise uses that
    mpl.rcParams['axes.titlesize'] = 14.0               # Can be set distinctly from font.size, otherwise uses that
    mpl.rcParams['xtick.labelsize'] = 14.0              # Can be set distinctly from font.size, otherwise uses that
    mpl.rcParams['ytick.labelsize'] = 14.0              # Can be set distinctly from font.size, otherwise uses that

    # LEGEND SETTINGS
    mpl.rcParams['legend.frameon'] = False              # Frames off by default! Better than alpha=0 for illustrator
    mpl.rcParams['legend.framealpha'] = 1               # If I explicitly turn frames on then make them full alpha

    # SAVE SETTINGS
    mpl.rcParams['savefig.transparent'] = True          # Easier manipulation inside of illustrator or presentations
    mpl.rcParams['savefig.bbox'] = 'tight'              # Doesn't cut off the supfigs or side-located legends
    mpl.rcParams['savefig.orientation'] = 'portrait'    # Better for default orientation in publication (the default)

    return None


def kinematic_var_lookup(variable, field):
    """Gets a field value for consistent easy plotting of kinematic variables

    Parameters
    ----------
    variable : str
        Kinematic variable to consider

    field : str
        Field to return.

    Returns
    -------
    field_value : str
        Field value for the specified kinematic variable.

    """
    if variable == 'yaw_vel':
        label = 'Yaw Velocity'
        label_tex = r'$v_{yaw}$'
        color = 'C0'
    elif variable == 'fwd_vel':
        label = 'Forward Vel.'
        label_tex = r'$v_{fwd}$'
        color = 'C1'
    elif variable == 'lat_vel':
        label = 'Lateral Vel'
        label_tex = r'$v_{lat}$'
        color = 'C2'
    elif variable == 'speed':
        label = 'Speed'
        label_tex = r'$\Sigma \vert v \vert$'
        color = 'C3'
    elif variable == 'yaw_lat_vel':
        label = 'Yaw+Lateral Vel'
        label_tex = r'$v_{yaw}+v_{lat}$'
        color = 'C4'
    elif variable == 'fwd_vel_turn_speed':
        label = 'Fwd+Turn Speed'
        label_tex = r'$v_{fwd}+|v_{yaw}+v_{lat}|$'
        color = 'C5'
    elif variable == 'turn_speed':
        label = 'Turn Speed'
        label_tex = r'$|v_{yaw}+v_{lat}|$'
        color = 'C6'
    else:
        raise NotImplementedError("Variable not accounted for")

    if field == 'label':
        return label
    elif field == 'label_tex':
        return label_tex
    elif field == 'color':
        return color
    else:
        raise NotImplementedError("Field not accounted for")


def filter_var_lookup(filter_type, field):
    """Gets a field value for consistent filter calculation and easy plotting of kinematic variables

    Parameters
    ----------
    filter_type : str
        Filter type to consider: "ball_to_vm", "vm_to_ball", "ball_to_fr", "fr_to_ball"

    field : str
        Field to return: "sig_in", "sig_out", "explanation", "title", "units"

    Returns
    -------
    field_value : str
        Field value for the specified kinematic variable.

    NOTE: there is NO explicit error checking here!
    """

    d = {'ball_to_fr':
             {'explanation': "Neuron FR Preceeds Ball <- 0 -> Neuron FR Lags Ball",
              'title': "Ball Activity to Firing Rate",
              'title_tex': r"$Motor \rightarrow Neuron_{spikes}$",
              'units_tex': r"$spikes\ /\ deg$",
              'units': "Neuron FR (Hz) / Ball Units",
              'in': 'kv',
              'out': 'fr'},
         'ball_to_vm':
             {'explanation': "Neuron Vm Preceeds Ball <- 0 -> Neuron Vm Lags Ball",
              'title': "Ball Activity to Membrane Voltage",
              'title_tex': r"$Motor \rightarrow Neuron_{voltage}$",
              'units_tex': r"$V_{mem}mV\ /\ deg$",
              'units': "Neuron Vm (mV) / Ball Units",
              'in': 'kv',
              'out': 'vm'},
         'fr_to_ball':
             {'explanation': 'Neuron FR Lags Ball <- 0 -> Neuron FR Preceeds Ball',
              'title': "Firing Rate to Ball Activity",
              'title_tex': r"$Neuron_{spikes} \rightarrow Motor$",
              'units_tex': r"$degs\ /\ spike$",
              'units': "Ball Units / Neuron FR (Hz)",
              'in': 'fr',
              'out': 'kv'},
         'vm_to_ball':
             {'explanation': "Neuron Vm Lags Ball <- 0 -> Neuron Vm Preceeds Ball",
              'title': "Membrane Voltage to Ball Activity",
              'title_tex': r"$Neuron_{voltage} \rightarrow Motor$",
              'units_tex': r"$degs\ /\ V_{mem}mV$",
              'units': " Ball Units / Neuron Vm (mV)",
              'in': 'vm',
              'out': 'kv'}
         }
    return d[filter_type][field]


def extract_active_times(variable, window_length, fraction_active=0.5, tolerance=0.2, n_zero_bins=5):
    """Select some time snippets from the given kinematic variable of length `window_length` with a
    fraction of activity `fraction_active` (within tolerance value `tolerance). Useful for plotting.
    Parameters
    ----------

    Returns
    -------
    active_inds : np.array

    """

    # Find the active segments
    activity = proc_utils.get_active(variable, n_zero_bins=n_zero_bins)

    # Find a region of window_length that has `fraction_active`
    activity_dm, activity_lags = filter_tools._create_time_lagged_matrix(activity, window_length)
    activity_frac = np.sum(activity_dm, axis=1) / activity_lags.shape[0]

    # Within `tolerance` of the fraction_active
    activity_close_bool = np.isclose(activity_frac, fraction_active, atol=tolerance)
    activity_close_inds = np.where(activity_close_bool == True)[0]

    # from (start,stop) at each window
    useful_ind_pairs = np.array([(i, i + activity_lags.shape[0]) for i in activity_close_inds])

    return useful_ind_pairs
