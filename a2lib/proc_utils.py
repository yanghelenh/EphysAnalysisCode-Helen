"""
Specialized helper functions for idiosyncratic processing
"""
import os, glob
import numpy as np
import scipy as sp

from scipy import io, signal, stats
from scipy import interpolate, convolve

def convert_ball_units(data,units):
    # conversion factor for the ball used is 0.0553 from deg / s to mm / s

    return data


def _fd_optimal_inact_bins(var, n_zero=1, use_mode=True):
    # if using mode, then n_zero is the number to zero out
    # if not, it takes a symmetrical +1 number around the central bin

    # fd method is very robust!
    var_counts, var_bin_edges = np.histogram(var, bins='fd')
    var_bin_cents = var_bin_edges + np.diff(var_bin_edges)[0] / 2
    if use_mode:
        # get bin that is closest to the mode
        modal_bins = np.argsort(var_counts)[::-1][:n_zero]
        sort_modal_bin_edges = np.sort(var_bin_edges[modal_bins])
        lower_bound = sort_modal_bin_edges[0]
        upper_bound = sort_modal_bin_edges[-1]
    else:
        # or get the bin that is closest to zero
        min_bin = np.argmin(-np.abs(var_bin_cents))
        lower_bound = var_bin_edges[min_bin - n_zero + 1]
        upper_bound = var_bin_edges[min_bin + n_zero]
    return lower_bound, upper_bound


def _get_nonzero(var, lower_bound, upper_bound):
    # Returns logical
    nonzero = np.zeros_like(var).astype(bool)
    nonzero[var < lower_bound] = True
    nonzero[var > upper_bound] = True
    return nonzero


def get_active(var, n_zero_bins=1):
    # Wrapper
    lb, ub = _fd_optimal_inact_bins(var, n_zero_bins)
    return _get_nonzero(var, lb, ub)


def extract_windows(x, inds, n_before, n_after):
    valid = ((inds - n_before >= 0) & (inds + n_after + 1 < x.size)).flatten()
    inds_valid = inds[valid]
    return np.stack([x[int(ind-n_before):int(ind+n_after+1)] for ind in inds_valid], axis=1)


def getpeaks(mat, thresh, fs, ms_prev_post, ms_max_window, verbose=False):
    # this needs to be cleaned up.. last two steps are redundant now
    if verbose:
        print(thresh)

    # Find threshold crossings
    changes = np.diff(mat > thresh, prepend=0)
    peak_starts = np.argwhere(changes == +1).flatten()
    peak_stops = np.argwhere(changes == -1).flatten()

    # remove a start if it happens at the end
    if peak_starts.size > peak_stops.size:
        peak_starts = peak_starts[:-1]

    # find peak index with argmax
    peak_locs_all = np.empty(0)
    for (curr_start, curr_stop) in zip(peak_starts, peak_stops):
        curr_loc = np.argmax(mat[curr_start:curr_stop + 1])
        peak_locs_all = np.append(peak_locs_all, curr_start + curr_loc).astype(int)
    if verbose:
        print("All: ", peak_locs_all.shape[0])

    # exclude peaks at the start of end of the recording
    samps_prev_post = np.int(ms_prev_post / 1000 * fs)
    peak_locs_clean = peak_locs_all[
        ~ ((peak_locs_all < samps_prev_post) | (peak_locs_all > mat.size - samps_prev_post))]
    if verbose:
        print("Clean 1: ", peak_locs_clean.shape[0])

    # exclude peaks that don't meet additional criteria
    max_window = np.int(ms_max_window / 1000 * fs)
    peak_locs = np.empty(0)
    for curr_loc in peak_locs_clean:
        if (mat.shape[0] > curr_loc + max_window) & (curr_loc - max_window >= 0):
            peak_locs = np.append(peak_locs, curr_loc).astype(int)
    if verbose:
        print("Final: ", peak_locs.shape[0])

    return peak_locs


def safe_interp_conv_smooth_ball_data(vec, pad_len, kern_len):
    """ A somewhat agressive smoothing protocol that makes later upsampling work well
    'Standard Method'
    1 - Pads the ends of the timeseries using windowed interpolation 
    2 - Upsamples by a factor of 10
    3 - Bidirectionally convolves with a 3.5 sigma truncated gaussian kernel 
    4 - Downsamples by a factor of 10 
    5 - Extracts the valid region of smoothed convolution timeseries.
    
    Parameters
    ----------
        vec : np.array

        pad_len : np.int

        kern_len : np.int

    Returns
    -------
        vec_conv : np.array

    """
    # get pads
    pre = interpolated_pad(vec, pad_len, 'start')
    post = interpolated_pad(vec, pad_len, 'end')
    
    # add pads
    vec_pad = np.concatenate([pre, vec, post])
    
    # Hardcoded for now...
    us_fac = 10
     
    # upsample 
    vec_pad_us = signal.resample_poly(vec_pad, us_fac, 1)
    
    # gaussian kernel for [kern_len] samples
    kern_us = stats.norm.pdf(np.linspace(-3.5, 3.5, us_fac*kern_len))
    kern_us = kern_us / np.sum(kern_us)
    
    # convolve in each direction and average
    vec_conv_us_fwd = convolve(vec_pad_us, kern_us, mode='same')
    vec_conv_us_rev = convolve(vec_pad_us[::-1], kern_us[::-1], mode='same')[::-1]
    vec_conv_us = (vec_conv_us_fwd + vec_conv_us_rev)/2
    
    # downsample 
    vec_conv = signal.resample_poly(vec_conv_us, 1, us_fac)
    
    # remove pads before returning
    return vec_conv[pad_len:][:-pad_len]


def interpolated_pad(vec, pad_len, pad_side):
    """ Pad by interpolation. Useful for convolution or pointwise differentiation.
    
    Only returns the padded region, reasonably robust for prepending and appending.
    
    Parameters
    ----------
    vec : np.array - no need to truncate

    pad_len : np.int

    pad_side : str - 'start' or 'end'

    Returns
    -------
    pad : np.array

    """

    if pad_side == 'start':
        pad_f = interpolate.CubicSpline(np.arange(pad_len), vec[:pad_len], extrapolate=True)
        pad = pad_f(np.arange(pad_len)-pad_len)
    elif pad_side == 'end':
        pad_f = interpolate.CubicSpline(np.arange(pad_len), vec[-pad_len:], extrapolate=True)
        pad = pad_f(np.arange(pad_len)+pad_len)
    else:
        raise Exception('pad_side must be either "start" or "end"')
    return pad


def safe_interp_diff_ball_data(vec, pad_len):
    """Differentiate ball data to get acceleration. Note this relies on the ball data already
    being smoothed to a large extent. i.e. `vec` needs to be smoothly varying. 
    
    1 - Pads the ends of the timeseries using windowed interpolation
    2 - Upsamples by a factor of 10
    3 - Takes the pointwise difference
    4 - Downsamples by a factor of 10 
    5 - Extracts the valid region of timeseries.
    
    Parameters
    ----------
    vec : np.array

    pad_len : np.int

    Returns
    -------
    vec_diff : np.array

    """
    # get pads
    pre = interpolated_pad(vec, pad_len, 'start')
    post = interpolated_pad(vec, pad_len, 'end')
    
    # add pads
    vec_pad = np.concatenate([pre, vec, post])
    
    # Hardcoded for now...
    us_fac = 10
    
    # upsample 
    vec_pad_us = signal.resample_poly(vec_pad, us_fac, 1)
    
    # differentiate
    vec_pad_us_diff = np.diff(vec_pad_us, n=1)
    
    # downsample
    vec_pad_diff = signal.resample_poly(vec_pad_us_diff, 1, us_fac)
    
    # remove padding and return
    return vec_pad_diff[pad_len:][:-pad_len]


def get_starts_stops_durs(vec):
    """Get the rising and falling edges and duration of a logical with some error checking
    Parameters
    ----------
    vec : np.array

    Returns
    -------
    starts : np.array

    stops : np.array

    diffs : np.array

    """
    vec_delta = np.diff(vec.astype(np.int), prepend=vec[0].astype(np.int))
    starts = np.argwhere(vec_delta == +1).flatten()
    stops = np.argwhere(vec_delta == -1).flatten()
    # add a stop at the final index if it never stopped
    if (stops[-1]-starts[-1]) <= 0:
        stops = np.append(stops, vec.shape[0])
    if starts.shape[0] < stops.shape[0]:
        starts = np.concatenate(([0], starts))
    diffs = stops - starts
    return starts, stops, diffs


# Automate input and output assignment over multiple data sources
def io_var_to_signals(in_var, out_var, neural_sig_type, neuron_list, ephys_dict, kinematic_var, ball_dict):
    sig = [[],[]]
    for i,var in enumerate([in_var,out_var]):
        if var == 'vm':
            # fr sum and diff are calculated here
            if neural_sig_type == 'sum':
                sig[i] = ephys_dict['vm'][neuron_list[0]] + ephys_dict['vm'][neuron_list[1]]
            elif neural_sig_type == 'diff':
                sig[i] = ephys_dict['vm'][neuron_list[0]] - ephys_dict['vm'][neuron_list[1]]
            else:
                sig[i] = ephys_dict['vm'][neural_sig_type]
        elif var == 'fr':
            # vm sum and diff are calculated here
            if neural_sig_type == 'sum':
                sig[i] = ephys_dict['fr'][neuron_list[0]] + ephys_dict['fr'][neuron_list[1]]
            elif neural_sig_type == 'diff':
                sig[i] = ephys_dict['fr'][neuron_list[0]] - ephys_dict['fr'][neuron_list[1]]
            else:
                sig[i] = ephys_dict['fr'][neural_sig_type]
        elif var == 'kv':
            # Now `speed` and `yaw_lat_vel` calculated ahead of time for consistency across notebooks
            sig[i] = ball_dict[kinematic_var]

    return sig[0], sig[1]


def get_non_stim_inds(ephys, ball, exclude_dur=0.25):
    """Nasty func. needs primitives
    """
    ds_fac = int(ephys['fs']/ball['fs'])
    updown_fac = 4
    stim_ds = sp.signal.upfirdn([1, 1, 1], ephys['stim'],
                                         up=updown_fac, down=updown_fac*ds_fac)[:ball['t_s'].shape[0]]
    stim_on = (stim_ds == 5).astype(np.float)
    stim_off = (stim_ds != 5).astype(np.float)

    if np.sum(stim_on) > 0:

        stim_starts = np.diff(stim_on, prepend=0) == 1
        stim_starts = np.argwhere(stim_starts).flatten()

        stim_stops = np.diff(stim_on, prepend=0) == -1
        stim_stops = np.argwhere(stim_stops).flatten()

        exclude_ind = np.int(exclude_dur * ball['fs'])

        # Just when the stimulus is on
        stim_only = stim_on.copy()

        # When the stimulus is on plus an excluded time
        stim_plus_exclude = np.concatenate(
            [np.arange(start, stop + exclude_ind, 1) for start, stop in zip(stim_starts, stim_stops)])
        stim_exclude = np.zeros_like(stim_only)
        stim_exclude[stim_plus_exclude] = 1

        # Just when the stimulus is off
        no_stim_only = stim_off.copy()

        # When the stimulus is off plus some included time
        no_stim_exclude = np.ones_like(stim_only)
        no_stim_exclude[stim_plus_exclude] = 0

    else:
        stim_only = np.zeros_like(stim_on)
        stim_exclude = np.zeros_like(stim_on)
        no_stim_only = np.ones_like(stim_on)
        no_stim_exclude = np.ones_like(stim_on)

    return stim_only.astype(bool), stim_exclude.astype(bool), no_stim_only.astype(bool), no_stim_exclude.astype(bool)


def get_ephys_ds(ephys, ball, updown_fac=2):
    """ Constructs a dictionary of resampled (downsampled) timeseries
    """

    # Construct a dictionary to hold resampled timeseries
    ephys_ds = {'vm': {}, 'fr': {}}
    ephys_ds['vm']['t_s'] = ball['t_s']
    ephys_ds['fr']['t_s'] = ball['t_s']

    ds_fac = int(ephys['fs']/ball['fs'])
    if 'a2_vm_l' in ephys:
        ephys_ds['vm']['a2_l'] = sp.signal.upfirdn([1, 1, 1], ephys['a2_vm_l'],
                                                   up=updown_fac, down=updown_fac*ds_fac)[:ball['t_s'].shape[0]]
    if 'a2_vm_r' in ephys:
        ephys_ds['vm']['a2_r'] = sp.signal.upfirdn([1, 1, 1], ephys['a2_vm_r'],
                                                   up=updown_fac, down=updown_fac*ds_fac)[:ball['t_s'].shape[0]]
    if 'a2_fr_l' in ephys:
        ephys_ds['fr']['a2_l'] = sp.signal.upfirdn([1, 1, 1], ephys['a2_fr_l'],
                                                   up=updown_fac, down=updown_fac*ds_fac)[:ball['t_s'].shape[0]]
    if 'a2_fr_r' in ephys:
        ephys_ds['fr']['a2_r'] = sp.signal.upfirdn([1, 1, 1], ephys['a2_fr_r'],
                                                   up=updown_fac, down=updown_fac*ds_fac)[:ball['t_s'].shape[0]]
    # fml here, get it done!
    if 'a1_vm_l' in ephys:
        ephys_ds['vm']['a1_l'] = sp.signal.upfirdn([1, 1, 1], ephys['a1_vm_l'],
                                                   up=updown_fac, down=updown_fac * ds_fac)[:ball['t_s'].shape[0]]
    if 'a1_vm_r' in ephys:
        ephys_ds['vm']['a1_r'] = sp.signal.upfirdn([1, 1, 1], ephys['a1_vm_r'],
                                                   up=updown_fac, down=updown_fac * ds_fac)[:ball['t_s'].shape[0]]
    if 'a1_fr_l' in ephys:
        ephys_ds['fr']['a1_l'] = sp.signal.upfirdn([1, 1, 1], ephys['a1_fr_l'],
                                                   up=updown_fac, down=updown_fac * ds_fac)[:ball['t_s'].shape[0]]
    if 'a1_fr_r' in ephys:
        ephys_ds['fr']['a1_r'] = sp.signal.upfirdn([1, 1, 1], ephys['a1_fr_r'],
                                                   up=updown_fac, down=updown_fac * ds_fac)[:ball['t_s'].shape[0]]
    # New field that should be in all flies
    ephys_ds['stim'] = sp.signal.upfirdn([1, 1, 1], ephys['stim'],
                                         up=updown_fac, down=updown_fac*ds_fac)[:ball['t_s'].shape[0]]
    return ephys_ds


class EphysSignalNames(object):
    """Deal with the a1 dict keys, this is kind of a hack, but was the fastest fix I could think of generating new figures...
    Because there are always left cells but not always right cells this works out okay :)
    """

    def __init__(self, neuron_list):
        self.neuron_list = neuron_list

    @property
    def l(self):
        return self.neuron_list[0]

    @property
    def r(self):
        if len(self.neuron_list) < 2:
            return ''
        else:
            return self.neuron_list[1]

    @property
    def vm_l(self):
        return self.neuron_list[0][:2] + '_vm_' + self.neuron_list[0][-1]

    @property
    def vm_r(self):
        if len(self.neuron_list) < 2:
            return ''
        else:
            return self.neuron_list[1][:2] + '_vm_' + self.neuron_list[1][-1]

    @property
    def fr_l(self):
        return self.neuron_list[0][:2] + '_fr_' + self.neuron_list[0][-1]

    @property
    def fr_r(self):
        if len(self.neuron_list) < 2:
            return ''
        else:
            return self.neuron_list[1][:2] + '_fr_' + self.neuron_list[1][-1]

    @property
    def spikes_l(self):
        return self.neuron_list[0][:2] + '_spikes_' + self.neuron_list[0][-1]

    @property
    def spikes_r(self):
        if len(self.neuron_list) < 2:
            return ''
        else:
            return self.neuron_list[1][:2] + '_spikes_' + self.neuron_list[1][-1]


def load_preproc_fly(save_npy_dir, skip_ephys=False):
    """
    Parameters
    ----------
    saved_npy_dir : str

    skip_ephys : bool

    Returns
    -------

    ball : dict

    ephys : dict, empty dict

    ephys_ds : dict

    """
    # Load in the preprocessed npy files
    ball = np.load(glob.glob(save_npy_dir + os.path.sep + 'ball.npy')[0], allow_pickle=True).item()

    ephys_ds = np.load(glob.glob(save_npy_dir + os.path.sep + 'ephys_ds.npy')[0], allow_pickle=True).item()

    # Time consuming, avoid if not needed
    if skip_ephys:
        ephys = {}
    else:
        ephys = np.load(glob.glob(save_npy_dir + os.path.sep + 'ephys.npy')[0], allow_pickle=True).item()

    return ball, ephys, ephys_ds


def load_sr_matfile_data(mat_filepath, neuron_labels):
    """ Import data from the .mat files (unprocessed experimental outcome)

    Parameters
    ----------
    mat_filepath : str
        Full filepath to the `.mat` file

    neuron_labels : list
        The labels that correspond to the mat fields [Ephys_A] or [Ephys_A,Ephys_B]
        i.e. single vs dual recordings

    Returns
    -------
    ball : dict, pd.DataFrame


    ephys : dict, pd.DataFrame


    """
    data = io.loadmat(mat_filepath)

    ball = {'fs': int(data['ball_SR']),
            't_s': data['t_ball'].reshape(-1).astype(np.float),
            'fwd': data['fwd'].reshape(-1).astype(np.float),
            'yaw': data['yaw'].reshape(-1).astype(np.float),
            'lat': data['lat'].reshape(-1).astype(np.float)}

    ephys = {'fs': int(data['ephys_SR']),
             't_s': data['t_ephys'].reshape(-1).astype(np.float),
             'stim':  data['stim'].reshape(-1).astype(np.float)}

    ephys[neuron_labels[0]] = data['ephys_A'].reshape(-1).astype(np.float)
    if len(neuron_labels) > 1:
        ephys[neuron_labels[1]] = data['ephys_B'].reshape(-1).astype(np.float)

    return ball, ephys
