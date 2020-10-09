"""
Tools to detect and quantify spiking activity in neural recordings
"""

import numpy as np
from scipy.stats import norm
from scipy.signal import convolve, detrend, medfilt, fftconvolve
from scipy.signal import buttord, butter, filtfilt, find_peaks, resample
from scipy.interpolate import interp1d

from tqdm import tqdm

from . stats_tools import mad_standardize, mad


def process_voltage_spike_detection(v, fs):
    """Process voltage timeseries to more easily extract spikes, detrends then uses a ~0.5s rolling estimate of noise

    Parameters
    ----------
    v : np. array
        Voltage data to process for spike detection

    fs : int
        Sampling rate for `v`

    Returns
    -------
    proc_v : np.array
        Processed voltage data

    """
    # Safely redundant with high pass filtering, but saved for using optionally later
    # Detrend signal at 30 and 120 seconds.
    v_dt = _detrend_physiology(v, fs=fs, timescales=[120, 30])

    # High pass filter
    v_hpf = _high_pass_filter(v_dt, fs=fs)

    # Normalize to a rolling estimate of median absolute deviation
    dm_0 = np.roll(v_hpf, -int(fs / 2)).reshape(1 * fs, -1)
    dm_1 = v_hpf.reshape(1 * fs, -1)
    dm_2 = np.roll(v_hpf, int(fs / 2)).reshape(1 * fs, -1)

    # Rolling median and MAD estimates
    mad_ts_0, med_ts_0 = mad(dm_0), np.median(dm_0, axis=0)
    mad_ts_1, med_ts_1 = mad(dm_1), np.median(dm_1, axis=0)
    mad_ts_2, med_ts_2 = mad(dm_2), np.median(dm_2, axis=0)

    mad_ts = (mad_ts_0 + mad_ts_1 + mad_ts_2) / 3
    mad_ts = resample(mad_ts, mad_ts.size*fs)

    med_ts = (mad_ts_0 + mad_ts_1 + mad_ts_2) / 3
    med_ts = resample(med_ts, med_ts.size*fs)

    # Standardize over time
    v_mad = (v_hpf-med_ts)/mad_ts

    return v_mad


def detect_spikes(v, method='findpeaks', prominence=6, window_len=10):
    """Detect spikes in a voltage trace.

    Parameters
    ----------
    v : np.array
        Array to detect spikes within.

    window_len : int, optional
        Window about which to compute prominence

    prominence : float, optional
        Threshold for event detection.

    Returns
    -------
    spikes : np.array (bool)
        Logical array with spikes as True.

    """

    if method == 'thresh':
        # Simple threshold
        v_thresh = v > prominence

        # Find threshold crossings
        starts = np.argwhere(np.diff(v_thresh.astype(int)) > 0).flatten()
        stops = np.argwhere(np.diff(v_thresh.astype(int)) < 0).flatten()

        # add a "stop" at the end if needed
        if starts[-1] > stops[-1]:
            stops = np.append(stops, v.shape[0])

        # very simply define peaks as center points
        peaks = starts + ((stops - starts) / 2).astype(np.int)

    elif method == 'findpeaks':
        peaks, _ = find_peaks(v, height=2.5, wlen=window_len, prominence=prominence)

    # Create spikes array
    spikes = np.zeros_like(v).astype(bool)
    spikes[peaks] = True

    return spikes


def bin_spikes(spikes, bin_len, fs):
    """Bin a timeseries of spiking activity

    Parameters
    ----------
    spikes : np.array (bool)
        Logical array with spike times set to True

    bin_len : float
        Time in seconds of bin length

    fs : float
        Sampling rate for `spikes`.

    Returns
    -------
    bins : np.array
        Binned spikes according to parameters.

    times : np.array
        Timepoint of each bin in a timeseries in seconds.
        Each timepoint represents the _beginning_ of the corresponding bin.

    """

    # Samples per bin
    n_bin_samples = int(bin_len * fs)

    # Fold spikes into bins and sum
    spikes_rs = spikes.copy().reshape(-1, n_bin_samples)
    bins = np.sum(spikes_rs, axis=1)

    # Make a timeseries as well (beginning of each bin)
    times = np.cumsum(np.ones_like(bins) * (n_bin_samples / fs))
    times -= 1 * (n_bin_samples / fs)

    return bins, times


def estimate_fr(binned_spikes, bin_times, sigma_s=0.1):
    """Estimate firing rate from uniformly binned spikes.

    Parameters
    ----------
    binned_spikes : np.array
        Array of spikes binned according to `bin_times`.

    bin_times : np.array
        Array of timepoints corresponding to each entry of `binned_spikes` in seconds.

    sigma_s : float
        Sigma value for gaussian kernel smoothing.

    Returns
    -------
    fr : np.array
        Estimated firing rate

    NOTE: this is after (almost the same as because it is so simple)
    https://github.com/baccuslab/pyret/blob/master/pyret/spiketools.py
    And in retrospect I should have probably used that library more
    """

    # assume all bins are the same length (it says so in the docstring!)
    bin_len = np.diff(bin_times[:2])

    # Make a normalized gaussian kernel
    tau = np.arange(-5 * sigma_s, 5 * sigma_s, bin_len)
    filt = np.exp(-0.5 * (tau / sigma_s) ** 2)
    filt = filt / np.sum(filt)
    size = int(np.round(filt.size / 2))

    # Convolve (non-fft) the binned spikes with the kernel
    fr = convolve(filt, binned_spikes, mode='full', method='direct')

    # Extract valid inds manually (paranoid)
    return fr[size:size + bin_times.size] / bin_len


def upsample_binned_fr(bins, bin_times, upsamp_times):
    """Upsamples and smooths bins from a firing rate estimate using cubic interpolation and rectification

    Parameters
    ----------
    bins : np.array
        Firing rate estimate from bins

    bin_times : np.array
        Array of timepoints corresponding to each entry of `binned_spikes` in seconds.

    upsamp_times : float
        Time points desired (e.g. original sampling times) in seconds.

    Returns
    -------
    resampled : np.array
        Upsampled and smoothed bins at `upsamp_rate`

    """

    # Cubic interpolation
    resampled = interp1d(bin_times,bins,kind='cubic',fill_value='extrapolate')(upsamp_times)

    # Make sure it is always >=0
    resampled[resampled<0] = 0

    # Smooth introduced discontinuity very lightly
    return convolve(resampled, np.array([1,1,1])/3, mode='same', method='direct')


def extract_membrane_voltage(trace, fs, kernel_len=30, detrend_trace=True):
    """Remove spikes from a physiology trace to look at "just membrane voltage". This is a VERY simple approximation
    that uses median filtering. The `kernel_len` needs to be approximately 1.5 times the width of the rising phase
    of the action potential. After median filtering a much smaller window is gaussian filtered to remove hard edges

    Parameters
    ----------
    trace : np.array
        Voltage trace with spiking activity to be removed

    fs : float, int
        Sampling rate for trace

    kernel_len : float
        Length in miliseconds for the median filtering operation. Around 30ms is typically fine

    detrend_trace : bool (optional)
        Flag to hierarchically detrend the data. Useful for calculating filters etc., or remove other drift
        Uses the defaults from `spike_tools._detrend_physiology()`

    Returns
    -------
    trace_vm : np.array
        Membrane voltage estimate after removing spikes

    TODO: convert to views and preallocate all arrays for speed
    """
    # Convert `kernel_len` to samples and make sure it is odd
    ks = int((kernel_len / 1000) * fs)
    ks = ks + 1 - (ks % 2)

    if ks < 13:
        raise Exception('The kernel length is likely too short to allow smoothing after median filtering.')

    # offset to keep for later, calculated as a median the second thirty seconds of the recording for robustness
    offset = np.median(trace[30*fs:60*fs])  # TODO: use a better heuristic
    trace = trace.copy() - offset

    # Median filtering takes a while, but works very well
    chunk_time_s = 1  # TODO: make this adaptive and not reliant on 1s integers
    chunk_size = fs * chunk_time_s
    trace_rs = trace.reshape(-1, chunk_size)

    trace_mf = np.zeros_like(trace_rs, dtype=np.float)
    for i in tqdm(range(trace_rs.shape[0]), desc="Median Filtering Chunked Voltage Traces"):
        trace_mf[i,:] = medfilt(trace_rs[i], ks)

    # stitch chunks back together
    trace_mf_st = np.array(trace_mf[0], dtype=np.float)

    for i in tqdm(range(len(trace_mf) - 1), desc="Stitching Chunked Voltage Traces"):
        # median filter the gaps
        stitch = np.concatenate((trace_mf_st[-ks:], trace_mf[i + 1][:ks]))
        stitch = medfilt(stitch, ks)

        # stitch them back together
        trace_mf_st = np.concatenate((trace_mf_st[:-ks], stitch, trace_mf[i + 1][ks:]))

    # gaussian kernel to clean up further
    gk = np.round(ks/4)
    gaus_kern = norm.pdf(np.linspace(-3.5, 3.5, gk))
    gaus_kern = gaus_kern / np.sum(gaus_kern)

    trace_smoothed = fftconvolve(trace_mf_st, gaus_kern, mode='same')

    if detrend_trace:
        trace_vm = _detrend_physiology(trace_smoothed, fs=fs)
    else:
        trace_vm = trace_smoothed

    return trace_vm + offset


def _high_pass_filter(v, fs, verbose=False):
    """High pass filter for intracellular recording spike detection

    Parameters
    ----------
    v : np.array
        The voltage trace to be filtered

    fs : int
        Sampling rate for `v` and the filter

    verbose : bool
        Print filter info.

    Returns
    -------
    v_filt : np.array
        High pass filtered array from `v`

    """

    # design a simple butterworth filter to remove low frequencies but preserve spike timing
    # it is usually possible to use a simple threshold after this, but `scipy.signal.find_peaks` works also
    order, wn = buttord(wp=100, ws=.05, gpass=3, gstop=30, fs=fs)
    b, a = butter(order, wn, btype='high', fs=fs)

    if verbose:
        # Stability determined by roots
        print(order, wn)
        print(np.abs(np.roots(a)), np.abs(np.roots(b)))

    # zero-phase filter the signal
    return filtfilt(b, a, v)


def _detrend_physiology(v, fs, timescales=None):
    """

    Parameters
    ----------
    v : np.array
        The voltage trace to be detrended.

    fs : int
        Sampling rate for `v` and detrend.

    timescales : float or array-like
        The timescales in seconds to apply detrending over, starts with the longest and goes progressively shorter
        The default value of `None` results in a 120s followed by 30s i.e. timescales=[120,30]
    Returns
    -------
    v_dt : np.array
        Detrended array from `v`

    """

    if timescales is None:
        timescales = [120,30]

    v_dt = np.zeros_like(v)

    # iteratively apply
    for ts in timescales:
        bp = np.arange(0, v.shape[0]-1, ts*fs)
        v_dt += detrend(v, bp=bp) * (1/len(timescales))

    return v_dt
