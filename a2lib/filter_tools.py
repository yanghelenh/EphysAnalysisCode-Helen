"""
Tools and utility functions for computing linear and nonlinear filters between neural data and stimuli.

Also contains functions for predicting linear or linear-nonlinear system outputs.
"""

import numpy as np
from scipy.fftpack import fft, ifft, fftshift, fftfreq
from scipy.fftpack.helper import next_fast_len
from scipy.signal import fftconvolve, convolve
from scipy.signal.windows import dpss
from scipy.stats import linregress


def get_linear_filter(signal_in, signal_out, window_length, sample_rate, center_input=True, freq_cut=None):
    """Return linear filter and lag time values indices.

    Note this function only removes `signal_in` autocorrelation and optionally a hard low-pass at `freq_cut`.
    Any other filtering is meant to be applied afterwards for flexibility.

    Parameters
    ----------
    signal_in : np.array
        The "stimulus" or system input. `signal_in` has its autocorrelation normalized in the resulting filter.
        Must share a sampling rate `sample_rate` and shape[0] time base with signal_out

    signal_out : np.array
        The "response" or system output. `signal_out` is cross correlated with `signal_in`.
        Must share a sampling rate `sample_rate` and shape[0] time base with signal_in

    window_length : float
        number of seconds before and after a point, resulting lags will be of size (sample_rate * window_length * 2)

    sample_rate : int
        rate for signal_in and signal_out

    center_input : bool (optional)
        Mean subtract (center) the windows where the FFT is being applied. Defaults to `True`

    freq_cut : int (optional)
        hard cutoff frequency applied in the frequency domain. If `None` there will be no cut applied.

    Returns
    ------
    linear_filt : np.array
        of shape (window_length * 2)

    lags : np.array
        in relative samples to center (zero) lag point

   """

    # Generate time-lagged matrices to simplify math and vectorization below
    window_length_samples = np.round(window_length * sample_rate).astype(int)
    signal_in_dm, lags_dm = _create_time_lagged_matrix(signal_in.copy(), window_length=window_length_samples)
    signal_out_dm, _ = _create_time_lagged_matrix(signal_out.copy(), window_length=window_length_samples)

    # Find the next fast length to zero pad fft calculations to
    fft_len = next_fast_len(lags_dm.shape[0])
    fft_kwargs = {'n': fft_len, 'axis': -1, 'overwrite_x': True}

    # Center the signal input if specified
    if center_input:
        signal_in_ffts = fft(signal_in_dm - np.expand_dims(np.mean(signal_in_dm, axis=1), axis=-1), **fft_kwargs)
        signal_out_ffts = fft(signal_out_dm - np.expand_dims(np.mean(signal_out_dm, axis=1), axis=-1), **fft_kwargs)
    else:
        signal_in_ffts = fft(signal_in_dm, **fft_kwargs)
        signal_out_ffts = fft(signal_out_dm, **fft_kwargs)

    # Compute the average input autocorrelation and the average input output crosscorrelation
    sig_in_autocorr_fft = np.mean(signal_in_ffts * np.conj(signal_in_ffts), axis=0)
    in_out_crosscorr_fft = np.mean(signal_out_ffts * np.conj(signal_in_ffts), axis=0)

    # remove high frequencies with an appropriate window
    fft_freqs = fftfreq(sig_in_autocorr_fft.shape[0], d=1 / sample_rate)
    window_fft = np.ones_like(fft_freqs)  # Default to no window

    if freq_cut is not None:
        window_fft *= np.finfo(float).eps
        window_fft[np.abs(fft_freqs) < freq_cut] = 1.0

    filter_fft = window_fft * in_out_crosscorr_fft / sig_in_autocorr_fft
    linear_filt = fftshift(np.real(ifft(filter_fft)))

    # empirically determined that this needs to be flipped to agree with lags_dm variable
    linear_filt = np.flip(linear_filt)

    return linear_filt, lags_dm/sample_rate


def lpf_linear_filter(linear_filter, bandwidth, fs, kind='dpss'):
    """Low pass filter a linear filter via convolution of with a window (e.g. DPSS window)

    Parameters
    ----------
    linear_filter : np.array
        Linear filter to be convolved

    bandwidth : float
        The preserved bandwidth desired in Hz

    fs : float
       The sampling rate of the linear filter

    kind : np.array (optional)
        The kind of window to apply convolution. By default uses a dpss window.

    Returns
    -------
    lpf_filter : np.array
        Window convolved linear filter
    """

    # Currently only have the dpss window set up
    if kind == 'dpss':
        window = _make_dpss_window(bandwidth_hz=bandwidth, fs=fs)
    else:
        raise NotImplementedError('window kind not yet implemented.')

    return fftconvolve(linear_filter, window, 'same')


def get_linear_filter_prediction(filter_, lags, input_, filter_length=None, flip_filter=True):
    """Get the predicted output of a filter `filter_` given an input `input_` by convolution

    Parameters
    ----------
    filter_ : np.array
        The linear filter to apply to `input_`.

    lags : np.array
        Timebase for linear filter in seconds.

    input_ : np.array
        Convolved by `filter_` to produce `output`. Must be at the same sampling rate as filter_

    filter_length : float
        Length of the filter to actually apply. Filter will extend `filter_length / 2` seconds before and after the zero
        lag point in `lags`.

    flip_filter :  bool
        Indicates if the filter needs to be flipped before convolution.

    Returns
    -------
    output : np.array
        The output of a convolution between `filter_` and `input_`. With the same timebase as input_

    times : np.array
        The time points corresponding to the output. Returned in case the convolution method changes
    """
    if flip_filter:
        filter_ = np.flip(filter_.copy())

    # If the filter length requested is the max just treat it this way
    if np.isclose(lags[-1], filter_length):
        filter_length = None

    if filter_length is not None:
        # This is simply symmetrical for now.
        inds = np.abs(lags) < (filter_length / 2)
        # Add back in the later edge to keep the applied filter an even number of samples
        edge_ind = np.argwhere(np.diff(inds))[0][-1]
        inds[edge_ind] = True
        filter_trunc = filter_.copy()[inds]
    else:
        filter_trunc = filter_.copy()

    # Convolve input and (modified) filter
    output = convolve(input_, filter_trunc, mode='full')[:input_.shape[0]]

    # Construct timebase
    isi = np.mean(np.diff(lags))  # inter sample interval
    times = np.cumsum(isi * np.ones_like(input_)) - isi

    return output, times


def variance_explained(prediction, observation):
    """Calculate variance explained by a linear regression model
    Parameters
    ----------
    prediction : np.array
        Predicted output of linear filter

    observation : np.array
        Actual system output measured.

    Returns
    -------
    var_explained : float
        Explained variance (R^2 value)
    """

    # X = prediction.reshape(-1,1)
    # y = observation.reshape(-1,1)
    # lr = LinearRegression(normalize=False).fit(X=X, y=y)
    # return lr.score(X=X, y=y)

    # Use straightforward linear regression
    lr = linregress(prediction,observation)
    return lr.rvalue


def get_naive_nonlinearity(linear_prediction, observation, n_bins='fd'):
    """Finds the nonlinearity in a Linear-Nonlinear model with filter output `linear_prediction` by a simple histogram
    lookup. i.e. the classic output nonlinearity

    Parameters
    ----------
    linear_prediction : np.array
        Signal resulting from convolution of an input with the linear filter.

    observation : np.array
        Actual output of the system. Used to determine the lookup function returned as `nonlinear_filter`.

    n_bins : int (optional)
        Number of bins over which to compute the nonlinearity lookup table. `fd` (default) it uses the numpy efficient
        estimate for the optimal number of histogram bins for this dataset size.

    Returns
    -------
    nonlinear_filter : np.array
        Nonlinear filter that makes up for the difference between linear filter output and system output. Contains Nans

    bin_edges : np.array
        values corresponding to bins edges in the `nonlinear_filter`. Edge conventions are from `np.histogram()`
    """

    # Determine bin edge and center locations
    if n_bins == 'fd':
        bin_edges = np.histogram_bin_edges(np.concatenate((observation, linear_prediction)), bins='fd')
    else:
        bin_edges = np.histogram_bin_edges(np.concatenate((observation, linear_prediction)), bins=n_bins)

    obs_hist = np.histogram(observation, bins=bin_edges)[0]
    lin_pred_hist = np.histogram(linear_prediction, bins=bin_edges)[0]

    nl_hist = obs_hist / lin_pred_hist

    # Remove infs
    nonlinear_filter = nl_hist.copy()
    nonlinear_filter[np.isinf(nl_hist)] = np.nan

    return nonlinear_filter, bin_edges


def apply_nonlinearity(input_, bin_values, bin_edges):
    """Convenience function for getting nonlinear output from histogram.

    Parameters
    ----------
    input_ : np.array
       The values to apply nonlinear lookup to. i.e. the ouput of a linear filter step (aka the linear prediction)

    bin_values : np.array
        Nonlinear filter to apply.

    bin_edges : np.array
        Edges correspoinding to `bin_values`. Edge conventions from `np.histogram()`

    Returns
    -------
    nl_pred : np.array
        Predicted output after

    """

    # numpy still has not fixed digitize to agree w/histogram edges
    def _fix(x):
        if x == 0:
            return 1
        else:
            return x
    input_digital = np.digitize(input_, bins=bin_edges[:-2], right=True)
    input_digital = np.array([_fix(x) for x in input_digital])

    # apply lookup table
    nl_pred = bin_values[input_digital]
    nl_pred[np.isnan(nl_pred) | np.isinf(nl_pred)] = 0
    return nl_pred


def truncate_filter(filter_, filter_lags, trunc_len_s):
    """
    Parameters
    ----------
    filter_ : np.array
        Linear filter to be truncated down to `trunc_len_s`

    filter_lags : np.array
        Monotonically increasing and signed timebase for `filter`. In seconds.

    trunc_len_s : float
        Full truncated filter length, half allocated to either side of zero.

    Returns
    -------
    filt_trunc : np.array
        Truncated linear filter.

    lags_trunc : np.array
        Timebase corresponding to the truncated filter `filt_trunc` in seconds.
    """
    # Return the linear filter truncated to an even number of samples close to filter_len_s
    inds = np.abs(filter_lags) < (trunc_len_s / 2)

    # Add back in the later edge to keep the applied filter an even number of samples
    if np.sum(inds) % 2 == 1:
        edge_ind = np.argwhere(np.diff(inds))[0][-1]
        inds[edge_ind] = True
    return filter_.copy()[inds], filter_lags.copy()[inds]


def _make_dpss_window(bandwidth_hz, fs):
    """Get the first Slepian (DPSS) window with a bandwidth `bandwidth_hz`.

    Parameters
    ----------
    bandwidth_hz : float, int
        Bandwidth of the window in Hz.

    fs : float, int
        Sampling rate for the window (and signal to be filtered).

    Returns
    -------
    window : np.array
        The first DPSS sequence (window) with the highest possible energy preservation in `bandwidth_hz`

    Notes
    -----
    To get a single sequence set NW = 1 (since there are 2*NW - 1 sequences with the criteria).
    The matlab documentation is somewhat helpful here: https://www.mathworks.com/help/signal/ref/dpss.html
    """

    # win_len/fs * BW = 2 =>  win_len = BW*fs/2
    window_length = int(1/bandwidth_hz*fs*0.5)
    return dpss(window_length, NW=1, Kmax=None)


def _create_time_lagged_matrix(input_, window_length):
    """Generate temporally overlapping segments of an input (e.g. stimulus).

    Parameters
    ----------
    input_ : np.array
        Note: avoids locally overloading python `input`

    window_length : int
        number of samples before and after a point, resulting lags will be of size (window_length * 2)

    Returns
    ------
    design_mat : np.array
        of shape (stimulus.shape[0] - lags.shape[0], window_length * 2)

    lags : np.array
        in relative samples to center (zero) lag point

    NOTE: this is probably better off as a numpy view, but allocation is not computationally limiting at the moment.
    """
    # lags.shape[0] == (window_length * 2)
    lags = np.arange(-window_length, window_length)

    # pad the stimulus vector on each side by window length
    padded_vec = np.concatenate((np.zeros(window_length), input_.copy(), np.zeros(window_length)))

    # time lagged and advanced versions of stimulus according to `lags` vector
    # the transpose is to generate typical design matrix format (samples x features)
    # i.e. padded_mat.shape == padded_vec.shape[0],lags.shape[0]
    padded_mat = np.stack([np.roll(padded_vec, shift=l) for l in lags]).T

    # remove a lags.shape[0] segment from the start and end. Resulting in:
    # design_mat.shape[0] == stimulus.shape[0] - lags.shape[0] + 1
    design_mat = padded_mat[lags.shape[0] - 1:-lags.shape[0] + 1, :]

    return design_mat, lags
