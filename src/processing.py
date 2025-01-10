import numpy as np
from scipy.optimize import curve_fit
from scipy.fft import rfft, rfftfreq, next_fast_len


def fft_freq(data):
    """Estimates the fundamental frequency and phase 
    of input signal through FFT peak search.

    Args:
        data: 1D signal data vector.

    Returns:
        Fundamental frequency and phase estimates.
    """

    n = next_fast_len(len(data))
    fdata = rfft(data, n)
    ffreq = rfftfreq(n)
    max_index = np.argmax(np.abs(fdata))

    return ffreq[max_index], np.angle(fdata)[max_index]


def cosine_curve(x, *params):
    """Generic sine helper function.

    Args:
        x:      X-axis point.
        params: Tuple containing sine frequency and phase.  

    Returns:
        Value of sin(freq*x + phase).
    """

    freq, phase = params

    return np.cos(2*np.pi*freq*x + phase)


def lsq_cosine_fit(data):
    """Performs least-squares fit of a sine curve onto input data

    Args:
        data: 1D signal data vector.

    Returns:
        Tuple, containing frequency estimate, its standard deviation and adjusted R2 coefficient
    """

    n = len(data)
    x = np.arange(n)
    p0 = fft_freq(data)

    popt, pcov = curve_fit(cosine_curve, x, data, p0, method='lm')
    perr = np.sqrt(np.diag(pcov))

    r2 = pcov[0, 1] * pcov[1, 0] / (pcov[0, 0] * pcov[1, 1])
    p = 2
    adj_r2 = 1 - (1-r2)*(n-1)/(n-p-1)

    return popt[0], perr[0], adj_r2
