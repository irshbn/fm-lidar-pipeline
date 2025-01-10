import numpy as np
from scipy.signal import hilbert


def hpf(data, dt, fc):
    """Performs high-pass filtering of input signal.

    Args:
        data:   1D signal data vector with uniform time step.
        dt:     Data vector time step (sec).
        fc:     Cutoff frequency (Hz).

    Returns:
        Filtered data vector.
    """

    filtered_data = np.zeros_like(data)
    filtered_data[0] = data[0]

    tau = (2*np.pi*fc)**-1
    alpha = tau/(tau+dt)

    for i in range(1, len(data)):
        filtered_data[i] = alpha * (filtered_data[i-1] + data[i] - data[i-1])

    return filtered_data


def lpf(data, dt, fc):
    """Performs low-pass filtering of input signal.

    Args:
        data:   1D signal data vector with uniform time step.
        dt:     Data vector time step (sec).
        fc:     Cutoff frequency (Hz).

    Returns:
        Filtered data vector.
    """

    smoothed_data = np.zeros_like(data)

    tau = (2*np.pi*fc)**-1
    alpha = dt / (dt + tau)

    smoothed_data[0] = alpha * data[0]

    for i in range(1, len(data)):
        smoothed_data[i] = alpha * data[i] + (1 - alpha) * smoothed_data[i-1]

    return smoothed_data


def kzf(data, m: int, k=4):
    """Performs Kolmogorov-Zurbenko (multiple moving average)
    filtering of input signal.

    Args:
        data:               1D signal data vector.
        m (int):            Moving average window size. Should be odd.
        k (int, optional):  Count of moving average passes. Defaults to 4.

    Returns:
        Filtered data vector.
    """

    if m % 2 == 0:
        return None

    l = (m - 1)//2

    temp_old = data
    temp_new = np.zeros_like(data)

    for i in range(k):
        inds = range(-l, l)
        temp_new[inds] = temp_old[inds]

        for i in range(l, len(data)-l):
            temp_new[i] = np.mean(temp_old[i-l:i+l])

        temp_old = temp_new

    return temp_new


def crop(data, prop):
    """Crops out a central portion of input signal.

    Args:
        data: 1D signal data vector.
        prop: Cropping proportion (0 < prop < 1).

    Returns:
        Central slice of data vector.
    """

    l_index = int(np.trunc(len(data)*(1 - prop) / 2))
    r_index = len(data) - l_index

    return data[l_index:r_index]


def equalize(data):
    """Equalizes input signal by dividing it by its envlope

    Args:
        data: 1D signal data vector.

    Returns:
        Equalized data vector.
    """

    envelope = np.abs(hilbert(data))
    return data / envelope


def precondition(data, dt, fc, prop, wsize: int, repcount=4):
    """Performs filtering, smoothing, cropping and equalization of signal data.

    Args:
        data: 1D signal data vector.
        dt: Data vector time step (sec).
        fc: Highpass filter cutoff frequency (Hz).
        prop: Cropping proportion.
        wsize (int): Moving average window size. Should be odd.
        repcount (int, optional): Count of moving average passes. Defaults to 4.

    Returns:
        Preconditioned data vector.
    """

    filtered_data = hpf(data, dt, fc)
    smoothed_data = kzf(filtered_data, wsize, repcount)
    cropped_data = crop(smoothed_data, prop)
    equalized_data = equalize(cropped_data)

    return equalized_data
