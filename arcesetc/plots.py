import h5py 
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from json import load
import os

__all__ = ['plot_order_counts', 'plot_order_sn', 'available_sptypes']

directory = os.path.dirname(__file__)
sptypes = load(open(os.path.join(directory, 'data', 'sptype_dict.json'), 'r'))
archive = h5py.File(os.path.join(directory, 'data', 'archive.hdf5'), 'r+')

sptype_to_temp = load(open(os.path.join(directory, 'data',
                                        'sptype_to_temp.json'), 'r'))
spectral_types = [key for key in sptype_to_temp.keys() if key in sptypes]
temps = np.array([sptype_to_temp[key] for key in spectral_types])


def closest_sptype(sptype):
    return spectral_types[np.argmin(np.abs(sptype_to_temp[sptype] - temps))]


def closest_target(sptype):
    closest_spectral_type = closest_sptype(sptype)
    return sptypes[closest_spectral_type]


def available_sptypes():
    return spectral_types


@u.quantity_input(exp_time=u.s, wavelength=u.Angstrom)
def plot_order_counts(sptype, wavelength, exp_time, V):

    target = closest_target(sptype)

    matrix = archive[target][:]
    template_vmag = archive[target].attrs['V'][0]

    closest_order = np.argmin(np.abs(matrix[:, 0] -
                                     wavelength.to(u.Angstrom).value))
    lam_0, delta_lam, n_lam = matrix[closest_order][:3]
    polynomial_coeffs = matrix[closest_order][3:]
    wave = np.arange(lam_0 - n_lam*delta_lam/2, lam_0 + n_lam*delta_lam/2,
                     delta_lam)
    magnitude_scaling = 10**(0.4 * (template_vmag - V))
    flux = (np.polyval(polynomial_coeffs, wave-lam_0) * exp_time.to(u.s).value *
            magnitude_scaling)

    fig, ax = plt.subplots()
    ax.plot(wave, flux)
    ax.set_xlabel('Wavelength [Angstrom]')
    ax.set_ylabel('Flux [DN]')
    for s in ['right', 'top']:
        ax.spines[s].set_visible(False)
    ax.grid(ls=':', color='silver')
    return fig, ax


@u.quantity_input(exp_time=u.s, wavelength=u.Angstrom)
def plot_order_sn(sptype, wavelength, exp_time, V):

    target = closest_target(sptype)

    matrix = archive[target][:]
    template_vmag = archive[target].attrs['V'][0]

    closest_order = np.argmin(np.abs(matrix[:, 0] -
                                     wavelength.to(u.Angstrom).value))
    lam_0, delta_lam, n_lam = matrix[closest_order][:3]
    polynomial_coeffs = matrix[closest_order][3:]
    wave = np.arange(lam_0 - n_lam*delta_lam/2, lam_0 + n_lam*delta_lam/2,
                     delta_lam)
    magnitude_scaling = 10**(0.4 * (template_vmag - V))
    flux = (np.polyval(polynomial_coeffs, wave-lam_0) * exp_time.to(u.s).value *
            magnitude_scaling)

    sn = flux / np.sqrt(flux)

    fig, ax = plt.subplots()
    ax.plot(wave, sn)
    ax.set_xlabel('Wavelength [Angstrom]')
    ax.set_ylabel('Signal/Noise')
    for s in ['right', 'top']:
        ax.spines[s].set_visible(False)
    ax.grid(ls=':', color='silver')
    return fig, ax