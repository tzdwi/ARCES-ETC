
from arcesetc import plot_order_sn, plot_order_counts
import astropy.units as u
import matplotlib.pyplot as plt

sptype = 'G4V'
wavelength = 6562
exp_time = 5 * u.min
V = 10 

plot_order_counts(sptype, wavelength, exp_time, V)
plt.show()
