import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.constants as const
import sirtipy
from scipy.integrate import trapz, quad
from scipy.interpolate import interp1d
from scipy.special import gamma as gammafunc


# PART A: SYNCHROTRON

# Figure out the normalization constant C such that
#  n = C int_gmin^gmax g^-p dg
def powerlaw_normalize(p, gmin, gmax, n):
    integral = (gmax*(1.-p) - gmin**(1.-p)) / (1.-p)
    return n / integral


# Synchrotron emission coefficient from a power-law electron distribution.
# Based on RL equation 6.36.
def j_synchrotron(frequency, location, Inu, B, p, gmin, gmax, ne):
    # Parameters: B=magnetic field (Gauss)  p=power law index
    #   gmin,gmax=gamma range for power law   ne=electron number density (cm^-3)
    
    # Figure out normalization C
    C = powerlaw_normalize(p, gmin, gmax, ne)
    # <|sin alpha|> = 2/pi
    sinalpha = 2. / np.pi
    # RL 6.36. Note that in RL it is per unit omega, so we need to multiply by 2pi (and also
    # use nu instead of omega in the formula)
    Ptot = np.sqrt(3) * const.e.esu.value**3 * C * B * sinalpha / (const.m_e.cgs.value * const.c.cgs.value**2 *
        (p+1.)) * gammafunc(0.25*p + 19./12) * gammafunc(0.25*p - 1./12) * (const.m_e.cgs.value * const.c.cgs.value *
        frequency * 2. * np.pi / (3. * const.e.esu.value * B * sinalpha))**(-(p-1.)/2)
    jnu = Ptot / (4. * np.pi)

    return jnu

# Synchrotron absorption coefficient from a power-law electron distribution.
# Based on RL equation 6.53.
# FIXME: Complete the alpha_synchrotron function. You can use the definition of j_synchrotron as
# a guide. Note that I have defined a function called "powerlaw_normalize" that will calculate
# the value of the normalization constant C. Also note that I hvae renamed the Gamma function
# scipy.special.gamma to "gammafunc" so that there is no danger of confusing it with the Lorentz factor.
#
#def alpha_synchrotron(frequency, location, Inu, B, p, gmin, gmax, ne):
    #...


# Create the region
medium_a = sirtipy.region()
# Physical parameters: power law of slope -2 from gamma=2 to 1e4, B=0.1G, ne=0.002 cm^-3.
# FIXME: Play with these parameter, especially p.
gmin = 2.
gmax = 1e4
p = 2.
B = 0.1
ne = 2e-3
# Create a tuple that contains all of the parameters we call the functions with.
func_params = (B, p, gmin, gmax, ne,)

# Add the synchrotron emission function to medium_a.
medium_a.add_emission_func(j_synchrotron, func_params)
# FIXME: Add the synchrotron absorption function to medium_a.
# medium_a....

# Create an empty initial logarithmically-spaced spectrum going from 200 MHz to 1e25 Hz.
freqaxis = sirtipy.frequency_axis('log', frange=[200*u.MHz, 1e25*u.Hz], numpts=200)
inspec_empty = sirtipy.spectrum('empty', frequencies=freqaxis)

# Transfer through a medium of depth 1e19 cm in steps of 1e16 cm.
s_a = [0., 1e19]
ds = 1e16

# 7. Calculate the radiative transfer to the final distance.
spec_a, tau_a, locations_a = sirtipy.radiative_transfer(medium_a, inspec_empty,
        s_a, ds, printevery=50)

# Now plot the results. First, what does the intensity spectrum look like at the
# different distances?
# 1. List the distances we care about, and figure out which elements they
# correspond to in the output arrays. For example, if distances[1]=1e9, and
# los_to_plot[1]=10, then spec_a[10] is the spectrum after 1e9 cm and tau_a[10,:] is the
# optical depth spectrum after 1e9 cm.
distances = [1e17, 3e17, 1e18, 3e18, 1e19]
los_to_plot = np.searchsorted(locations_a, distances)
# 2. Create a figure.
plt.figure()
# 3. Loop through the different distances we care about so we can plot the
# spectrum at each point.
for l in los_to_plot:
    # This actually plots the spectrum. If a spectrum is already plotted, this
    # overplots.
    spec_a[l].plot(xunit=u.Hz, label=('$s=%s$ cm' % sirtipy.latex_float(locations_a[l])))
# 4. Label the plot and make the axes logarithmic and run over a useful range.
plt.legend(loc='upper right')
# FIXME: Change title if appropriate.
plt.title('Synchrotron (Emission Only)')
plt.xscale('log')
plt.yscale('log')
# 5. Save it to a file.
# FIXME: Change filename.
plt.savefig('3a-sync_emission.png')



## PART B: COMPTON SCATTERING
#
## Compton scattering loss -- sigma x ne, but use Klein-Nishina cross section RL 7.5
## for the cross-section instead of the Thompson value.
## Use the approximations in RL 7.6, and assume a sharp break at x=1 between the two
## limits.
#def alpha_Compton(frequency, location, Inu, B, p, gmin, gmax, ne):
#    x = const.h.cgs.value * frequency / (const.m_e.cgs.value * const.c.cgs.value**2)
#    # The following lines define arrays that distinguish the non-relativistic regime from the
#    # relativistic regime. These will be arrays that are True for every frequency where the
#    # expression is true, and False for every frequency where the expression is false.
#    nonrel = (x <= 1)
#    ultrarel = (x > 1)
## FIXME: Write an expression to compute the Klein-Nishina cross section here.
## NOTE: When you multiply by boolean (True/False) values, True=1 and False=0.
## Therefore, you can split an array into two approximations by multiplying the value
## of one approximation by an array that tells you when that approximation is valid, plus
## the value of the other approximation by an array that tells you when the other approximation
## is valid.
## For example: Consider the effective number of scatterings through a medium, which is:
##   N_scatter = {  tau_s    if tau_s <= 1
##               {  tau_s^2  if tau_s > 1
## We could calculate this as follows:
##   thin_regime = (tau_s <= 1)    # True when tau_s <= 1,   False otherwise
##   thick_regime = (tau_s > 1)    # True when tau_s > 1,    False otherwise
##   N_scatter = thin_regime * tau_s    +     thick_regime * tau_s**2
## When tau_s <=1 this evaluates to:
##                         1 * tau_s    +                0 * tau_s**2
## and when tau_s > 1 this evaluates to:
##                         0 * tau_s    +                1 * tau_s**2
## Use this same logic, and the arrays "nonrel" and "ultrarel" defined above, to stitch together
## the two Klein-Nishina approximations in RL 7.6.
## FIXME: Fill in this line:
##    sigma_KN = ...
#
#    return ne * sigma_KN
#
#
#


# Inverse Compton emitted spectrum.
# From Blumenthal & Gould (1970) equations 2.76-77 
# (like RL 7.29 but including first-order corrections at high photon energy)
def j_Compton(frequency, location, Inu, B, p, gmin, gmax, ne):
    # Figure out normalization C
    C = powerlaw_normalize(p, gmin, gmax, ne)
    # Ap from 7.29b
    Ap = 2.**(p+3) * (p**2 + 4.*p + 11) / ( (p+3.)**2 * (p+5.) * (p+1.) )

    # Do the integral in 7.29a based on Inu
    # The integral is supposed to be over all photon energies, but we only know
    # the radiation field over the specified frequency axis, so we will do that.
    # Convert frequency to energy, and convert intensity to photon number density
    # per unit energy.
    # Intensity: erg / s / Hz / cm2
    #   n: # / erg / cm^3
    # -> n = Inu / c h epsilon
    epsilon = const.h.cgs.value * frequency
    neps = Inu / (const.c.cgs.value * const.h.cgs.value * epsilon)
    # Force non-negativity. Due to finite difference scheme, Inu can end up
    # negative if the absorption coefficient is large relative to the step
    # size, which plays hell with the integral!!
    neps[neps < 0] = 0.
    # Blumenthal & Gould eq 2.77:
    Gp = (p*p * 6*p + 16) * (p+1) * (p+3.)**2 * (p+5.) / ( (p**2 + 4.*p + 11.) * \
        (p+4.)**2 * (p+6.) )
    Gpfac = np.array([Gp * np.sqrt(epsilon*epsilon1) / (const.m_e.cgs.value * const.c.cgs.value**2) \
        for epsilon1 in epsilon])
    Gpfac[Gpfac > 1.] = 1.
    integral = np.array([trapz(epsilon**(0.5*(p-1.)) * neps * (1. - Gpfi), epsilon)  for Gpfi in Gpfac])
    #integral = trapz(epsilon**(0.5*(p-1)) * neps, epsilon)

    # classical electron radius
    r0 = const.e.esu.value**2 / (const.m_e.cgs.value * const.c.cgs.value**2)

    # RL 7.29a. Note that it is expressed per unit energy and we want per unit
    # frequency, so we multiply by an extra factor of h.
    jnu = np.pi * const.c.cgs.value * r0**2 * C * Ap * epsilon**(-0.5*(p-1.)) * integral * const.h.cgs.value

    #pdb.set_trace()

    return jnu

#medium_b = sirtipy.region()
#gmin = 1.
#gmax = 1e6
#p = 2.
#B = 0.1
#ne = 1e-3
#func_params = (B, p, gmin, gmax, ne,)
#s_b = [0., 1e26]
#ds = 1e23
## CMB input
## FIXME: Create a 2.73K blackbody spectrum to use as background. See Experiment 2 as a guide.
##inspec_CMB = ...
#
## Add Compton scattering processes
#medium_b.add_absorption_func(alpha_Compton, func_params)
#medium_b.add_emission_func(j_Compton, func_params)
#
## Perform the radiative transfer.
#spec_b, tau_b, locations_b = sirtipy.radiative_transfer(medium_b, inspec_CMB, s_b, ds, printevery=50)
#
#distances = [0, 1e23, 1e24, 1e25, 1e26]
#los_to_plot = np.searchsorted(locations_b, distances)
#toobig = (los_to_plot > len(locations_b)-1)
#if np.sum(toobig)>0:
#    los_to_plot[toobig] = len(locations_b)-1
## plot
#plt.figure()
## spectrum at each point.
#for l in los_to_plot:
#    # This actually plots the spectrum. If a spectrum is already plotted, this
#    # overplots.
#    spec_b[l].plot(xunit=u.Hz, label=('$s=%s$ cm' % sirtipy.latex_float(locations_b[l])))
## 4. Label the plot and make the axes logarithmic and run over a useful range.
#plt.legend(loc='upper right')
#plt.title('CMB + Inverse Compton')
#plt.xscale('log')
#plt.yscale('log')
## FIXME: Adjust this line, and add a similar plt.xlim(...) line to focus in on particular
## regions of the spectrum.
#plt.ylim(1e-55, 1e-13)
## FIXME: Change filename if necessary.
#plt.savefig('3b.png')

