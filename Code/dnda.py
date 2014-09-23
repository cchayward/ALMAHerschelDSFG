"""

Shane Bussmann

2014 August 20

Plot dN/dA as a function of angular separation from the center of light.  dN =
number of objects between radius 1 and radius 2.  dA = area between radius 1
and radius 2.

"""

from astropy.table import Table
import matplotlib
import matplotlib.pyplot as plt
from pylab import savefig
import numpy
import sep_util


# set font properties
font = {'family' : 'Arial',
        'weight' : 'normal',
        'size'   : 12}
matplotlib.rc('font', **font)
matplotlib.rcParams['axes.linewidth'] = 1.2

fig = plt.figure(figsize=(5.0, 4.5))

# Plotting parameters
hodgecolor = 'LightPink'
hodgesimcolor = 'LightPink'
hodgems = 4
hodgefmt = 'D'

nbins = 15
binwidth = 1.0
bin_edges = numpy.arange(0, nbins + 1, binwidth)

# Load the data
fluxcomponent_file = 'hodge2013.dat'
fluxcomponent = Table.read(fluxcomponent_file, format='ascii')

# filter out single source systems
fluxcomponent = sep_util.rmSingles(fluxcomponent, targetstring='lessid')
nmultiples = len(fluxcomponent)

avgsep_hodge, wmeansep_hodge = sep_util.getSeparation(fluxcomponent, 
        targetstring='lessid')

deltasep = avgsep_hodge.max() - avgsep_hodge.min()
#nbins = deltasep / binwidth
sep_util.histArea(avgsep_hodge, nbins, color=hodgecolor, fmt=hodgefmt,
        ms=hodgems, norm=nmultiples)

# plot simulated positions
nsim = 100
sep_util.simArea(fluxcomponent, nsim, bin_edges, targetstring='lessid',
        edgecolor=hodgesimcolor, facecolor='none', hatch='//', norm=nmultiples)

# ***********
# ALMA sample
# ***********

# plotting parameters
acolor = 'green'
asimcolor = 'green'
ams = 5
afmt = 's'

fluxcomponent_file = 'table_positions.dat'
fluxcomponent = Table.read(fluxcomponent_file, format='ascii')

# filter out single source systems
fluxcomponent = sep_util.rmSingles(fluxcomponent, targetstring='target')
nmultiples = len(fluxcomponent)

avgsep_alma, wmeansep_alma = sep_util.getSeparation(fluxcomponent,
        fluxstring='S_870_observed')

sep_util.histArea(avgsep_alma, nbins, color=acolor, fmt=afmt, ms=ams,
        norm=nmultiples)

sep_util.simArea(fluxcomponent, nsim, bin_edges, fluxstring='S_870_observed',
        edgecolor=asimcolor, facecolor='none', hatch='\\', norm=nmultiples)

xmin = 0
ymin = 0
xmax = 6
ymax = 0.25
plt.axis([xmin, xmax, ymin, ymax])

plt.xlabel(r'${\rm Angular\,Separation\,from\,Centroid\,(arcsec)}$', fontsize='large')
plt.ylabel(r'$dN/dA \, ({\rm arcsec}^{-2}$)', fontsize='large')
plt.minorticks_on()
plt.tick_params(width=1.2, which='both')
plt.tick_params(length=2, which='minor')
plt.tick_params(length=4, which='major')

fake = numpy.arange(2) + 1e5
plt.plot(fake, color=hodgecolor, label='Hodge+13 Observed')
plt.plot(fake, color=hodgesimcolor, linestyle='--', 
        label='Hodge+13 Uniform Distribution')
plt.plot(fake, color=acolor, label='ALMA Sample Observed')
plt.plot(fake, color=asimcolor, linestyle='--', 
        label='ALMA Sample Uniform Distribution')
plt.legend(loc='upper right', numpoints=1, handletextpad=0.35, borderpad=0.4,
        labelspacing=0.18, handlelength=1.0)
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
#plt.setp(ltext, fontsize='medium')
plt.subplots_adjust(left=0.14, right=0.95, top=0.97, bottom=0.13, wspace=0.39)

savefig('dNdA.pdf')
import pdb; pdb.set_trace()
