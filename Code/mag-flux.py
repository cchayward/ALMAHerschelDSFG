"""

Shane Bussmann

2014 August 22

Plot magnification factor vs. S500.  Compare with Schechter function and broken
power-law predictions for the same.

"""

from astropy.table import Table
import matplotlib.pyplot as plt
from pylab import savefig
import matplotlib

# set font properties
font = {'family' : 'Arial',
        'weight' : 'normal',
        'size'   : 12}
matplotlib.rc('font', **font)
matplotlib.rcParams['axes.linewidth'] = 1.5

fig = plt.figure(figsize=(5.0, 4.5))

# ALMA plotting parameters
acolor = 'green'
ams = 4
afmt = 's'

# SMA plotting parameters
bcolor = 'red'
bms = 3
bfmt = 'o'

mudat = Table.read('source_observed.dat', format='ascii')
smadat = Table.read('bussmann2013_muflux.dat', format='ascii')

yyy = mudat['mu870']
xxx = mudat['f870']
yyyerr = mudat['e_mu870']
xxxerr = mudat['e_f870']

plt.errorbar(xxx, yyy, yerr=yyyerr, xerr=xxxerr, fmt=',', ecolor='grey',
        capsize=0)
plt.plot(xxx, yyy, afmt, color=acolor, ms=ams, label='ALMA sample')

yyy = smadat['mu']
xxx = smadat['fnu']
yyyerr = smadat['e_mu']
xxxerr = smadat['e_fnu']

plt.errorbar(xxx, yyy, fmt=',', yerr=yyyerr, xerr=xxxerr, ecolor='grey',
        capsize=0)
plt.plot(xxx, yyy, bfmt, color=bcolor, ms=bms, label='SMA sample')

plt.loglog()

xmin = 0.8
xmax = 2e2
ymin = 0.9
ymax = 35
plt.axis([xmin, xmax, ymin, ymax])

plt.xlabel(r'$S_{\rm 870-observed}\,({\rm mJy})$', fontsize='x-large')
plt.ylabel(r'$\mu_{870}$', fontsize='x-large')
plt.minorticks_on()
plt.tick_params(width=1.2, which='both')
plt.tick_params(length=2, which='minor')
plt.tick_params(length=4, which='major')

plt.legend(loc='upper left', numpoints=1, handletextpad=0.35, borderpad=0.4,
        labelspacing=0.18, handlelength=1.0)
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize='medium')
plt.subplots_adjust(left=0.13, right=0.97, top=0.97, bottom=0.14, wspace=0.39)

savefig('s870_mu.pdf')
