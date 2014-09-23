"""

2014 August 21
Shane Bussmann

Plot the distribution of observed and intrinsic properties for the ALMA sample;
Flux, size, axial ratio, and angular separation.

"""

import matplotlib.pyplot as plt
import numpy
from pylab import savefig
from astropy.table import Table
import matplotlib
import sep_util
from scipy.stats import ks_2samp


# set font properties
font = {'family' : 'Arial',
        'weight' : 'normal',
        'size'   : 12}
matplotlib.rc('font', **font)
matplotlib.rcParams['axes.linewidth'] = 1.5

def makehist(variable, binwidth, xlabel, ylabel, savefile):

    # clear the plotting window
    plt.clf()

    #changed = observed[variable] / intrinsic[variable] != 1
    goodobserved = observed#[changed]
    goodintrinsic = intrinsic#[changed]

    # distribution for observed values
    if variable == 'logf870':
        values = numpy.log10(goodobserved[variable])
        tmpvalues = numpy.log10(goodintrinsic[variable])
    elif variable == 'separation':
        #goodobserved = sep_util.rmSingles(goodobserved)
        #goodintrinsic = sep_util.rmSingles(goodintrinsic)
        values, wmean = sep_util.getSeparation(goodobserved, fluxstring='f870')
        tmpvalues, wm = sep_util.getSeparation(goodintrinsic, fluxstring='f870')
    else:
        values = goodobserved[variable]
        tmpvalues = goodintrinsic[variable]
    m1 = 0
    m2 = values.max()
    tmpm2 = tmpvalues.max()
    if tmpm2 > m2:
        m2 = tmpvalues.max()
    bins = numpy.arange(m1, m2*5, binwidth)
    plt.hist(values, bins = bins, histtype='stepfilled', edgecolor='black',
            facecolor='grey', label='Observed', cumulative=True)
    
    # distribution for intrinsic values
    #m1 = 0
    #m2 = math.ceil(values.max())
    bins = numpy.arange(m1, m2*5, binwidth)
    plt.hist(tmpvalues, bins = bins, histtype='stepfilled', edgecolor='green',
            facecolor='none', label='Intrinsic', cumulative=True, hatch='//',
            linewidth=2)

    xmin = 0
    xmax = m2 * 1.1
    ymin = 0
    ymax = values.size * 1.1
    print(values.size, tmpvalues.size)
    plt.axis([xmin, xmax, ymin, ymax])

    plt.ylabel(ylabel, fontsize='x-large')
    plt.xlabel(xlabel, fontsize='x-large')

    plt.minorticks_on()
    plt.tick_params(width=1.5, which='both')
    plt.tick_params(length=2, which='minor')
    plt.tick_params(length=4, which='major')
    plt.subplots_adjust(left=0.13, right=0.97, top=0.97, bottom=0.13)

    plt.legend(loc='lower right', handletextpad=0.35, borderpad=0.4,
            labelspacing=0.18, handlelength=1.0)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize='large')
    weakly = (values / tmpvalues < 3) & (values/tmpvalues > 1)
    #print(numpy.median(values), numpy.median(tmpvalues))
    print(numpy.median(values[weakly]) / numpy.median(tmpvalues[weakly]))

    # 2-sided KS test
    kstest = ks_2samp(values[weakly], tmpvalues[weakly])
    print(kstest)

    savefig(savefile)

intrinsic = Table.read('source_intrinsic.dat', format='ascii')
observed = Table.read('source_observed.dat', format='ascii')

fig = plt.figure(figsize=(5.0, 4.5))

# 870-micron flux density
#makehist('logf870', 0.01, r'${\rm log}_{10}(S_{870}/{\rm mJy)}$', 'f870_delens.pdf')
makehist('f870', 0.01, r'$S_{870}\,{\rm (mJy)}$', r'$N (S < S_{870})$', 
        'f870_delens.pdf')

# effective radius
makehist('reff', 0.002, r'$r_s \, {\rm (arcsec)}$', r'$N (r < r_s)$', 
        'reff_delens.pdf')

# axial ratio
makehist('q', 0.005, r'$q_s$', r'$N (q < q_s)$', 'q_delens.pdf')

# angular separation from image centroid
makehist('separation', 0.005, 
        r'${\rm Angular\,Offset\,from\,Centroid\,(arcsec)}$', 
        r'$N ({\rm Offset < Offset}_s)$', 'offset_delens.pdf')
