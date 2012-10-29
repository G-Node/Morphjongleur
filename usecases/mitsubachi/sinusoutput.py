# -*- coding: utf-8 -*-
'''
@author: stransky
'''
import morphjongleur.model.morphology
import morphjongleur.model.neuron_passive
import morphjongleur.model.clamp
import morphjongleur.model.experiment
import morphjongleur.util.pattern_generator
import numpy
import matplotlib.pyplot
import scipy.optimize
import numpy

def experiment(morphology, compartment, frequency, amplitude=-1e-9):#, amplitude
    recordingpoints  = [morphjongleur.model.experiment.RecordingPoint(compartment=recordingpoint) for recordingpoint in morphology.terminal_tips]
   
    iclamp  = morphjongleur.util.pattern_generator.SinusClamp(compartment=compartment,
                amplitude=amplitude,
                frequency=frequency, duration=3./frequency
            )
    neuron_passive_parameter    = morphjongleur.model.neuron_passive.Neuron_passive_parameter(Ra=35.4,g=0.001)
    experiment  = morphjongleur.model.experiment.Experiment(
                    morphology=morphology, 
                    recordingpoints=recordingpoints, clamps=[iclamp], 
                    neuron_passive_parameter=neuron_passive_parameter, 
                    duration=5e-3, dt=1e-4,
                    description = "TauExperiment @ %i\t%s " % (morphology.root.compartment_id, morphology.name),
                )
    
    experiment.run_simulation()

    return recordingpoints


exponential_decay   = lambda x,p: p[0] * numpy.exp(-p[1]*x)
exponential_inverse = lambda x,p: p[0] * numpy.exp(-p[1]/x)
err = lambda p, f, x, y: ( f(x,p) - y )**2

def plot_distance_amplitude_scatter(distances, amplitudes, morphology_name, frequency, color='black', ratio=None, picture_file=None, picture_formats=['png', 'pdf', 'svg']):  
    #matplotlib.pyplot.title('amplitude in %s @ %i Hz' % (morphology_name, frequency))
    matplotlib.pyplot.grid(True, color='lightgrey')

    distances   = numpy.array(distances)
    amplitudes  = numpy.array(amplitudes)
    p   = numpy.array([8.,0.0001])
    fp  = exponential_decay
    p, success   = scipy.optimize.leastsq(err, p, args=(fp, distances, amplitudes))
    #errs    = e(p, y, x)
    #avg = numpy.average( errs )
    #mxm = numpy.max( errs )
    #mnm = numpy.min( errs )
    sorted_distances    = numpy.array(sorted(distances))
    matplotlib.pyplot.plot(sorted_distances, fp(sorted_distances,p), color='black')
    matplotlib.pyplot.legend(["fit: $%f \cdot e^{-%f \cdot x}$" % (p[0], p[1])], loc='best')

    matplotlib.pyplot.scatter(distances, amplitudes, c=color, marker='.', edgecolors=color)#'. o
    #matplotlib.pyplot.axes().set_aspect('equal', 'datalim')
    matplotlib.pyplot.ylabel(u'ΔU [mV]')
    matplotlib.pyplot.ylim(ymin=0)
    matplotlib.pyplot.xlabel(u'distance [µm]')
    matplotlib.pyplot.xlim(xmin=0)
 
    if ratio != None:#matplotlib.figure.figaspect(arg)
        fig = matplotlib.pyplot.gcf()
        fig.set_size_inches(ratio[0],ratio[1])
    if(picture_file != None):
        for picture_format in picture_formats:
            matplotlib.pyplot.savefig(picture_file+'.'+picture_format,format=picture_format, transparent=False)
    else:
        matplotlib.pyplot.show()
    matplotlib.pyplot.close()

def plot_frequency_amplitude_scatter(frequencies, amplitudes, compartment_name, distance, morphology_name, color='black', ratio=None, picture_file=None, picture_formats=['png', 'pdf', 'svg']):  
    #matplotlib.pyplot.title(u'amplitude in %s @ %s distance %f µm' % (morphology_name, str(compartment_name), distance))
    matplotlib.pyplot.grid(True, color='lightgrey')

    frequencies = numpy.array(frequencies)
    amplitudes  = numpy.array(amplitudes)
    p   = numpy.array([8.,0.0001])
    fp  = exponential_decay
    p, success   = scipy.optimize.leastsq(err, p, args=(fp, frequencies, amplitudes))
    #errs    = e(p, y, x)
    #avg = numpy.average( errs )
    #mxm = numpy.max( errs )
    #mnm = numpy.min( errs )
    sorted_frequencies  = numpy.array(sorted(frequencies))
    matplotlib.pyplot.plot(sorted_frequencies,fp(sorted_frequencies,p), color='black')
    matplotlib.pyplot.legend(["fit: $%f \cdot e^{-%f \cdot x}$" % (p[0], p[1])], loc='best')

    matplotlib.pyplot.scatter(frequencies, amplitudes, c=color, marker='.', edgecolors=color)#'. o
    #matplotlib.pyplot.axes().set_aspect('equal', 'datalim')
    matplotlib.pyplot.ylabel(u'ΔU [mV]')
    matplotlib.pyplot.ylim(ymin=0)
    matplotlib.pyplot.xlabel('frequency [Hz]')
    #matplotlib.pyplot.xlim(xmin=0)

    if ratio != None:#matplotlib.figure.figaspect(arg)
        fig = matplotlib.pyplot.gcf()
        fig.set_size_inches(ratio[0],ratio[1])
    if(picture_file != None):
        for picture_format in picture_formats:
            matplotlib.pyplot.savefig(picture_file+'.'+picture_format,format=picture_format, transparent=False)
    else:
        matplotlib.pyplot.show()
    matplotlib.pyplot.close()



def plot_amplitude_histogramm(amplitudes, morphology_name, frequency, color='black', ratio=None, picture_file=None, picture_formats=['png', 'pdf', 'svg']):  
    #matplotlib.pyplot.title('amplitude in %s @ %i Hz' % (morphology_name, frequency))
    matplotlib.pyplot.grid(True, color='lightgrey')

    matplotlib.pyplot.hist(amplitudes, bins=20, normed=0, color=color)
    matplotlib.pyplot.ylabel('\#')
    #matplotlib.pyplot.ylim(0, 130)#0.3
    matplotlib.pyplot.xlabel(u'ΔU [mV]')
    #matplotlib.pyplot.xlim(0, 45)

    if ratio != None:#matplotlib.figure.figaspect(arg)
        fig = matplotlib.pyplot.gcf()
        fig.set_size_inches(ratio[0],ratio[1])
    if(picture_file != None):
        for picture_format in picture_formats:
            matplotlib.pyplot.savefig(picture_file+'.'+picture_format,format=picture_format, transparent=False)
    else:
        matplotlib.pyplot.show()
    matplotlib.pyplot.close()

def plot_duration_histogramm(durations, morphology_name, frequency, color='black', ratio=None, picture_file=None, picture_formats=['png', 'pdf', 'svg']):  
    #matplotlib.pyplot.title('duration in %s @ %i Hz ' % (morphology_name, frequency))
    matplotlib.pyplot.ylabel('\#')
    #matplotlib.pyplot.ylim(0, 70)
    matplotlib.pyplot.xlabel('duration [ms]')
    #matplotlib.pyplot.xlim(0, 2.5)
    matplotlib.pyplot.xticks()
    matplotlib.pyplot.grid(True, color='lightgrey')

    matplotlib.pyplot.hist(durations, bins=20, normed=0, color=color)

    if ratio != None:#matplotlib.figure.figaspect(arg)
        fig = matplotlib.pyplot.gcf()
        fig.set_size_inches(ratio[0],ratio[1])
    if(picture_file != None):
        for picture_format in picture_formats:
            matplotlib.pyplot.savefig(picture_file+'.'+picture_format, format=picture_format, transparent=False)
    else:
        matplotlib.pyplot.show()
    matplotlib.pyplot.close()

def plot_delay_histogramm(delays, morphology_name, frequency, color='black', ratio=None, picture_file=None, picture_formats=['png', 'pdf', 'svg']):  
    #matplotlib.pyplot.title('delay in %s @ %i Hz ' % (morphology_name, frequency))
    matplotlib.pyplot.xticks()
    matplotlib.pyplot.grid(True, color='lightgrey')

    matplotlib.pyplot.hist(delays, bins=20, normed=0, color=color)
    matplotlib.pyplot.ylabel('\#')
    matplotlib.pyplot.ylim(0, 300)
    matplotlib.pyplot.xlabel('delay [ms]')
    matplotlib.pyplot.xlim(0, 1000./frequency)

    if ratio != None:#matplotlib.figure.figaspect(arg)
        fig = matplotlib.pyplot.gcf()
        fig.set_size_inches(ratio[0],ratio[1])
    if(picture_file != None):
        for picture_format in picture_formats:
            matplotlib.pyplot.savefig(picture_file+'.'+picture_format, format=picture_format, transparent=False)
    else:
        matplotlib.pyplot.show()
    matplotlib.pyplot.close()

def plot_phaseshift_histogramm(phaseshifts, morphology_name, frequency, color='black', ratio=None, picture_file=None, picture_formats=['png', 'pdf', 'svg']):  
    #matplotlib.pyplot.title('phaseshift in %s @ %i Hz ' % (morphology_name, frequency))
    matplotlib.pyplot.grid(True, color='lightgrey')

    matplotlib.pyplot.hist(phaseshifts, bins=20, normed=0, color=color)
    matplotlib.pyplot.ylabel('\#')
    #matplotlib.pyplot.ylim(0, 300)#
    matplotlib.pyplot.xlabel('phase angle [full circles]')
    matplotlib.pyplot.xlim(0, 1)

    if ratio != None:#matplotlib.figure.figaspect(arg)
        fig = matplotlib.pyplot.gcf()
        fig.set_size_inches(ratio[0],ratio[1])
    if(picture_file != None):
        for picture_format in picture_formats:
            matplotlib.pyplot.savefig(picture_file+'.'+picture_format, format=picture_format, transparent=False)
    else:
        matplotlib.pyplot.show()
    matplotlib.pyplot.close()

#===============================================================================
# def plot_candlestick(results, morphology_name='', color='black', ratio=None, yline=None, fit_function=lambda x,p: p[0] * numpy.exp(-p[1]*x), error_function=lambda p, f, x, y: ( f(x,p) - y )**2, picture_file=None, picture_formats=['png', 'pdf', 'svg']):
#    matplotlib.rc('text', usetex=True)
# 
#    x = []
#    y = []
# 
#    for x_value in sorted(results.keys()):
#        for y_value in durations[x_value]:
#            x.append(x_value)
#            y.append(y_value)
#    x       = numpy.array(x)
#    y       = numpy.array(y)
#    p   = numpy.array([8.,0.0001])
#    fp  = exponential_decay
#    p, success   = scipy.optimize.leastsq(err, p, args=(fp, x, y))#, maxfev=10000*len(x)
#    print "duration %s" % (str(p))
#    #err = e(p, y, x)
#    #avg = numpy.average( errs )
#    #mxm = numpy.max( errs )
#    #mnm = numpy.min( errs )
#    matplotlib.pyplot.plot(x,fp(x,p), color=color)#numpy.arange(85,510,1)
# 
#    titles=[]
#    titles.append( "$%.0fV \cdot e^{-%.8fs \cdot \mathrm{f}}$" % (p[0], p[1]) )
#    matplotlib.pyplot.legend(titles, loc='best')
# 
#    #matplotlib.pyplot.title('Duration in '+str(morphology))
#    matplotlib.pyplot.xticks(rotation=90)#  xrange(1+len(data)),titles
#    matplotlib.pyplot.grid(True, color='lightgrey')
# 
#    bp = matplotlib.pyplot.boxplot(durations.values(), positions=durations.keys(), widths=min_dist*2./3)
#    matplotlib.pyplot.ylabel('duration [ms]')
#    #matplotlib.pyplot.ylim(0, 3)
#    matplotlib.pyplot.xlabel('frequency [Hz]')
#    #matplotlib.pyplot.xlim(85, 510)
#    #matplotlib.pyplot.setp(bp['whiskers'], color='k',  linestyle='-' )
#    #matplotlib.pyplot.setp(bp['fliers'], color='k')
# 
#    if yline != None:#matplotlib.figure.figaspect(arg)
#            matplotlib.pyplot.axvline(x=yline, color='darkgray')
#    if ratio != None:#matplotlib.figure.figaspect(arg)
#        fig = matplotlib.pyplot.gcf()
#        fig.set_size_inches(ratio[0],ratio[1])
#    if(picture_file != None):
#        for picture_format in picture_formats:
#            matplotlib.pyplot.savefig(picture_file+'.'+picture_format, format=picture_format, transparent=False)
#    else:
#        matplotlib.pyplot.show()
#    matplotlib.pyplot.close()
#    return p
#===============================================================================

def plot_duration_candlestick(durations, morphology_name='', color='black', ratio=None, yline=None, picture_file=None, picture_formats=['png', 'pdf', 'svg']):
    matplotlib.rc('text', usetex=True)

    x = []
    y = []
    last_value   = float('+inf')
    min_dist    = float('+inf') 
    for frequency in sorted(durations.keys()):
        if abs(frequency - last_value) < min_dist:
            min_dist    =  abs(frequency - last_value)
        last_value    = frequency
        for delay in durations[frequency]:
            x.append(frequency)
            y.append(delay)
    x       = numpy.array(x)
    y       = numpy.array(y)
    p   = numpy.array([8.,0.0001])
    fp  = exponential_decay
    p, success   = scipy.optimize.leastsq(err, p, args=(fp, x, y))#, maxfev=10000*len(x)
    print "duration %s" % (str(p))
    #err = e(p, y, x)
    #avg = numpy.average( errs )
    #mxm = numpy.max( errs )
    #mnm = numpy.min( errs )
    matplotlib.pyplot.plot(x,fp(x,p), color=color)#numpy.arange(85,510,1)

    titles=[]
    titles.append( "$%.0fV \cdot e^{-%.8fs \cdot \mathrm{f}}$" % (p[0], p[1]) )
    matplotlib.pyplot.legend(titles, loc='best')

    #matplotlib.pyplot.title('Duration in '+str(morphology))
    matplotlib.pyplot.xticks(rotation=90)#  xrange(1+len(data)),titles
    matplotlib.pyplot.grid(True, color='lightgrey')
 
    bp = matplotlib.pyplot.boxplot(durations.values(), positions=durations.keys(), widths=min_dist*2./3)
    matplotlib.pyplot.ylabel('duration [ms]')
    #matplotlib.pyplot.ylim(0, 3)
    matplotlib.pyplot.xlabel('frequency [Hz]')
    #matplotlib.pyplot.xlim(85, 510)
    #matplotlib.pyplot.setp(bp['whiskers'], color='k',  linestyle='-' )
    #matplotlib.pyplot.setp(bp['fliers'], color='k')

    if yline != None:#matplotlib.figure.figaspect(arg)
            matplotlib.pyplot.axvline(x=yline, color='darkgray')
    if ratio != None:#matplotlib.figure.figaspect(arg)
        fig = matplotlib.pyplot.gcf()
        fig.set_size_inches(ratio[0],ratio[1])
    if(picture_file != None):
        for picture_format in picture_formats:
            matplotlib.pyplot.savefig(picture_file+'.'+picture_format, format=picture_format, transparent=False)
    else:
        matplotlib.pyplot.show()
    matplotlib.pyplot.close()
    return p

def plot_delay_candlestick(delays, morphology_name='', color='black', ratio=None, yline=None, picture_file=None, picture_formats=['png', 'pdf', 'svg']):
    matplotlib.rc('text', usetex=True)

    x = []
    y = []
    last_value   = float('+inf')
    min_dist    = float('+inf') 
    for frequency in sorted(delays.keys()):
        if abs(frequency - last_value) < min_dist:
            min_dist    =  abs(frequency - last_value)
        last_value    = frequency
        for delay in delays[frequency]:
            x.append(frequency)
            y.append(delay)
    x       = numpy.array(x)
    y       = numpy.array(y)
    p   = numpy.array([8.,0.0001])
    fp  = exponential_decay
    p, success   = scipy.optimize.leastsq(err, p, args=(fp, x, y))#, maxfev=10000*len(x)
    print "delays %s" % (str(p))
    #err = e(p, y, x)
    #avg = numpy.average( errs )
    #mxm = numpy.max( errs )
    #mnm = numpy.min( errs )
    matplotlib.pyplot.plot(x,fp(x,p), color=color)#numpy.arange(85,510,1)

    titles=[]
    titles.append( "$%.2fV \cdot e^{-%.8fs \cdot \mathrm{f}}$" % (p[0], p[1]) )
    matplotlib.pyplot.legend(titles, loc='best')

    #matplotlib.pyplot.title('Delay in '+str(morphology))
    matplotlib.pyplot.xticks(rotation=90)#  xrange(1+len(data)),titles
    matplotlib.pyplot.grid(True, color='lightgrey')
 
    bp = matplotlib.pyplot.boxplot(delays.values(), positions=delays.keys(), widths=min_dist*2./3)
    matplotlib.pyplot.ylabel('delays [ms]')
    #matplotlib.pyplot.ylim(0, 3)
    matplotlib.pyplot.xlabel('frequency [Hz]')
    #matplotlib.pyplot.xlim(85, 510)
    #matplotlib.pyplot.setp(bp['whiskers'], color='k',  linestyle='-' )
    #matplotlib.pyplot.setp(bp['fliers'], color='k')

    if yline != None:#matplotlib.figure.figaspect(arg)
            matplotlib.pyplot.axvline(x=yline, color='darkgray')
    if ratio != None:#matplotlib.figure.figaspect(arg)
        fig = matplotlib.pyplot.gcf()
        fig.set_size_inches(ratio[0],ratio[1])
    if(picture_file != None):
        for picture_format in picture_formats:
            matplotlib.pyplot.savefig(picture_file+'.'+picture_format, format=picture_format, transparent=False)
    else:
        matplotlib.pyplot.show()
    matplotlib.pyplot.close()
    return p

def plot_phaseshift_candlestick(phaseshifts, morphology_name='', color='black', ratio=None, yline=None, picture_file=None, picture_formats=['png', 'pdf', 'svg']):
    matplotlib.rc('text', usetex=True)

    x = []
    y = []
    last_value   = float('+inf')
    min_dist    = float('+inf') 
    for frequency in sorted(phaseshifts.keys()):
        if abs(frequency - last_value) < min_dist:
            min_dist    =  abs(frequency - last_value)
        last_value    = frequency
        for phase_angle in phaseshifts[frequency]:
            x.append(frequency)
            y.append(phase_angle)
    x       = numpy.array(x)
    y       = numpy.array(y)
    p   = numpy.array([8.,0.0001])
    fp  = exponential_decay
    p, success   = scipy.optimize.leastsq(err, p, args=(fp, x, y))#, maxfev=10000*len(x)
    print "phase shift %s" % (str(p))
    #err = e(p, y, x)
    #avg = numpy.average( errs )
    #mxm = numpy.max( errs )
    #mnm = numpy.min( errs )
    #matplotlib.pyplot.plot(x,fp(x,p), color=color)#numpy.arange(85,510,1)

    #titles=[]
    #titles.append( "$%.2fV \cdot e^{-%.8fs \cdot \mathrm{f}}$" % (p[0], p[1]) )
    #matplotlib.pyplot.legend(titles, loc='best')


    #matplotlib.pyplot.title('Phaseshift in '+str(morphology))
    matplotlib.pyplot.xticks(rotation=90)#  xrange(1+len(data)),titles
    matplotlib.pyplot.grid(True, color='lightgrey')
 
    bp = matplotlib.pyplot.boxplot(phaseshifts.values(), positions=phaseshifts.keys(), widths=min_dist*2./3)
    matplotlib.pyplot.ylabel('phase angle [full circles]')
    matplotlib.pyplot.ylim(0, 1)
    matplotlib.pyplot.xlabel('frequency [Hz]')
    #matplotlib.pyplot.xlim(85, 510)
    #matplotlib.pyplot.setp(bp['whiskers'], color='k',  linestyle='-' )
    #matplotlib.pyplot.setp(bp['fliers'], color='k')

    if yline != None:#matplotlib.figure.figaspect(arg)
            matplotlib.pyplot.axvline(x=yline, color='darkgray')
    if ratio != None:#matplotlib.figure.figaspect(arg)
        fig = matplotlib.pyplot.gcf()
        fig.set_size_inches(ratio[0],ratio[1])
    if(picture_file != None):
        for picture_format in picture_formats:
            matplotlib.pyplot.savefig(picture_file+'.'+picture_format, format=picture_format, transparent=False)
    else:
        matplotlib.pyplot.show()
    matplotlib.pyplot.close()
    return p
 
def plot_amplitude_s_candlestick(amplitudes, morphology_name='', color='black', ratio=None, yline=None, picture_file=None, picture_formats=['png', 'pdf', 'svg']):
    matplotlib.rc('text', usetex=True)

    x = []
    y = []
    last_value   = float('+inf')
    min_dist    = float('+inf') 
    for frequency in sorted(amplitudes.keys()):
        if abs(frequency - last_value) < min_dist:
            min_dist    =  abs(frequency - last_value)
        last_value    = frequency
        for voltage in amplitudes[frequency]:
            x.append(1./frequency)
            y.append(voltage)
    x       = numpy.array(x)
    y       = numpy.array(y)
    p   = numpy.array([1.,0.5])
    fp  = exponential_inverse
    p, success   = scipy.optimize.leastsq(err, p, args=(fp, x, y))
    print "amplitude_s %s" % (str(p))
    #err = e(p, y, x)
    #avg = numpy.average( errs )
    #mxm = numpy.max( errs )
    #mnm = numpy.min( errs )
    matplotlib.pyplot.plot(x,fp(x,p), color=color)

    titles=[]
    titles.append( "fit: $%.2f \cdot e^{-%.2f / \mathrm{t}}$" % (p[0], p[1]) )
    matplotlib.pyplot.legend(titles, loc='best')

    matplotlib.pyplot.title('amplitude in '+str(morphology_name.replace('_','\_')))
    matplotlib.pyplot.xticks(rotation=90)#  xrange(1+len(data)),titles
    matplotlib.pyplot.grid(True, color='lightgrey')
    
    bp = matplotlib.pyplot.boxplot(amplitudes.values(), positions=[1./frequency for frequency in amplitudes.keys()], widths=0.0001)
    matplotlib.pyplot.ylabel(u'ΔU [mV]')
    #matplotlib.pyplot.ylim(0, 40)
    matplotlib.pyplot.xlabel('time [s]')
    matplotlib.pyplot.xlim(0, 1./100)
    #matplotlib.pyplot.setp(bp['whiskers'], color='k',  linestyle='-' )
    #matplotlib.pyplot.setp(bp['fliers'], color='k')

    if yline != None:#matplotlib.figure.figaspect(arg)
            matplotlib.pyplot.axvline(x=yline, color='darkgray')
    if ratio != None:#matplotlib.figure.figaspect(arg)
        fig = matplotlib.pyplot.gcf()
        fig.set_size_inches(ratio[0],ratio[1])
    if(picture_file != None):
        for picture_format in picture_formats:
            matplotlib.pyplot.savefig(picture_file+'.'+picture_format, format=picture_format, transparent=False)
    else:
        matplotlib.pyplot.show()
    matplotlib.pyplot.close('all')
    return p

def plot_amplitude_candlestick(amplitudes, morphology_name='', color='black', ratio=None, yline=None, picture_file=None, picture_formats=['png', 'pdf', 'svg']):
    matplotlib.rc('text', usetex=True)

    x = []
    y = []
    last_value   = float('+inf')
    min_dist    = float('+inf') 
    for frequency in sorted(amplitudes.keys()):
        if abs(frequency - last_value) < min_dist:
            min_dist    =  abs(frequency - last_value)
        last_value    = frequency
        for voltage in amplitudes[frequency]:
            x.append(frequency)
            y.append(voltage)
        print "%f %s" % (frequency,amplitudes[frequency])
    x       = numpy.array(x)
    y       = numpy.array(y)
    p   = numpy.array([8.,0.0001])
    fp  = exponential_decay
    p, success   = scipy.optimize.leastsq(err, p, args=(fp, x, y))#, maxfev=10000*len(x)
    print "amplitude %s" % (str(p))
    #err = e(p, y, x)
    #avg = numpy.average( errs )
    #mxm = numpy.max( errs )
    #mnm = numpy.min( errs )
    matplotlib.pyplot.plot(x,fp(x,p), color=color)#numpy.arange(85,510,1)

    titles=[]
    titles.append( "$%.0fV \cdot e^{-%.8fs \cdot \mathrm{f}}$" % (p[0], p[1]) )
    matplotlib.pyplot.legend(titles, loc='best')

    #matplotlib.pyplot.title('amplitude in '+str(morphology_name.replace('_','\_')))
    matplotlib.pyplot.xticks(rotation=90)#  xrange(1+len(data)),titles
    matplotlib.pyplot.grid(True, color='lightgrey')

    bp = matplotlib.pyplot.boxplot(amplitudes.values(), positions=amplitudes.keys(), widths=min_dist*2./3)
    matplotlib.pyplot.ylabel(u'voltage amplitude [mV]')
    #matplotlib.pyplot.ylim(0, 40)
    matplotlib.pyplot.xlabel('frequency [Hz]')
    #matplotlib.pyplot.xlim(x[0]-10, x[-1]+10)
    #matplotlib.pyplot.setp(bp['whiskers'], color='k',  linestyle='-' )
    #matplotlib.pyplot.setp(bp['fliers'], color='k')

    if yline != None:#matplotlib.figure.figaspect(arg)
            matplotlib.pyplot.axvline(x=yline, color='darkgray')
    if ratio != None:#matplotlib.figure.figaspect(arg)
        fig = matplotlib.pyplot.gcf()
        fig.set_size_inches(ratio[0],ratio[1])
    if(picture_file != None):
        for picture_format in picture_formats:
            matplotlib.pyplot.savefig(picture_file+'.'+picture_format, format=picture_format, transparent=False)
    else:
        matplotlib.pyplot.show()
    matplotlib.pyplot.close('all')
    return

def plot_frequency_distance_amplitude_scatter3d(frequencies, distances, amplitudes, morphology_name, colors=['black'], picture_file=None, picture_formats=['png', 'pdf', 'svg']):  
    import matplotlib.pyplot
    from mpl_toolkits.mplot3d import Axes3D
    
    fig = matplotlib.pyplot.figure()
    ax = fig.add_subplot(111, projection='3d')
    i = 0
    #TODO: assert equal dimension
    for key in sorted(distances, key=distances.get):
        ax.scatter(frequencies[key], len(frequencies[key]) * [distances[key]], amplitudes[key], c=range(len(frequencies[key])))#, marker="$%i$" % (i)
        i = i+1

    ax.set_zlabel('voltage amplitude [mV]')
    #ax.set_zlim(0, 45)
    ax.set_ylabel(u'distance [µm]')
    #ax.set_ylim(0, 900)
    ax.set_xlabel('frequency [Hz]')
    #ax.set_xlim(100, 500)

    
    ratio=(5,6)
    if ratio != None:#matplotlib.figure.figaspect(arg)
        fig = matplotlib.pyplot.gcf()
        fig.set_size_inches(ratio[0],ratio[1])
    if(picture_file != None):
        for picture_format in picture_formats:
            matplotlib.pyplot.savefig(picture_file+'.'+picture_format, format=picture_format, transparent=False)
    else:
        matplotlib.pyplot.show()
    matplotlib.pyplot.close('all')
    #return p

if __name__ == '__main__':
    '''
    Parameter: files, not directories

    cd data/amira/
    pathto/metric_analysis.py *.swc
    '''
    import sys
    import numpy
    import morphjongleur.util.parser.swc
    import morphjongleur.util.transformations
    import morphjongleur.model.morphology
    picture_formats = ['png']#,'pdf','svg'
    c = 'black'
    colors  = ['black','black','#00ff00','#008000', '#0000ff','#000080']
    frequencies = xrange(100,501,15)#xrange(490,9,-15)   xrange(1,10000,1)    [1/t for t in xrange(1e-4,1e-4, 1e0)]
    for swc in sys.argv[1:]:
        print swc, 
        amplitues = {}
        delays    = {}
        durations = {}
        sinus_clamps    = []
        for f in frequencies:
            print f
            m   = morphjongleur.model.morphology.Morphology.swc_parse(swc, verbose=False)
            amplitues[f]    = [ ]
            delays[f]       = [ ]
            durations[f]    = [ ]
            distances       = [ ]#TODO: calculate _once_
            recordingpoints = experiment(morphology=m, compartment=m.root, frequency=f)
            assert len(recordingpoints[0].experiment.iclamps) == 1
            sinus_clamp =   recordingpoints[0].experiment.iclamps[0]
            sinus_clamp.plot(picture_file='/tmp/SinusClamp_%i' % (sinus_clamp.frequency), picture_formats=['png'])
            sinus_clamps.append( sinus_clamp )
            for recordingpoint in recordingpoints:
                voltage_trace   = recordingpoint.get_voltage_trace()
                amplitues[f].append( voltage_trace.amplitude )
                delays[f].append(    voltage_trace.delay )
                durations[f].append( voltage_trace.duration )
                distances.append( m.root.distance_path( recordingpoint.compartment ) )

            plot_distance_amplitude_scatter(distances=distances, amplitudes=amplitues[f], morphology_name=m.name, frequency=f, color=c, ratio=(5,6), picture_file='/tmp/distance_amplitude_scatter_%s%iHz' %(m.name,f), picture_formats=picture_formats )
            plot_amplitude_histogramm(amplitudes=amplitues[f], morphology_name=m.name, frequency=f, color=c, ratio=(5,6), picture_file='/tmp/voltage_histogramm_%s%iHz'  %(m.name,f), picture_formats=picture_formats )
            plot_duration_histogramm(durations=durations[f], morphology_name=m.name, frequency=f, color=c, ratio=(5,6), picture_file='/tmp/duration_histogramm_%s%iHz' %(m.name,f), picture_formats=picture_formats )
            plot_delay_histogramm(delays=delays[f], morphology_name=m.name, frequency=f, color=c, ratio=(5,6), picture_file='/tmp/delay_histogramm_%s%iHz'    %(m.name,f), picture_formats=picture_formats )
            plot_phaseshift_histogramm(phaseshifts=[f * delay for delay in delays[f]], morphology_name=m.name, frequency=f, color=c, ratio=(5,6), picture_file='/tmp/phase_histogramm_%s%iHz'    %(m.name,f), picture_formats=picture_formats )

        print ''
        morphjongleur.model.clamp.PatternClamp.plots(pattern_clamps=sinus_clamps, ratio=(100,10), picture_file='/tmp/SinusClamp', picture_formats=['png'])
            #redundat but show instant process
        plot_amplitude_candlestick( amplitudes=amplitues, morphology_name=m.name, color='black', ratio=(5,6), yline=265, picture_file='/tmp/amplitude_out_candlestick_%s' % (m.name),    picture_formats=picture_formats )
        plot_duration_candlestick(  durations=durations, morphology_name=m.name, color='black', ratio=(5,6), yline=265, picture_file='/tmp/duration_candlestick_%s' %(m.name),      picture_formats=picture_formats )
        plot_delay_candlestick(     delays=delays, morphology_name=m.name, color='black', ratio=(5,6), yline=265, picture_file='/tmp/delay_candlestick_%s' %(m.name),         picture_formats=picture_formats )
        
        ps  = {}
        for f,ds in delays.iteritems():
            ps[f]   = [f * delay for delay in ds]
        plot_phaseshift_candlestick(phaseshifts=ps, morphology_name=m.name, color='black', ratio=(5,6), yline=265, picture_file='/tmp/phase_candlestick_%s' %(m.name),         picture_formats=picture_formats )
        plot_amplitude_s_candlestick( amplitudes=amplitues, morphology_name=m.name, color='black', ratio=(5,6), yline=265, picture_file='/tmp/amplitude_s_candlestick_%s' % (m.name),  picture_formats=picture_formats )
        #plot_frequency_distance_amplitude_scatter3d(frequencies=frequencies, distances=distances, amplitudes=amplitues, morphology_name=m.name, color=c, ratio=(5,6), picture_file='/tmp/frequency_distance_amplitude_scatter3d_%s'%(m.name), picture_formats=picture_formats )
