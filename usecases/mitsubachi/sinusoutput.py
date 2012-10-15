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
   
   
    iclamp  = morphjongleur.model.clamp.IClamp(compartment=morphology.root,
                amplitude=-1e-9, delay=0e-3, duration=3./frequency
            )
   
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
    amplitues = [recordingpoint.get_voltage_trace().amplitude for recordingpoint in recordingpoints]
    return amplitues


exponential_decay   = lambda x,p: p[0] * numpy.exp(-p[1]*x)

err = lambda p, f, x, y: ( f(x,p) - y )**2

def plot_distance_amplitude_scatter(distances, amplitudes, morphology_name, frequency, color='black', picture_file=None, picture_formats=['png', 'pdf', 'svg']):  
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
 
    
    ratio=(5,6)
    if ratio != None:#matplotlib.figure.figaspect(arg)
        fig = matplotlib.pyplot.gcf()
        fig.set_size_inches(ratio[0],ratio[1])
    if(picture_file != None):
        for picture_format in picture_formats:
            matplotlib.pyplot.savefig(picture_file+'.'+picture_format,format=picture_format, transparent=False)
    else:
        matplotlib.pyplot.show()
    matplotlib.pyplot.close()

def plot_frequency_amplitude_scatter(frequencies, amplitudes, compartment_name, distance, morphology_name, color='black', picture_file=None, picture_formats=['png', 'pdf', 'svg']):  
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
 
    
    ratio=(5,6)
    if ratio != None:#matplotlib.figure.figaspect(arg)
        fig = matplotlib.pyplot.gcf()
        fig.set_size_inches(ratio[0],ratio[1])
    if(picture_file != None):
        for picture_format in picture_formats:
            matplotlib.pyplot.savefig(picture_file+'.'+picture_format,format=picture_format, transparent=False)
    else:
        matplotlib.pyplot.show()
    matplotlib.pyplot.close()



def plot_amplitude_histogramm(amplitudes, morphology_name, frequency, color='black', picture_file=None, picture_formats=['png', 'pdf', 'svg']):  
    #matplotlib.pyplot.title('amplitude in %s @ %i Hz' % (morphology_name, frequency))
    matplotlib.pyplot.grid(True, color='lightgrey')

    matplotlib.pyplot.hist(amplitudes, bins=20, normed=0, color=color)
    matplotlib.pyplot.ylabel('\#')
    matplotlib.pyplot.ylim(0, 130)#0.3
    matplotlib.pyplot.xlabel(u'ΔU [mV]')
    matplotlib.pyplot.xlim(0, 45)
 
    
    ratio=(5,6)
    if ratio != None:#matplotlib.figure.figaspect(arg)
        fig = matplotlib.pyplot.gcf()
        fig.set_size_inches(ratio[0],ratio[1])
    if(picture_file != None):
        for picture_format in picture_formats:
            matplotlib.pyplot.savefig(picture_file+'.'+picture_format,format=picture_format, transparent=False)
    else:
        matplotlib.pyplot.show()
    matplotlib.pyplot.close()

def plot_duration_histogramm(durations, morphology_name, frequency, color='black', picture_file=None, picture_formats=['png', 'pdf', 'svg']):  
    #matplotlib.pyplot.title('duration in %s @ %i Hz ' % (morphology_name, frequency))
    matplotlib.pyplot.ylabel('\#')
    #matplotlib.pyplot.ylim(0, 70)
    matplotlib.pyplot.xlabel('duration [ms]')
    #matplotlib.pyplot.xlim(0, 2.5)
    matplotlib.pyplot.xticks()
    matplotlib.pyplot.grid(True, color='lightgrey')

    matplotlib.pyplot.hist(durations, bins=20, normed=0, color=color)

    
    ratio=(5,6)
    if ratio != None:#matplotlib.figure.figaspect(arg)
        fig = matplotlib.pyplot.gcf()
        fig.set_size_inches(ratio[0],ratio[1])
    if(picture_file != None):
        for picture_format in picture_formats:
            matplotlib.pyplot.savefig(picture_file+'.'+picture_format, format=picture_format, transparent=False)
    else:
        matplotlib.pyplot.show()
    matplotlib.pyplot.close()

def plot_delay_histogramm(delays, morphology_name, frequency, color='black', picture_file=None, picture_formats=['png', 'pdf', 'svg']):  
    #matplotlib.pyplot.title('delay in %s @ %i Hz ' % (morphology_name, frequency))
    matplotlib.pyplot.xticks()
    matplotlib.pyplot.grid(True, color='lightgrey')

    matplotlib.pyplot.hist(delays, bins=20, normed=0, color=color)
    matplotlib.pyplot.ylabel('\#')
    matplotlib.pyplot.ylim(0, 300)
    matplotlib.pyplot.xlabel('delay [ms]')
    matplotlib.pyplot.xlim(0, 1000./frequency)

    
    ratio=(5,6)
    if ratio != None:#matplotlib.figure.figaspect(arg)
        fig = matplotlib.pyplot.gcf()
        fig.set_size_inches(ratio[0],ratio[1])
    if(picture_file != None):
        for picture_format in picture_formats:
            matplotlib.pyplot.savefig(picture_file+'.'+picture_format, format=picture_format, transparent=False)
    else:
        matplotlib.pyplot.show()
    matplotlib.pyplot.close()

def plot_phaseshift_histogramm(phaseshifts, morphology_name, frequency, color='black', picture_file=None, picture_formats=['png', 'pdf', 'svg']):  
    #matplotlib.pyplot.title('phaseshift in %s @ %i Hz ' % (morphology_name, frequency))
    matplotlib.pyplot.grid(True, color='lightgrey')

    matplotlib.pyplot.hist(phaseshifts, bins=20, normed=0, color=color)
    matplotlib.pyplot.ylabel('\#')
    matplotlib.pyplot.ylim(0, 300)#
    matplotlib.pyplot.xlabel('phase angle [full circles]')
    matplotlib.pyplot.xlim(0, 1)

    
    ratio=(5,6)
    if ratio != None:#matplotlib.figure.figaspect(arg)
        fig = matplotlib.pyplot.gcf()
        fig.set_size_inches(ratio[0],ratio[1])
    if(picture_file != None):
        for picture_format in picture_formats:
            matplotlib.pyplot.savefig(picture_file+'.'+picture_format, format=picture_format, transparent=False)
    else:
        matplotlib.pyplot.show()
    matplotlib.pyplot.close()


def plot_duration_candlestick(sinus_results, morphology, color='black', picture_file=None, picture_formats=['png', 'pdf', 'svg']):
    f   = {}
    for r in sinus_results:
        if not f.has_key(r.frequency):
            f[r.frequency] = []
        f[r.frequency].append( r.duration )
    x = []
    y = []
    for frequency in sorted(f.keys()):
        for delay in f[frequency]:
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
 
    bp = matplotlib.pyplot.boxplot(f.values(), positions=f.keys(), widths=10)
    matplotlib.pyplot.ylabel('duration [ms]')
    #matplotlib.pyplot.ylim(0, 3)
    matplotlib.pyplot.xlabel('frequency [Hz]')
    #matplotlib.pyplot.xlim(85, 510)
    #matplotlib.pyplot.setp(bp['whiskers'], color='k',  linestyle='-' )
    #matplotlib.pyplot.setp(bp['fliers'], color='k')

    
    ratio=(5,6)
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

def plot_delay_candlestick(sinus_results, morphology, color='black', picture_file=None, picture_formats=['png', 'pdf', 'svg']):
    f   = {}
    for r in sinus_results:
        if not f.has_key(r.frequency):
            f[r.frequency] = []
        f[r.frequency].append( r.delay )
    x = []
    y = []
    for frequency in sorted(f.keys()):
        for delay in f[frequency]:
            x.append(frequency)
            y.append(delay)
    x       = numpy.array(x)
    y       = numpy.array(y)
    p   = numpy.array([8.,0.0001])
    fp  = exponential_decay
    p, success   = scipy.optimize.leastsq(err, p, args=(fp, x, y))#, maxfev=10000*len(x)
    print "delay %s" % (str(p))
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
 
    bp = matplotlib.pyplot.boxplot(f.values(), positions=f.keys(), widths=10)
    matplotlib.pyplot.ylabel('delay [ms]')
    matplotlib.pyplot.ylim(0, 3)
    matplotlib.pyplot.xlabel('frequency [Hz]')
    matplotlib.pyplot.xlim(85, 510)
    #matplotlib.pyplot.setp(bp['whiskers'], color='k',  linestyle='-' )
    #matplotlib.pyplot.setp(bp['fliers'], color='k')

    
    ratio=(5,6)
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

def plot_phaseshift_candlestick(sinus_results, morphology, color='black', picture_file=None, picture_formats=['png', 'pdf', 'svg']):
    f   = {}
    for r in sinus_results:
        if not f.has_key(r.frequency):
            f[r.frequency] = []
        f[r.frequency].append(r.phase_angle)
    x = []
    y = []
    for frequency in sorted(f.keys()):
        for phase_angle in f[frequency]:
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
 
    bp = matplotlib.pyplot.boxplot(f.values(), positions=f.keys(), widths=10)
    matplotlib.pyplot.ylabel('phase angle [full circles]')
    matplotlib.pyplot.ylim(0, 1)
    matplotlib.pyplot.xlabel('frequency [Hz]')
    #matplotlib.pyplot.xlim(85, 510)
    #matplotlib.pyplot.setp(bp['whiskers'], color='k',  linestyle='-' )
    #matplotlib.pyplot.setp(bp['fliers'], color='k')

    
    ratio=(5,6)
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



def plot_amplitude_candlestick(results, morphology, color='black', picture_file=None, picture_formats=['png', 'pdf', 'svg']):
    matplotlib.rc('text', usetex=True)

    x = []
    y = []
    for frequency in sorted(results.keys()):
        for voltage in results[frequency]:
            x.append(frequency)
            y.append(voltage)
        print frequency, 
        print results[frequency]
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

    #matplotlib.pyplot.title('amplitude in '+str(morphology.replace('_','\_')))
    matplotlib.pyplot.xticks(rotation=90)#  xrange(1+len(data)),titles
    matplotlib.pyplot.grid(True, color='lightgrey')

    bp = matplotlib.pyplot.boxplot(results.values(), positions=results.keys(), widths=10)
    matplotlib.pyplot.ylabel(u'voltage amplitude [mV]')
    matplotlib.pyplot.ylim(0, 40)
    matplotlib.pyplot.xlabel('frequency [Hz]')
    matplotlib.pyplot.xlim(x[0]-10, x[-1]+10)
    #matplotlib.pyplot.setp(bp['whiskers'], color='k',  linestyle='-' )
    #matplotlib.pyplot.setp(bp['fliers'], color='k')

    
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
    return p


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
    picture_formats = ['png','svg', 'pdf']
    c = 'black'
    for swc in sys.argv[1:]:#['../../data/test.swc']:#
        m   = morphjongleur.model.morphology.Morphology.swc_parse(swc, verbose=False)
        print m.name, 
        resultes    = {}
        for f in range(10,501,15):
            amplitues = experiment(morphology=m, compartment=m.root, frequency=f)
            resultes[f] = amplitues
            try:
                plot_distance_amplitude_scatter(
                    [r.compartment.path_distance( r.morphology.root_biggest_child ) for r in results_f],
                    amplitues, 
                    m.name, f, c, '/tmp/distance_amplitude_scatter_%s%iHz' %(m,f), picture_formats
                )
                plot_amplitude_histogramm(amplitues, m.name, f, c, '/tmp/voltage_histogramm_%s%iHz'  %(m,f), picture_formats )
                #plot_duration_histogramm([r.duration for r in results_f], m.name, f, c, '/tmp/duration_histogramm_%s%iHz' %(m,f), picture_formats )
                #plot_delay_histogramm([r.delay for r in results_f], m.name, f, c, '/tmp/delay_histogramm_%s%iHz'    %(m,f), picture_formats )
                #plot_phaseshift_histogramm([r.phase_angle for r in results_f], m.name, f, c, '/tmp/phase_histogramm_%s%iHz'    %(m,f), picture_formats )
            except IndexError, te:
                import traceback
                print traceback.format_exc()
                pass
            except Exception, e:
                import traceback
                print traceback.format_exc()
        print ''
        #plot_frequency_distance_amplitude_scatter3d(fs, ds, aps, m, c, '/tmp/frequency_distance_amplitude_scatter3d_%s'%(m), picture_formats )

        #plot_duration_candlestick(              results, m_name, 'black', '/tmp/duration_candlestick_%s' %(m),      picture_formats )
        #plot_delay_candlestick(                 results, m_name, 'black', '/tmp/delay_candlestick_%s' %(m),         picture_formats )
        plot_amplitude_candlestick(    resultes, m.name, 'black', '/tmp/amplitude_out_candlestick_%s' % (m.name),    picture_formats )
        #print "%s\t%f\t%f\n" % (m_name, u, t)
        #plot_amplitude_s_candlestick(  results, m_name, 'black', '/tmp/amplitude_s_candlestick_%s' % (m),  picture_formats )
        #plot_phaseshift_candlestick(            results, m_name, 'black', '/tmp/phase_candlestick_%s' %(m),         picture_formats )