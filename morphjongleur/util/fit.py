# -*- coding: utf-8 -*-
'''
@author: stransky

@see http://en.wikipedia.org/wiki/Curve_fitting
@see http://linuxgazette.net/115/andreasen.html
@see http://www.aps.anl.gov/bcda/dataVis/fit.py.html
'''

import scipy.optimize
import numpy #import scipy.io.array_import

def erf(x,p):
    import scipy.special
    return - 0.5 * scipy.special.erf( p[1] *(x + p[0]) ) + 0.5 #http://de.wikipedia.org/wiki/Fehlerfunktion

def logistic(x,p):
    return 1 / (1 + numpy.exp( p[1] * (x + p[0]) ))

def tanh(x,p):  #-> logistische
    return -0.5 * numpy.tanh( p[1] * (x + p[0]) / 2.0 ) + 0.5

def square(x,p): 
    return - 0.5 * ( p[1] * (x + p[0]) ) / numpy.sqrt( 1 + ( p[1] * (x + p[0]) )**2 ) + 0.5

def arctan(x,p):
    return 0.5 * 2.0 /numpy.pi * numpy.arctan( numpy.pi / 2.0 *( p[1] * (x + p[0]) ) ) + 0.5

def abs(x,p):
    return 0.5 * ( p[1] * (x + p[0]) ) / (1 + numpy.abs( p[1] * (x + p[0]) )) + 0.5

def polynom(x,p):
    return p[1] + p[2] * (x + p[0]) + p[3] * (x + p[0])**2 + p[4] * (x + p[0])**3
    #return p[0] * x**p[1] + p[2]

def log(x,p):
    return p[3] * numpy.log ( p[1] * (x + p[0]) ) + p[2]

def l_log(x,p):    #
    return 1.0  / (1 + numpy.exp( p[3] * numpy.log( p[1] * (x + p[0]) ) ) )

def l_l_p(x,p):    # overfitting
    return p[4] / (1 + numpy.exp( p[3] * numpy.log( p[1] * (x + p[0]) ) ) )

def l_p(x,p):    # total overfitting
    return p[4] / (1 + numpy.exp( p[3] / ( p[1] * x **p[2]  + p[0]) ) )#??? != ( p[1] * (x+ p[0] )**p[2]  ) 


def plot_fit(data = numpy.ndarray((11,2), buffer=numpy.array([[1.0,0.0],[0.9,0.01],[0.8,0.04],[0.7,0.09],[0.6,0.16],[0.5,0.25],[0.4,0.36],[0.3,0.49],[0.2,0.64],[0.1,0.81],[0.0,1.0]]),dtype=float), functions=[erf], picture_file=None):
    def residuals(p, y, x): 
        err = ( y - peval(x,p) )**2
        return err

    x = data[:,0]
    y = data[:,1]
    titles=[]
    ranking = {}
    print "function: err (x-Achsenverschiebung, innere Steigung, y-Achsenverschiebung, ausere Steigung)"

    import matplotlib.pyplot    #    from scipy import gplt
    matplotlib.rc('text', usetex=True) # http://matplotlib.sourceforge.net/users/usetex.html http://stackoverflow.com/questions/5408862/matplotlib-unicode-axis-labels-using-the-cairo-renderer http://matplotlib.sourceforge.net/users/usetex.html
    matplotlib.pyplot.xlim(0, int(numpy.round(x[len(x)-1])))
    matplotlib.pyplot.ylim(0, 1)
    matplotlib.pyplot.xlabel("$\sigma$ of $f_{WN}$")
    matplotlib.pyplot.ylabel("vector strength $r$")
    matplotlib.pyplot.plot(x,y, color='black')#
    #matplotlib.pyplot.scatter(x,y, marker='x')#, color='black'
    titles.append('samples')    # must be last by scatterplot???
    for peval in functions:
        p   = numpy.array([0.5,0.5,0.5,0.5,0.9,0.5])
        plsq= scipy.optimize.leastsq(residuals, p, args=(y, x))
        errs= residuals(plsq[0], y, x)
        avg = numpy.sum( errs ) / len(x)
        mxm = numpy.max( errs )
        mnm = numpy.min( errs )

        ranking[avg]  = (peval, plsq, mxm,mnm)
        print "%20s: avg=%f,max=%f,min=%f\t(" %(peval.__name__, avg,mxm,mnm) , 
        for i in range(len(p)):
            print "\t%f" % (plsq[0][i]) ,
        print "\t)"

    for key in sorted(ranking.keys()):
        (peval, plsq, mxm,mnm)  = ranking[key]
        matplotlib.pyplot.plot(x,peval(x,plsq[0]))
        titles.append( "%s: %f" % (peval.__name__.replace('_','\_'), key) )

    matplotlib.pyplot.legend(titles, loc='best')
    matplotlib.pyplot.grid()
    if(picture_file != None):
        matplotlib.pyplot.savefig(picture_file+'')   # png
        matplotlib.pyplot.savefig(picture_file+'.pdf',format='pdf')
        matplotlib.pyplot.savefig(picture_file+'.svg',format='svg')
    else:
        matplotlib.pyplot.show()
    matplotlib.pyplot.close()

if __name__ == "__main__":
    #plot_fit()

    #import mp.model.synapse
    #plot_fit(mp.model.synapse.plot_r())

    data = numpy.loadtxt( '../../../Data/0-3_10000.dat')
    plot_fit(data=data, functions=[],                                 picture_file='../../../doc/fit_data_0-3_100000')
    plot_fit(data=data, functions=[erf, tanh, square, arctan, abs],   picture_file='../../../doc/fit_sigmoid_0-3_100000')
    plot_fit(data=data, functions=[erf, logistic, l_log, l_l_p, l_p], picture_file='../../../doc/fit_final_0-3_100000')
