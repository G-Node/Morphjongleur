#!/usr/bin/python
# -*- coding: utf-8 -*-
'''

@author: stransky
'''
import math
import numpy
import matplotlib.pyplot
from morphjongleur.model.morphology import MorphologyInfo,Morphology

def euclid_distance(x, y):
    if len(x) != len(y):
        raise Exception, 'different dimension';
    dist = 0;
    for i in xrange(len(x)):
        dist += (x[i] - y[i])**2;
    return math.sqrt(dist);


class Sholl3DAnalysis(MorphologyInfo):

    def __init__(self, morphology, center):
        '''
        '''
        self.morphology = morphology
        self.center     = center
        self.branches   = []

        branch_ranges   = {}
        for key,branch in morphology.branches.iteritems():
            distances   = []
            for compartment in branch:
                d   = center.distance_euclid(compartment)
                distances.append( d + compartment.radius )
                distances.append( d - compartment.radius )
                compartment   = compartment.parent
                branch_ranges[key] = (numpy.min(distances),numpy.max(distances))
        
        #print [(c.compartment_id, branch_ranges[c]) for c in branch_ranges.keys()]
        begins   = sorted([begin for (begin,_) in branch_ranges.values()])
        ends     = sorted([end for (_,end) in branch_ranges.values()])
        self.max_dist   = max(begins[len(begins)-1],ends[len(ends)-1])
        #print begins
        #print ends
        
        self.ds = sorted(list(set([center.distance_euclid(c) for c in morphology.compartments])))
        b   = 0
        e   = 0
        for d in range(len(self.ds)):
            while b < len(begins) and begins[b] < self.ds[d]:
                b += 1
            while e < len(ends) and ends[e] < self.ds[d]:
                e += 1
            self.branches.append(int(b-e))
            
        #print [(d,b) for d,b in zip(self.ds,self.branches)]
        



    @staticmethod
    def branches_in_distance(branchpoints, terminal_tips, center):
        '''
        obsolete
        return [(# branches,distance)], distances
        '''
        distance_sort   = lambda compartment_distance: compartment_distance[1]
        distance_branchpoints    = sorted([(c, center.distance_euclid(c)) for c in branchpoints], key=distance_sort)
        distance_terminaltips    = sorted([(c, center.distance_euclid(c)) for c in terminal_tips], key=distance_sort)

        b = 0
        t = 0
        n = 1   #center
        distances_branches  = [ (0,n) ]
        if len(distance_branchpoints) > 0 and distance_branchpoints[b][1] == 0:
            b = 1
        assert len(distance_branchpoints) == 0 or distance_branchpoints[b][1] != 0
        if distance_terminaltips[t][1] == 0:
            t = 1
        assert distance_terminaltips[t][1] != 0
        while b < len(distance_branchpoints) and t < len(distance_terminaltips):
            if distance_branchpoints[b][1] < distance_terminaltips[t][1]:
                n = n + len(distance_branchpoints[b][0].children) -1
                distances_branches.append((distance_branchpoints[b][1], n))
                b = b + 1
            elif distance_branchpoints[b][1] > distance_terminaltips[t][1]:
                n = n - 1
                distances_branches.append((distance_terminaltips[t][1], n))
                t = t + 1
            else: # ==
                b = b + 1
                t = t + 1
                #print >> sys.stderr, "branchpoint %s and terminatltip %s in equal distance" % (distance_branchpoints[b] , distance_terminaltips[t])
#            if n < 0:
#                raise AttributeError, 'less then 0 branches'

        while b < len(distance_branchpoints):
            assert distance_branchpoints[b][1] <= distance_terminaltips[t-1][1]
            n = n + len(distance_branchpoints[b][0].children) -1
            distances_branches.append((distance_branchpoints[b][1], n))
            b = b + 1
        
        while t < len(distance_terminaltips):
            assert len(distance_branchpoints) == 0 or distance_branchpoints[b-1][1] <= distance_terminaltips[t][1]
            n = n - 1
            distances_branches.append((distance_terminaltips[t][1], n))
            t = t + 1
#            if n < 0:
#                raise AttributeError, 'less then 0 branches'
            
        assert b == len(distance_branchpoints)
        assert t == len(distance_terminaltips)
        return distances_branches

    def number(self, distance=7):
        if distance > self.ds[len(self.ds) -1]:
            return 0
        index   = numpy.searchsorted(self.ds, distance)
        return self.branches[index]

    def normalized(self, distance=7):
        return self.number_of_branches(distance) / ( 4./3. * numpy.pi * distance**3. )


    def numbers(self, distances=[]):
        #assert sorted
        bs  = []
        d    = 0
        for distance in distances:
            while d < len(self.ds) and self.ds[d] < distance:
                d += 1
            if d > 0:
                d -= 1
            if distance > self.ds[-1]:
                bs.append(0)
            else:
                bs.append(self.branches[d])
        assert len(distances) == len(bs)
        return bs

    def normalizeds(self, distances=[]):
        return [branch / ( 4./3. * numpy.pi * distance**3. ) for branch,distance in zip(self.numbers(distances),distances)]
                

    def plot(self, color='black', ratio=None, picture_file=None, picture_formats=['png', 'pdf', 'svg']):
        matplotlib.pyplot.step(self.ds, self.branches, color)

        imaxb   = 0
        maxb    = 0
        for b in range(len(self.branches)):
            if self.branches[b] > maxb:
                imaxb   = b
                maxb    = self.branches[b]

        matplotlib.pyplot.axvline(x=self.ds[imaxb], color='black')
        matplotlib.pyplot.axhline(y=self.branches[imaxb], color='black')
    
        #titles=[]
        #titles.append( "Maximum (%.2f, %i)" % (self.ds[imaxb], self.branches[imaxb]) )
        #matplotlib.pyplot.legend(titles, loc='best')
        print "(%.2f, %i)" % (self.ds[imaxb],self.branches[imaxb])
        matplotlib.pyplot.annotate("(%.2f, %i)" % (self.ds[imaxb],self.branches[imaxb]), xy=(self.ds[imaxb], self.branches[imaxb]))
        #matplotlib.pyplot.annotate("%.2f" % (self.ds[imaxb]), xy=(self.ds[imaxb], 0))
        #matplotlib.pyplot.annotate("%i" % (self.branches[imaxb]), xy=(0, self.branches[imaxb]))
        

        matplotlib.pyplot.ylabel('#branches')
        matplotlib.pyplot.ylim(ymin=0)
        matplotlib.pyplot.xlabel(u'distance [µm]')
        matplotlib.pyplot.xlim(xmin=0)
        
        if ratio != None:#matplotlib.figure.figaspect(arg)
            fig = matplotlib.pyplot.gcf()
            fig.set_size_inches(ratio[0],ratio[1])

        if(picture_file != None):
            for picture_format in picture_formats:
                matplotlib.pyplot.savefig(picture_file+'.'+picture_format, format=picture_format, transparent=True)
        else:
            matplotlib.pyplot.show()
        matplotlib.pyplot.close()


    def plot_normalized(self, color='black', picture_file=None, x_log=True, y_log=True, ratio=None, picture_formats=['png', 'pdf', 'svg']):
        import scipy.optimize
        xs  = numpy.arange(1, self.max_dist, .1*self.max_dist/self.morphology.number_of_compartments )#filter(lambda x: x > 0, self.ds)    
        ys  = self.normalizeds(xs)
        matplotlib.pyplot.plot(xs, ys, color)

        p   = numpy.array([0.5,-2.5])
        polynomial_decay    = lambda x,p: p[0] * x ** p[1]
        err = lambda p, f, x, y: ( f(xs,p) - ys )**2
        fp  = polynomial_decay
        p, success   = scipy.optimize.leastsq(err, p, args=(fp, xs, ys))#, maxfev=10000*len(x)
        matplotlib.pyplot.plot(xs,fp(xs,p), color='black')

        matplotlib.pyplot.ylabel(u'#branches / distance^3 [1/µm^3]')
        if y_log:
            matplotlib.pyplot.yscale('log')
        #matplotlib.pyplot.ylim(ymin=0)
        matplotlib.pyplot.xlabel(u'distance [µm]')
        if x_log:
            matplotlib.pyplot.xscale('log')
        #matplotlib.pyplot.xlim(xmin=0)
    
        titles=[]
        #titles.append( "LOG(S/N)")
        titles.append( "$%.2f \cdot d^{%.2f}$" % (p[0], p[1]) )
        matplotlib.pyplot.legend(titles, loc='best')
        print "(%.2f, %.2f)" % (p[0], p[1])
        
        if ratio != None:#matplotlib.figure.figaspect(arg)
            fig = matplotlib.pyplot.gcf()
            fig.set_size_inches(ratio[0],ratio[1])

        if(picture_file != None):
            for picture_format in picture_formats:
                matplotlib.pyplot.savefig(picture_file+'.'+picture_format, format=picture_format, transparent=True)
        else:
            matplotlib.pyplot.show()
        matplotlib.pyplot.close()


if __name__ == '__main__':
    '''
    Parameter: files, not directories

    cd data/swc_files/
    path_to/metric_analysis.py *.swc
    '''
    import sys
    import morphjongleur.util.parser.swc
    picture_formats = ['png','svg', 'pdf']#
    colors = ['#ffccff', '#00ff00','#008000', '#0000ff','#000080']# H060602DB_10_2(whole).swc H060607DB_10_2(whole).swc H060602VB_10_2(whole).swc H060607VB_10_2(whole).swc
    with_head   = True
    i = 0
    for swc in sys.argv[1:]:#['../../data/test.swc']:#
        color = colors[i % len(colors)]
        i += 1

        morphology  = Morphology.swc_parse(swc, verbose=False)
        s           = Sholl3DAnalysis(morphology, morphology.root)
        s.plot(color=color, ratio=(4,4.5), picture_file='/tmp/%s_sholl3d' % (morphology.name), picture_formats=picture_formats)
        s.plot_normalized(color=color, picture_file='/tmp/%s_sholl3d_normalized' % (morphology.name), x_log=False, y_log=False, ratio=(4,4.5), picture_formats=picture_formats)
        s.plot_normalized(color=color, picture_file='/tmp/%s_sholl3d_normalized_x' % (morphology.name), x_log=True, y_log=False, ratio=(4,4.5), picture_formats=picture_formats)
        s.plot_normalized(color=color, picture_file='/tmp/%s_sholl3d_normalized_y' % (morphology.name), x_log=False, y_log=True, ratio=(4,4.5), picture_formats=picture_formats)
        s.plot_normalized(color=color, picture_file='/tmp/%s_sholl3d_normalized_xy' % (morphology.name), x_log=True, y_log=True,  ratio=(4,4.5), picture_formats=picture_formats)

