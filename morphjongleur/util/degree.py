#!/usr/bin/python
# -*- coding: utf-8 -*-
'''

@author: stransky
'''
import morphjongleur.model.morphology
import morphjongleur.util.metric_analysis
import matplotlib.pyplot


def plot(tuppels=[], colors=[], titles=[], ylabel=None, ratio=None, picture_file=None, picture_formats=['png', 'pdf', 'svg']):
    matplotlib.rc('text', usetex=True)
#    assert len(tuppels) == len(colors)
    for i in range(len(tuppels)):
        xs  = []
        ys  = []
        for x,y in tuppels[i].iteritems():
            xs.append( x )
            ys.append( y )
        matplotlib.pyplot.plot(xs, ys, colors[i%len(colors)])

    matplotlib.pyplot.legend(titles, loc='best')

    matplotlib.pyplot.ylabel(ylabel)
    matplotlib.pyplot.ylim(ymin=0)
    matplotlib.pyplot.xlabel(u'branch\\textdegree')
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



if __name__ == '__main__':
    '''
    Parameter: files, not directories

    cd data/swc_files/
    path_to/metric_analysis.py *.swc
    '''
    import sys
    import morphjongleur.util.parser.swc
    picture_formats = ['png','svg', 'pdf']#
    colors  = ['#00ff00','#008000', '#0000ff','#000080']# H060602DB_10_2(whole).swc H060607DB_10_2(whole).swc H060602VB_10_2(whole).swc H060607VB_10_2(whole).swc
    titles  = ['nurse dorsal branch', 'forager dorsal branch', 'nurse ventral branch', 'forager ventral branch']
    ratio   = (4,4.5)#(9,4.5)
    picture_formats=['png', 'pdf', 'svg']
    with_head   = True
    i = 0
    branches_degrees_mean_lengths       = []
    branches_degrees_mean_volumes       = []
    branches_degrees_mean_lateral_areas = []
    branches_degrees_mean_compactnesss  = []
    for swc in sys.argv[1:]:#['../../data/test.swc']:#
        color = colors[i % len(colors)]
        i += 1

        morphology  = morphjongleur.model.morphology.Morphology.swc_parse(swc, verbose=False)
        a   = morphjongleur.util.metric_analysis.MetricAnalysis(morphology)

        (ks, vs)    = a.variable_table(['name', 'compartments', #'datetime_recording', 
        'total_length', 'slant_length', 
        'number_of_branch_points', 'number_of_terminal_tips', 
        
        'branches_degrees_mean_length', 'branches_degrees_mean_volume', 'branches_degrees_mean_lateral_area', 'branches_degrees_mean_compactness'
        ])
        
        if with_head:
            print ks
            with_head   = False
        print vs
        branches_degrees_mean_lengths.append( a.branches_degrees_mean_length )
        branches_degrees_mean_volumes.append( a.branches_degrees_mean_volume )
        branches_degrees_mean_lateral_areas.append( a.branches_degrees_mean_lateral_area )
        branches_degrees_mean_compactnesss.append( a.branches_degrees_mean_compactness )
        #plot(tuppels=[a.branches_degrees_mean_compactness], colors=colors, titles=titles, ylabel=u'compactness [$\mu m$]', picture_file='/tmp/branches_degrees_mean_compactness_%s' % (morphology.name), picture_formats=picture_formats)
    

    plot(tuppels=branches_degrees_mean_lengths,         colors=colors, titles=titles, ylabel=u'length [$\mu m$]',           ratio=ratio, picture_file='/tmp/branches_degrees_mean_length', picture_formats=picture_formats)
    plot(tuppels=branches_degrees_mean_volumes,         colors=colors, titles=titles, ylabel=u'volume [$\mu m^3$]',         ratio=ratio, picture_file='/tmp/branches_degrees_mean_volume', picture_formats=picture_formats)
    plot(tuppels=branches_degrees_mean_lateral_areas,   colors=colors, titles=titles, ylabel=u'lateral area [$\mu m^2$]',   ratio=ratio, picture_file='/tmp/branches_degrees_mean_lateral_area', picture_formats=picture_formats)
    plot(tuppels=branches_degrees_mean_compactnesss,    colors=colors, titles=titles, ylabel=u'compactness [$\mu m$]',      ratio=ratio, picture_file='/tmp/branches_degrees_mean_compactness', picture_formats=picture_formats)
    
