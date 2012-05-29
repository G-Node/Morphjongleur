# -*- coding: utf-8 -*-
#!/usr/bin/env python
'''
@author: stransky
'''

from mrj.io.database import Database
import mrj.orm.experiment_mso
import mrj.model.morphology
import time

class Plotter(object):
    '''
    http://matplotlib.sourceforge.net/api/pyplot_api.html#matplotlib.pyplot.boxplot
    http://stackoverflow.com/questions/6546602/matplotlib-boxplot-color
    http://www.sqlalchemy.org/docs/orm/tutorial.html#common-filter-operators
    '''
    def __init__(self, db = Database(
                db_name='postgresql://hal08.g-node.pri/lehmann',
                exec_role='morphjokey_admin',
                exec_path='morphjokey',
                echo=True
            ),
            save_only=None
        ):
        self.db = db

        import matplotlib
        if not save_only==None :
            matplotlib.use(save_only) # Force matplotlib to not use any Xwindows backend: http://matplotlib.sourceforge.net/faq/howto_faq.html

        # must be mapped before Object is created
        mapping = mrj.orm.experiment_mso.Mapper( db.engine )
        mapping.orm_map()

    def verify(self, datas, x_column='age', maxlength=100):
        start_time = time.time()
        import numbers, types
        print str(datas)
        columns = {}

        i = 0
        for d in datas:
            i += 1
            if not vars(d).has_key('age') or d.age == None:   #TODO: obsolete
                print "WARNING: morphology has no age info"
                d.age = i
            d.age = int(d.age)
            for c in vars(d):
                if not c.startswith('_') and not columns.has_key(c):
                    columns[c]  = {}
                if isinstance(d.__getattribute__(c), numbers.Number) or isinstance(d.__getattribute__(c), types.StringTypes):
                    columns[c][d.__getattribute__(c)]  = True
                else:
                    pass #list, ... 

                if not columns.has_key(x_column):
                    columns[x_column]  = {}
                columns[x_column][d.__getattribute__(x_column)]   = True

        for c in columns.keys():
            columns[c]  = sorted(columns[c].keys())
            if len(columns[c]) <= maxlength:
                print "%s #%i:\t%s" % (c, len(columns[c]), ", ".join(map(str, columns[c])) )
        print "= %i datalines" % (i)

        print "preprocess %f" %(time.time()-start_time) ,
        return columns

    def morphs(self, select_morphs=None, picture_file=None):
        start_time  = time.time()
        datas   = self.db.session.query(mrj.model.morphology.Morphology_info).filter(
            mrj.model.morphology.Morphology_info.morphology_key.in_(select_morphs)
        )
        end_time = time.time()

        columns = self.verify( datas )
        assert len(columns['morphology_key']) == len(select_morphs)

        print "morphs: db.query %f" %(end_time-start_time)

        p.morph(datas, quantity="cylindric_volume",                                                                   picture_file=picture_file)
        p.morph(datas, quantity="frustum_volume",       yaxis_description='neuron volume [$\mu$m$^3$]',               picture_file=picture_file)
        p.morph(datas, quantity="cylindric_surface_area",                                                             picture_file=picture_file)
        p.morph(datas, quantity="frustum_surface_area", yaxis_description='neuron surface area [$\mu$m$^2$]',         picture_file=picture_file)
        p.morph(datas, quantity="branches",             yaxis_description='\#branches in neuron',                     picture_file=picture_file)
        p.morph(datas, quantity="path_length",          yaxis_description='length of all path in neuron [$\mu$m]',    picture_file=picture_file)
        p.morph(datas, quantity="cylindric_mcse",                                                                     picture_file=picture_file)
        p.morph(datas, quantity="frustum_mcse",         yaxis_description="neuron's mean cross-section area [$\mu$m]",picture_file=picture_file)
        #TODO: p.morph(datas, quantity="spatial_strech", yaxis_description='spatial_strech [$\mu$m]',                  picture_file=picture_file)

    def morph(self, datas, quantity='volumes', yaxis_description=None, picture_file=None, picture_formats=['png']):
        import matplotlib.pyplot
        start_time = time.time()
        ages        = {}
        i = 1
        for d in datas:
            if not vars(d).has_key('age') or d.age == None:   #TODO: should have because of morphs:verify !
                print "WARNING: morphology has no age info"
                d.age = int(d.file_origin[d.file_origin.find('P')+1:d.file_origin.find('P')+3])

            if not ages.has_key(d.age):
                ages[d.age]   = []
            ages[d.age].append( d.__getattribute__(quantity) )

        for a in sorted(ages.keys()):
                print "#%s: %i" % (a, len(ages[a])) ,
        print ''

        end_time = time.time()
        print "preprocess %f" %(end_time-start_time) ,
        start_time = end_time

        matplotlib.rc('text', usetex=True) # http://matplotlib.sourceforge.net/users/usetex.html http://stackoverflow.com/questions/5408862/matplotlib-unicode-axis-labels-using-the-cairo-renderer http://matplotlib.sourceforge.net/users/usetex.html
        matplotlib.pyplot.xlabel('age [d]');
        if yaxis_description == None:
            matplotlib.pyplot.ylabel(quantity.replace('_',' '))
        else:
            matplotlib.pyplot.ylabel(yaxis_description)
        #matplotlib.pyplot.title('MSO')
        #matplotlib.pyplot.text(60, .025, r'$Beschreibung$')
        #matplotlib.pyplot.grid(True);
        matplotlib.pyplot.boxplot(ages.values(),positions=ages.keys())
        xmin, xmax, ymin, ymax  = matplotlib.pyplot.axis()
        matplotlib.pyplot.xlim(0, 6);
        matplotlib.pyplot.ylim(0, ymax)
        #matplotlib.pyplot.xticks( range(0,42,5) );
        if(picture_file != None):
            for picture_format in picture_formats:
                try:
                    matplotlib.pyplot.savefig(picture_file+quantity+'.'+picture_format,format=picture_format)
                except Exception, e:
                    import traceback
                    print picture_format 
                    print traceback.format_exc()
        else:
            matplotlib.pyplot.show()
        matplotlib.pyplot.close()
        print "plot %f" %(time.time()-start_time)

#http://matplotlib.sourceforge.net/faq/howto_faq.html
#http://matplotlib.sourceforge.net/examples/pylab_examples/boxplot_demo2.html
#http://matplotlib.sourceforge.net/api/pyplot_api.html

if __name__ == '__main__':
    p = Plotter(db = Database(
                db_name='postgresql://hal08.g-node.pri/morphjongleur',
                exec_role='morphjokey_admin',
                exec_path='mitsubachi'
            ),
            save_only=None
        )#Agg (for PNGs), PDF, SVG or PS 
    ms          = [3,4,5,6]
    p.morphs(  ms, picture_file='../../../doc/mitsubachi_')#TODO: compartments #37, rest 38???
