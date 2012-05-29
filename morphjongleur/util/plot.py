# -*- coding: utf-8 -*-
#!/usr/bin/env python
'''
@author: stransky
'''

from mrj.io.database import Database
import sqlalchemy
import numpy
import time
import mrj.orm.experiment_mso
import mrj.model.morphology

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
                d.age = int(d.file_origin[d.file_origin.find('P')+1:d.file_origin.find('P')+3])
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

    def bars_plot(self, datas, columns, x_column='age', bars_column='r', y_column='v_max', x_label='age [d]', y_label='maximal voltage change [mV]', colors=['blue','green','yellow','orange','red'], picture_file=None, picture_formats=['png', 'pdf', 'svg']):
        import matplotlib.pyplot
        start_time = time.time()
        v       = {}
        v_means = {}
        v_stds  = {}
        bars    = sorted( columns[bars_column] )
        xs      = sorted( columns[x_column] )
        for b in bars:
            v[b]        = {}
            for x in xs:
                v[b][x] = []
            v_means[b]  = []
            v_stds[b]   = []

        for d in datas:
            b = d.__getattribute__(bars_column)
            x = d.__getattribute__(x_column)
            y = d.__getattribute__(y_column) 
            y = abs( y - d.neuronal_e )#TODO: e
            v[b][x].append( y )

        print ''
        for x in xs:
            for b in bars:
                print "#%s,%s: %i" % (b,x,len(v[b][x])) ,
                v_means[b].append( numpy.mean( v[b][x] ) )
                v_stds[b].append(  numpy.std(  v[b][x] ) )
            print ''

        ind = numpy.arange( len(xs) )  # the x locations for the groups
        width = 1.0/(len(bars)+1)      # the width of the bars

        end_time = time.time()
        print "preprocess %f" %(end_time-start_time) ,
        start_time = end_time

        matplotlib.pyplot.subplot(111)
        rects   = []
        for i in range(len(bars)):
            rect    = matplotlib.pyplot.bar(ind+i*width, 
                                v_means[ bars[i] ],#TODO: (numpy.mean( v[b][x] ) for x in xs) 
                            width,
                            color=colors[i],
                            yerr=v_stds[ bars[i] ],#(numpy.std( v[b][x] ) for x in xs)
                            #error_kw=dict(elinewidth=6, ecolor='red')
                            )
            rects.append(rect[0])

        # add some
        matplotlib.pyplot.xlabel(x_label)
        matplotlib.pyplot.ylabel(y_label)
        #matplotlib.pyplot.title('titel')
        #matplotlib.pyplot.legend([""], loc='best')
        matplotlib.pyplot.xticks(ind + (len(bars)-1)*width, xs )

        matplotlib.pyplot.legend( tuple(rects), tuple(map(str, bars)), loc='best')
        if(picture_file != None):
            for picture_format in picture_formats:
                matplotlib.pyplot.savefig(picture_file+'.'+picture_format,format=picture_format)
        else:
            matplotlib.pyplot.show()
        matplotlib.pyplot.close()

        print "plot %f" %(time.time()-start_time)

    def dendrite(self, select_morphs=None, quantity='v_max', n=10, description=['exitatory'], picture_file=None, picture_formats=['png', 'pdf', 'svg']):
        start_time = time.time()
        datas   = self.db.session.query(mrj.orm.experiment_mso.MSO_Result_View).filter(
            sqlalchemy.and_(mrj.orm.experiment_mso.MSO_Result_View.r==1
            , mrj.orm.experiment_mso.MSO_Result_View.morphology_key.in_(select_morphs)
            , sqlalchemy.or_(mrj.orm.experiment_mso.MSO_Result_View.axon == 'without', mrj.orm.experiment_mso.MSO_Result_View.morph_name.like('ohne%') )
            , mrj.orm.experiment_mso.MSO_Result_View.n == n if n != None else mrj.orm.experiment_mso.MSO_Result_View.n >= 0
            , mrj.orm.experiment_mso.MSO_Result_View.exp_description.in_(description)
            )
        )
        columns = self.verify( datas )
        assert len(columns['r']) == 1
        assert len(columns['axon']) == 1
        assert n == None or len(columns['n']) == 1
        assert len(columns['synaptic_e']) == 1
        assert len(columns['neuronal_e']) == 1

        self.bars_plot(datas, columns, x_column='age', bars_column='dendrite', y_column=quantity, colors=['black','blue','orange'], picture_file=picture_file, picture_formats=picture_formats)


    def r(self, select_morphs=None, quantity='v_max', n=10, description=['exitatory'], picture_file=None, picture_formats=['png', 'pdf', 'svg']):
        start_time = time.time()
        datas   = self.db.session.query(mrj.orm.experiment_mso.MSO_Result_View).filter(
            sqlalchemy.and_(
            mrj.orm.experiment_mso.MSO_Result_View.n == n if n != None else mrj.orm.experiment_mso.MSO_Result_View.n >= 0
            , mrj.orm.experiment_mso.MSO_Result_View.morphology_key.in_(select_morphs)
            , sqlalchemy.or_(mrj.orm.experiment_mso.MSO_Result_View.axon == 'without', mrj.orm.experiment_mso.MSO_Result_View.morph_name.like('ohne%') )#TODO:in DB
            , mrj.orm.experiment_mso.MSO_Result_View.exp_description.in_(description)
            )
        )
        end_time = time.time()
        print "morphs: db.query %f" %(end_time-start_time) ,

        columns = self.verify( datas )
        assert n == None or len(columns['n']) == 1
        assert len(columns['axon']) == 1
        assert len(columns['synaptic_e']) == 1
        assert len(columns['neuronal_e']) == 1

        self.bars_plot(datas, columns, x_column='age', bars_column='r', y_column=quantity, colors=['blue','green','yellow','orange','red'], picture_file=picture_file, picture_formats=picture_formats)


    def age(self, select_morphs=None, quantity='v_max', n=10, description=['exitatory'], picture_file=None, picture_formats=['png', 'pdf', 'svg']):
        import matplotlib.pyplot
        start_time = time.time()
        datas   = self.db.session.query(mrj.orm.experiment_mso.MSO_Result_View).filter(
            sqlalchemy.and_(mrj.orm.experiment_mso.MSO_Result_View.r==1
            , mrj.orm.experiment_mso.MSO_Result_View.morphology_key.in_(select_morphs)
            , sqlalchemy.or_(mrj.orm.experiment_mso.MSO_Result_View.dendrite == 'both', mrj.orm.experiment_mso.MSO_Result_View.dendrite == 'None')#TODO:in DB
            , sqlalchemy.or_(mrj.orm.experiment_mso.MSO_Result_View.axon == 'without', mrj.orm.experiment_mso.MSO_Result_View.morph_name.like('ohne%') )
            , mrj.orm.experiment_mso.MSO_Result_View.n == n if n != None else mrj.orm.experiment_mso.MSO_Result_View.n >= 0
            , mrj.orm.experiment_mso.MSO_Result_View.exp_description.in_(description)
            )
        )
        end_time = time.time()

        columns = self.verify( datas )
        assert len(columns['r']) == 1
        assert len(columns['dendrite']) == 1
        assert len(columns['axon']) == 1
        assert n == None or len(columns['n']) == 1
        assert len(columns['synaptic_e']) == 1
        assert len(columns['neuronal_e']) == 1

        ages        = {}
        for d in datas:
            if not ages.has_key(d.age):
                ages[d.age]   = []
            ages[d.age].append( abs(d.__getattribute__(quantity) - d.neuronal_e ))#TODO: e

        print ''
        for a in sorted(ages.keys()):
                print "#%s: %i" % (a, len(ages[a])) ,
        print ''

        print "age: db.query %f" %(end_time-start_time) ,
        end_time = time.time()
        print "preprocess %f" %(end_time-start_time) ,
        start_time = end_time

        matplotlib.pyplot.xlabel('age [d]');
        matplotlib.pyplot.ylabel('maximal voltage change [mV]')#quantity.replace('_',' ')
        #matplotlib.pyplot.title('MSO')
        #matplotlib.pyplot.text(60, .025, r'$Beschreibung$')
        #matplotlib.pyplot.grid(True);
        matplotlib.pyplot.boxplot(ages.values(), positions=ages.keys())
        xmin, xmax, ymin, ymax  = matplotlib.pyplot.axis()
        matplotlib.pyplot.xlim(0, 40);
        matplotlib.pyplot.ylim(0, ymax)
        #matplotlib.pyplot.xticks( range(0,42,5) );
        if(picture_file != None):
            for picture_format in picture_formats:
                matplotlib.pyplot.savefig(picture_file+'.'+picture_format,format=picture_format)
        else:
            matplotlib.pyplot.show()
        matplotlib.pyplot.close()

        print "plot %f" %(time.time()-start_time)

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

    def morph(self, datas, quantity='volumes', yaxis_description=None, picture_file=None, picture_formats=['png', 'pdf', 'svg']):
        import matplotlib.pyplot
        start_time = time.time()
        ages        = {}
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
        matplotlib.pyplot.xlim(0, 40);
        matplotlib.pyplot.ylim(0, ymax)
        #matplotlib.pyplot.xticks( range(0,42,5) );
        if(picture_file != None):
            for picture_format in picture_formats:
                matplotlib.pyplot.savefig(picture_file+quantity+'.'+picture_format,format=picture_format)
        else:
            matplotlib.pyplot.show()
        matplotlib.pyplot.close()
        print "plot %f" %(time.time()-start_time)

if __name__ == '__main__':
    p = Plotter()#Agg (for PNGs), PDF, SVG or PS 
    ms          = [500, 501, 505, 506, 507, 512, 513, 514, 515, 516, 517, 520, 521, 522, 523, 524, 525, 526, 527, 528, 529, 530, 531, 532, 534, 535, 536, 537, 542, 543, 544, 545, 546, 547, 548, 549, 550, 551] # range(500, 551+1) # range(392, 443+1) # [244] #
#    p.morphs(  range(500, 551+1), picture_file='../../../doc/morph_')#TODO: compartments #37, rest 38???
#    p.age(     ms, quantity='v_max', n=10, description=['exitatory','exitatory r'], picture_file='../../../doc/age_10')
#    p.age(     ms, quantity='v_max', n=30, description=['exitatory','exitatory r'], picture_file='../../../doc/age_30')
#    p.age(     ms, quantity='v_max', n=50, description=['exitatory','exitatory r'], picture_file='../../../doc/age_50')
    p.age(     ms, quantity='v_max', n=None,  description=['all exitatory'], picture_file='../../../doc/age_all')
#    p.dendrite(ms, quantity='v_max', n=10, description=['exitatory','exitatory r'], picture_file='../../../doc/dendrite_10')
#    p.dendrite(ms, quantity='v_max', n=30, description=['exitatory','exitatory r'], picture_file='../../../doc/dendrite_30')#few data
#    p.dendrite(ms, quantity='v_max', n=50, description=['exitatory','exitatory r'], picture_file='../../../doc/dendrite_50')#few data
    p.dendrite(ms, quantity='v_max', n=None,   description=['all exitatory'], picture_file='../../../doc/dendrite_all')
#    p.r(       ms, quantity='v_max', n=10, description=['exitatory','exitatory r'], picture_file='../../../doc/r_10')
#    p.r(       ms, quantity='v_max', n=30, description=['exitatory','exitatory r'], picture_file='../../../doc/r_30')
#    p.r(       ms, quantity='v_max', n=50, description=['exitatory','exitatory r'], picture_file='../../../doc/r_50')
    p.r(       ms, quantity='v_max', n=None, description=['all exitatory'], picture_file='../../../doc/r_all')
