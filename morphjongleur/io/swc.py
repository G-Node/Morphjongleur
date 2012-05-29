# -*- coding: utf-8 -*-
'''
@author: stransky
'''
import morphjongleur.model.morphology

@staticmethod
def swc_parse_stream(stream, name='', file_origin='', description='', time=None ):
    m = morphjongleur.model.morphology.Morphology(
      name                  = name,
      file_origin           = file_origin,
      description           = '',
      datetime_recording    = time
    )
    assert len(m.compartments) == 0
    for row in stream:
        if len(row) < 1:
            print row
            continue
        if row[0].isdigit():
            (compartment_id, type, x,y,z, radius, compartment_parent_id)  = row #TODO > 7 spalten -> error
            c = morphjongleur.model.morphology.Compartment(compartment_id, compartment_parent_id, radius, x,y,z)
            m.add_compartment(c)
        else:
            m.description  += ''.join(row)+'\n'
            print (''.join(row))  #?str.

    return m

@staticmethod
def swc_parse(file_name):
    import csv
    import os
    import datetime
    reader = csv.reader(open(file_name, 'r'), delimiter=' ', quotechar='"')#csv.Dialect.skipinitialspace
    return morphjongleur.model.morphology.Morphology.swc_parse_stream(
            reader, 
            name            = os.path.splitext( os.path.basename(file_name) )[0], 
            file_origin     = file_name, 
            time            = datetime.datetime.fromtimestamp( os.stat(file_name).st_mtime ).isoformat(' ')
        )

@staticmethod
def swcs_parse(file_names):
    morphologies    = []
    for file_name in file_names:
        morphologies.append( morphjongleur.model.morphology.Morphology.swc_parse(file_name) )
    return morphologies
        

def swc_write_stream(self, stream):
    '''
    write swc
    '''
    for c in self.compartments: #TODO: type
        stream.write("%i %i %f %f %f %f %i\n" % (c.compartment_id, 1, c.x, c.y, c.z, c.radius, c.compartment_parent_id))

def swc_write(self, file):
    '''
    write swc to file
    '''
    stream  = open(file, "w")
    self.swc_write_stream(stream)
    stream.close()

morphjongleur.model.morphology.Morphology.swc_parse_stream   = swc_parse_stream
morphjongleur.model.morphology.Morphology.swc_parse          = swc_parse
morphjongleur.model.morphology.Morphology.swcs_parse         = swcs_parse
morphjongleur.model.morphology.Morphology.swc_write_stream   = swc_write_stream
morphjongleur.model.morphology.Morphology.swc_write          = swc_write

if __name__ == "__main__":
    print __debug__
    m = morphjongleur.model.morphology.Morphology.swc_parse('../../data/test.swc')
    print m
    m.swc_write('/tmp/test.swc')
    m = morphjongleur.model.morphology.Morphology.swc_parse('/tmp/test.swc')
    print m
