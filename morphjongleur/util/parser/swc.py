# -*- coding: utf-8 -*-
'''
@author: stransky
'''
import morphjongleur.model.morphology

@classmethod
def swc_parse_stream(cls, stream, name='', file_origin='', description='', time=None, verbose=True ):
    m = cls(
      name                  = name,
      file_origin           = file_origin,
      description           = '',
      datetime_recording    = time
    )
    assert len(m._compartments) == 0
    unparsables  = 0
    errors  = 0
    for row in stream:
        if len(row) < 1:
            continue
        if row[0].startswith('#'):
            m.description  += ''.join(row)+'\n'
            if verbose:
                print (''.join(row))  #?str.
        elif row[0].isdigit():
            try:
                (compartment_id, type, x,y,z, radius, compartment_parent_id)  = row #TODO > 7 spalten -> error
                c = cls._type_compartment(compartment_id, compartment_parent_id, radius, x,y,z)
                m.add_compartment(c)
            except ValueError:
                errors += 1
        else: 
            unparsables += 1
    if (unparsables > 0 or errors > 0) and verbose:
        import sys
        print >> sys.stderr, "%i unparsable lines and %i errors in %s" % (unparsables, errors, m.name)

    return m

@classmethod
def swc_parse(cls, file_name, verbose=True):
    import csv
    import os
    import datetime
    reader = csv.reader(open(file_name, 'r'), delimiter=' ', quotechar='"')#csv.Dialect.skipinitialspace
    return cls.swc_parse_stream(
            reader, 
            name            = os.path.splitext( os.path.basename(file_name) )[0], 
            file_origin     = file_name, 
            time            = datetime.datetime.fromtimestamp( os.stat(file_name).st_mtime ).isoformat(' '),
            verbose         = verbose
        )

@classmethod
def swcs_parse(cls, file_names):
    morphologies    = []
    for file_name in file_names:
        morphologies.append( cls.swc_parse(file_name) )
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
