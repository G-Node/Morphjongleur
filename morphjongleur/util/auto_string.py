# -*- coding: utf-8 -*-
'''
@author: stransky
'''
class Auto_string(object):
    def _return_dict(self):
        import numbers, types
        dictionary = {}
        filtered = []
        for v in vars(self.__class__):  #SQLAlchemy:    m._sa_instance_state.dictionary
            if (not v.startswith('_') #or isinstance(self.__dict__[v], types.MethodType) or isinstance(self.__dict__[v], types.BuiltinMethodType) ) #type __module__ =str, method instance
                and self.__dict__.has_key(v) ):# may not be mapped from DB
                if isinstance(self.__dict__[v], numbers.Number) or isinstance(self.__dict__[v], types.StringTypes) : #http://docs.python.org/reference/datamodel.html#objects-values-and-types http://docs.python.org/library/types.html#types.FunctionType
                    dictionary[v] = self.__dict__[v]
                else:   # filtered because of confenience
                    filtered.append(v)
        return (dictionary,filtered)
    def __str__(self):
        (dictionary, _)    = self._return_dict()
        s   = ''
        #s   = '<'
        #s  += self.__class__.__name__ 
        #s  += '(\n'
        for key,value in dictionary.iteritems():
                s += "%23s = %s\n" % (key, value)
        #s  += '\n)>'
        return s

    def variable_table(self, keys):
        (dictionary, _)    = self._return_dict() 
        ks  = ''
        vs  = ''
        if keys == None:
            keys    = dictionary.keys()
        for key in keys:
            ks += "%s\t" % (key)
            vs += "%s\t" % (dictionary[key])
        return (ks,vs)
