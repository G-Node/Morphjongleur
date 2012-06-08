# -*- coding: utf-8 -*-
'''
@author: stransky
TODO: as decorator, hidden, list of types

'''

def _return_dict(self):
    dictionary = {}
    filtered = []
    for v in vars(self.__class__):  #SQLAlchemy:    m._sa_instance_state.dictionary
        if (
            self.__dict__.has_key(v) # may not be mapped from DB
            and
            (not v.startswith('_') or self.show_hidden) #or isinstance(self.__dict__[v], types.MethodType) or isinstance(self.__dict__[v], types.BuiltinMethodType) ) #type __module__ =str, method instance
            ):
            filtered.append(v)
            for allowed_type in self.list_of_types:
                if isinstance(self.__dict__[v], allowed_type): #http://docs.python.org/reference/datamodel.html#objects-values-and-types http://docs.python.org/library/types.html#types.FunctionType
                    dictionary[v] = self.__dict__[v]
                    filtered.remove(v)   # not filtered because of confenience
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

def auto_string(cls):
    import numbers, types
    cls._show_hidden  = False
    cls.list_of_types   =[numbers.Number, types.StringTypes]
    cls._return_dict    = _return_dict
    cls.__str__    = __str__
    cls.variable_table    = variable_table
    return cls