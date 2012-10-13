# -*- coding: utf-8 -*-
'''
@author: stransky
@see http://www.skymind.com/~ocrow/python_string/
'''

def variable_map(self):
    dictionary = {}
    filtered = []
    for v in vars(self):  #self.__class__    #SQLAlchemy:    m._sa_instance_state.dictionary
        #self.__dict__.has_key(v) and # may not be mapped from DB
        if (not v.startswith('_') or self.__show_hidden): #or isinstance(self.__dict__[v], types.MethodType) or isinstance(self.__dict__[v], types.BuiltinMethodType) ) #type __module__ =str, method instance
            filtered.append(v)
            for allowed_type in self.__list_of_types:
                if isinstance(self.__dict__[v], allowed_type): #http://docs.python.org/reference/datamodel.html#objects-values-and-types http://docs.python.org/library/types.html#types.FunctionType
                    dictionary[v] = self.__dict__[v]
                    filtered.remove(v)   # not filtered because of confenience
    for v in vars(self.__class__):  #self.__class__    #SQLAlchemy:    m._sa_instance_state.dictionary
        if (not v.startswith('_') or self.__show_hidden): #or isinstance(self.__dict__[v], types.MethodType) or isinstance(self.__dict__[v], types.BuiltinMethodType) ) #type __module__ =str, method instance
            if isinstance(self.__class__.__dict__[v], property): #http://docs.python.org/reference/datamodel.html#objects-values-and-types http://docs.python.org/library/types.html#types.FunctionType
                dictionary[v] = getattr(self, v)
            else: 
                filtered.append(v)
    return (dictionary,filtered)

def variable_table(self, keys=None):
    (dictionary, _)    = self.variable_map() 
    if keys == None:
        keys    = sorted(dictionary.keys())
    ks  = "\t".join(["%s" % (str(key))                for key in keys])
    vs  = "\t".join(["%s" % (str(dictionary[key]))    for key in keys])
    return (ks,vs)

def __repr__(self):
    (dictionary, _)    = self.variable_map()
    str_list = ['<', self.__class__.__name__ , '(']
    for key in sorted(dictionary.keys()):
            str_list.append("%s = %s, " % (key, dictionary[key]))
    str_list.append(')>')
    return ''.join(str_list)

def __str__(self):
    (dictionary, _)    = self.variable_map()
    max_key_length  = 0
    for key in dictionary.keys():
        if len(key) > max_key_length:
            max_key_length  = len(key)
    format_string  = "%%%is = %%s\n" % (max_key_length)
    return ''.join([format_string % (key, dictionary[key]) for key in sorted(dictionary.keys())])

def auto_string(cls):
    import numbers, types
    cls.__show_hidden    = False
    cls.__list_of_types   =[numbers.Number, types.StringTypes, types.TupleType, types.DictType, property]
    cls.variable_map    = variable_map
    cls.__repr__        = __repr__
    cls.__str__         = __str__
    cls.variable_table  = variable_table
    return cls
