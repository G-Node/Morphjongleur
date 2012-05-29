'''
Created on 23.03.2012
sudo easy_install quantities sqlalchemy sqlamp
@author: stransky
'''

if __name__ == '__main__':
    import sqlalchemy
    print sqlalchemy.__version__
    import datajongleur
    from datajongleur.beanbags.neuro.pq_based import *
    session = datajongleur.get_session()
    
    iq = InfoQuantity([1,2.5,3], 'ms',
            info={'author': 'Mustermann', 'age':13}
        )
    iq.save()
    
    print iq
    
    uuid = iq.uuid
    iq2 = InfoQuantity.load(uuid)
    
    print iq2

    iq2.info['age'] = 17
    print iq.info
    
    print iq + iq
    print type(iq + iq)

    print iq.amount
    print iq.units
