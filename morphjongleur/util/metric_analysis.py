'''
frustrum v.s. cylinder

@author: stransky
'''

def euclid_distance(x, y):
    import math
    if len(x) != len(y):
        raise Exception, "different dimension";
    dist = 0;
    for i in range(len(x)):
        dist += (x[i] - y[i])**2;
    return math.sqrt(dist);


def compartments():
    return float('nan')

def leafs():
    return float('nan')

def branching_points():
    #TODO: assert leafes -1  
    return float('nan')

def path_length():
    '''
    total cell length = sum of all paths (based on center of each compartment)
    '''
    dist = 0;
    for x in compartments:
        for y in compartments:
            if(x[0] < y[0]): #o.B.d.A
                d = euclid_distance(x[1:4],y[1:4]);
                if d > dist:
                    dist = d;
    return dist;

def surface_length_frustum():
    '''
    total cell length running over the surface of an frustrum
    '''
    return float('nan')

def pca_length():
    x = None;
    for r in result:
        c = scipy.array([r['x'], r['y'], r['z']]);
        if(x == None):
            x = c;
        else:
            x = scipy.vstack((x,c));
    
    p = mdp.pca(x);
    x_min   = float("inf");
    x_max   = float("-inf");
    for x,y,z in p:
        if x < x_min:
            x_min = x; 
        if x > x_max:
            x_max = x;
    
    return x_max - x_min

def polyeder():
    '''
    http://python.net/~gherman/convexhull.html
http://www.scipy.org/Cookbook
http://michelanders.blogspot.de/2012/02/3d-convex-hull-in-python.html
http://code.activestate.com/recipes/66527-finding-the-convex-hull-of-a-set-of-2d-points/ 
    '''
    pass


def cylindric_volume():
    return float('nan')

def frustum_volume():
    return float('nan')

def cylindric_lateral_area():
    '''
    without side circles of terminal tips
    '''
    return float('nan')

def frustum_lateral_area():
    '''
    without side circles of terminal tips
    '''
    return float('nan')

def cylindric_surface_area():
    return float('nan')

def frustum_surface_area():
    return float('nan')

def arithmetic_mean_branchpoint_radius():
    '''
    mean cross section area
    '''
    return float('nan')

def geometric_mean_branchpoint_radius():
    '''
    mean cross section area
    '''
    return float('nan')

def harmonic_mean_branchpoint_radius():
    '''
    mean cross section area
    '''
    return float('nan')

def arithmetic_mean_branchpoint_distance():
    '''
    mean cross section area
    '''
    return float('nan')

def geometric_mean_branchpoint_distance():
    '''
    mean cross section area
    '''
    return float('nan')

def harmonic_mean_branchpoint_distance():
    '''
    mean cross section area
    '''
    return float('nan')

def cylindric_arithmetic_mean_cross_section_area():
    '''
    mean cross section area
    '''
    return float('nan')

def frustum_arithmetic_mean_cross_section_area():
    '''
    mean cross section area
    '''
    return float('nan')

    #TODO: name
    #mean abstand
    #pca
    #konvexen polyeder zur not emal in mathesoftware