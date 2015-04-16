import numpy as np

def get_size(values, xc=-1, yc=-1, aperture_size=-1, weight_func = lambda r : 1):
    
    size = np.shape(values)
    
    # Determine length of x/y components of array
    x_length = size[0]
    y_length = size[1]
    
    # Check/determine the centre points
    if(xc<=0):
        xc = (x_length-1)/2.
    else:
        xc = float(xc)
    if(yc<=0):
        yc = (y_length-1)/2.
    else:
        yc = float(yc)
        
    # Check/determine aperture_size
    if(aperture_size<=0):
        aperture_size = int(max(x_length,y_length)/2)
    else:
        aperture_size = int(aperture_size)
    
    num_rings = aperture_size
    
    num_points = np.zeros(num_rings, dtype=int)
    total_values = np.zeros(num_rings, dtype=float)
    
    
    # Now loop through all points and bin them
    for ix in xrange(0,x_length):
        dx = ix-xc
        for iy in xrange(0,y_length):
            dy = iy-yc
            
            r = np.sqrt(dx*dx+dy*dy)
            
            i = int(r+0.5)
            
            if(i<num_rings):
                num_points[i] += 1
                total_values[i] += values[ix,iy]
    
    # Check if num_points[0] is 0
    
    if(num_points[0]==0):
        num_points[0] = 1 # to prevent a divide-by-zero error
        init_i = 1
    else:
        init_i = 0
    
    mean_values = np.divide(total_values,num_points)
        
    
    # Get the mean contained values
    mean_contained_values = np.zeros_like(mean_values)
    tot = 0.0
    count = 0
    for i in xrange(init_i,num_rings-1):
        tot += num_points[i]*mean_values[i]
        count += num_points[i]
        mean_contained_values[i+1] = tot/count
    
    # Get the dropoff from the difference between mean contained and mean    
    dropoff = np.subtract(mean_contained_values,mean_values)
    
    # And now get the mean distance from zero of the midpoints, using dropoff*w as weight
    tot = 0.0
    tot_weight = 0.0
    for i in xrange(init_i+1,num_rings):
        weight = dropoff[i]*weight_func(i)
        tot += i*weight
        tot_weight += weight
        
    size = tot/tot_weight
    
    return size
    