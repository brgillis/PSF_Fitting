import numpy as np

def get_2d_multipole_moments(struct, xc=-1, yc=-1, weight_power=0, aperture_size=0):
    """Function to get the multipole moments (through quadrupole)
       of a grid of values.
       
       Requires: struct <HDUList> (fits data structure for image)
       Optional: xc <int> (pixel position of x centre)
                 yc <int> (pixel position of y centre)
                 weight_power <float> (power of r in weighting if negative or zero.
                                       eg. -2 weights by r^-2.
                                       If positive, scale for exponential weighting in pixels:
                                       5 weights by e^(-r (in pix)/5))
                 aperture_size <float> (size of circular aperture to measure moments within)
                 
       Returns: mp <float>,
                dp_x <float>,
                dp_y <float>,
                qp_xx <float>, (using definition: qp_xx = sum(weight*v*(x^2-y^2))
                qp_xy <float>, (using definition: qp_xy = sum(weight*v*2*x*y))
                               (qp_yy and qp_yx are omitted since qp_yx == qp_xy,
                                and in 2D with this definition qp_yy = -qp_xx)
                qp2_xx <float>, (using definition: qp2_xx = sum(weight*v*x^2)
                qp2_yy <float>, (using definition: qp2_yy = sum(weight*v*y^2)
                qp2_xy <float>, (using definition: qp2_xy = sum(weight*v*x*y))
                size <float>, ( = abs(det(qp2)/mp**2)**(1/4) * sign(det(qp2)). Negative values allowed
                               so we don't truncate the noise distribution, but this may need to be
                               altered.)
                e1 <float>, ( = (qp2_xx - qp2_yy) / (qp2_xx + qp2_yy + 2. * sqrt(det(qp2))) )
                e2 <float>  ( = 2. * qp2_xy / (qp2_xx + qp2_yy + 2. * sqrt(det(qp2))) )
                            ( If det(qp2) < 0, e1 = e2 = 0 instead. This too may need to be altered.)
                                
       Side-effects: (none)
       
       Except: Will raise an exception if:
                   values is improperly formatted
                   
    """
    
    values = struct[0].data
    size = np.shape(values)
    
    # Determine length of x/y components of array
    x_length = size[0]
    y_length = size[1]
    
    # Check/determine the centre points
    if(xc<=0):
        xc = (x_length-1)/2
    else:
        xc = int(np.floor(xc))
    if(yc<=0):
        yc = (y_length-1)/2
    else:
        yc = int(np.floor(yc))
        
    # Get the maximum size of a square that can fit around the centre point
    # within the value array. Note that this gives half the side-length
    if(xc<yc):
        sq_size = int(xc)
    else:
        sq_size = int(yc)
    if(sq_size>x_length-xc-1):
        sq_size = int(x_length-xc-1)
    if(sq_size>y_length-yc-1):
        sq_size = int(y_length-yc)
        
    # Initialise the moments to zero
    mp = 0
    dp_x = 0
    dp_y = 0
    qp_xx = 0
    qp_xy = 0
    qp2_xx = 0
    qp2_yy = 0
    
    weight = 0

    # Loop over the values array, summing up each pixel's contribution to the moment
    for ix in xrange(xc-sq_size,xc+sq_size+1):
        dx = ix-xc
        for iy in xrange(yc-sq_size,yc+sq_size+1):
            dy = iy-yc
            d = np.sqrt(dx*dx+dy*dy)
            
            if((d<=aperture_size) or (aperture_size<=0)):
                
                # Determine the weighted value for this point
                if(weight_power==0):
                    v = values[ix,iy]
                    weight += 1
                else:
                    if(d<=0):
                        v = 0
                    else:
                        v = values[ix,iy] * np.power(d,weight_power)
                        weight += np.power(d,weight_power)
    
                mp += v
                
                dp_x += v * dx
                dp_y += v * dy
                
                qp_xx += v * ( 2*dx**2 - d**2 )
                qp_xy += v * ( 2*dx*dy )
                
                qp2_xx += v * dx**2
                qp2_yy += v * dy**2
        
    # Normalise by the total weight
    if(weight > 0):
        mp /= weight
        dp_x /= weight
        dp_y /= weight
        qp_xx /= weight
        qp_xy /= weight
        qp2_xx /= weight
        qp2_yy /= weight
        
    qp2_xy = qp_xy/2
        
    det = qp2_xx*qp2_yy-qp2_xy**2
    absdet = np.abs(det)
    
    if((det==0)or(mp==0)):
        size = 0
        e1 = 0
        e2 = 0
    else:
        size = np.power(absdet/(mp*mp),0.25)*absdet/det
        if(det<0):
            e1 = 0
            e2 = 0
        else:
            e1 = (qp2_xx - qp2_yy) / (qp2_xx + qp2_yy + 2. * np.sqrt(absdet))
            e2 = 2. * qp2_xy / (qp2_xx + qp2_yy + 2. * np.sqrt(absdet))
    
    return mp, dp_x, dp_y, qp_xx, qp_xy, qp2_xx, qp2_yy, qp2_xy, size, e1, e2

def get_size_and_shape(struct, xc=-1, yc=-1, weight_power=0, total_flux=0, aperture_size=0):
    """Function to estimate the size of an image through a weighted quadrupole
       
       Requires: struct <HDUList> (fits data structure for image)
       Optional: xc <int> (pixel position of x centre)
                 yc <int> (pixel position of y centre)
                 weight_power <float> (power of r in weighting if negative or zero.
                                       eg. -2 weights by r^-2.
                                       If positive, scale for exponential weighting in pixels:
                                       5 weights by e^(-r (in pix)/5))
                 total_flux <float> (normalisation factor for total flux of star. If zero, monopole
                                     moment of the image is calculated and used for normalisation)
                 aperture_size <float> (size of circular aperture to measure moments within)
                 
       Returns: size <float>, ( = abs(det(qp2)/mp**2)**(1/4) * sign(det(qp2)). Negative values allowed
                               so we don't truncate the noise distribution, but this may need to be
                               altered.)
                e1 <float>, ( = (qp2_xx - qp2_yy) / (qp2_xx + qp2_yy + 2. * sqrt(det(qp2))) )
                e2 <float>  ( = 2. * qp2_xy / (qp2_xx + qp2_yy + 2. * sqrt(det(qp2))) )
                            ( If det(qp2) < 0, e1 = e2 = 0 instead. This too may need to be altered.)
                                
       Side-effects: (none)
       
       Except: Will raise an exception if:
                   values is improperly formatted
                   
       TO-DO: A lot of code here is duplicated with get_2d_multipole_moments. This function is meant to
              be for the case when only size and shape are needed, so the other moments don't have to
              be calculated. Perhaps a way should be found to minimize code duplication.
    """
    values = struct[0].data
    size = np.shape(values)
    
    # Determine length of x/y components of array
    x_length = size[0]
    y_length = size[1]
    
    # Check/determine the centre points
    if(xc<=0):
        xc = (x_length-1)/2
    else:
        xc = int(np.floor(xc))
    if(yc<=0):
        yc = (y_length-1)/2
    else:
        yc = int(np.floor(yc))
        
    # Get the maximum size of a square that can fit around the centre point
    # within the value array. Note that this gives half the side-length
    if(xc<yc):
        sq_size = int(xc)
    else:
        sq_size = int(yc)
    if(sq_size>x_length-xc-1):
        sq_size = int(x_length-xc-1)
    if(sq_size>y_length-yc-1):
        sq_size = int(y_length-yc)
        
    # Initialise the moments to zero
    qp2_xx = 0
    qp2_yy = 0
    qp2_xy = 0
    mp = 0
    
    weight = 0

    # Loop over the values array, summing up each pixel's contribution to the moment
    for ix in xrange(xc-sq_size,xc+sq_size+1):
        dx = ix-xc
        for iy in xrange(yc-sq_size,yc+sq_size+1):
            dy = iy-yc
            d = np.sqrt(dx*dx+dy*dy)
            
            if((d<=aperture_size) or (aperture_size<=0)):
            
                # Determine the weighted value for this point
                w = get_weight(d,weight_power)
                v = values[ix,iy] * w
                weight += w
                
                mp += v
                qp2_xx += v * dx**2
                qp2_yy += v * dy**2
                qp2_xy += v * dx*dy
        
    # Normalise
    
    if(total_flux != 0):
        weight *= total_flux
    else:
        weight *= mp
        mp *= mp # Since we divide by weight later
    
    if(weight != 0):
        mp /= weight
        qp2_xx /= weight
        qp2_yy /= weight
        qp2_xy /= weight
        
    det = qp2_xx*qp2_yy-qp2_xy**2
    absdet = np.abs(det)
    
    if(det<=0):
        size = 0
        e1 = 0
        e2 = 0
    else:
        size = np.power(absdet,0.25)*absdet/det
        if(det<0):
            e1 = 0
            e2 = 0
        else:
            e1 = (qp2_xx - qp2_yy) / (qp2_xx + qp2_yy + 2. * np.sqrt(absdet))
            e2 = 2. * qp2_xy / (qp2_xx + qp2_yy + 2. * np.sqrt(absdet))
    
    return size, e1, e2, mp

def append_moments( file_name, ID, x_pix, y_pix, star_mag, res_mm, star_size, psf_size,
                 star_e1,star_e2,psf_e1,psf_e2,star_flux ):
    """ Appends a line to the moments file, containing the moments (and other information) passed
        to this function.
        
        Requires: file_name <string> (name of the moments file)
                  ID <int> (ID of this star - designed to be int, but can be anything str-convertible)
                  x_pix <int> (x-pixel coordinate of this star)
                  y_pix <int> (y-pixel coordinate of this star)
                  star_mag <float> (apparent magnitude of this star)
                  res_mm <tuple> (tuple containing residual multipole moments. Should be formatted
                              as it would come from get_2d_multipole_moments)
                  star_size <float> (size measurement for the star)
                  psf_size <float> (size measurement for the psf)
                  star_e1 <float> (e1 for the star)
                  star_e2 <float> (e2 for the star)
                  psf_e1 <float> (e1 for the psf)
                  psf_e2 <float> (e2 for the psf)
                  star_flux <float> (flux for the star)
                  
        Returns: (nothing)
        
        Side-effects: Appends a line to the moments file (file_name) with moments data
        
        Except: Will raise an exception if file_name cannot be accessed and written to.
    
    """
    
    
    with open(file_name, "a") as fo:
        fo.write(str(ID) + "\t" + str(x_pix) + "\t" + str(y_pix) + "\t" +
                   str(res_mm[0]) + "\t" + str(res_mm[1]) + "\t" +
                   str(res_mm[2]) + "\t" + str(res_mm[3]) + "\t" + 
                   str(res_mm[4]) + "\t" + str(star_mag) + "\t" +
                   str(star_size) + "\t" + str(psf_size) + "\t" +
                   str(star_e1) + "\t" + str(star_e2) + "\t" +
                   str(psf_e1) + "\t" + str(psf_e2) + "\t" +
                   str(star_size-psf_size) + "\t" + 
                   str(star_e1-psf_e1) + "\t" + str(star_e2-psf_e2) + "\t" +
                   str(star_flux) + "\n")
        
def get_weight( d, weight_factor ):
    """ Determines a weight value for the multipole moments functions given a distance and the
        passed weight factor.
        
        Requires: d <float> (distance)
                  weight_factor <float> (weight factor passed to the functions)
                  
        Returns: weight <float>
        
        Side-effects: (none)
    
    """ 
    
    if(weight_factor<=0):
        # Power weighting
        if(weight_factor==0):
            return 1
        else:
            if(d<=0):
                return 0
            else:
                return np.power(d,weight_factor)
    else:
        # Exponential weighting
        if(d<=0):
            return 1
        else:
            return np.exp(-d/weight_factor)
        
