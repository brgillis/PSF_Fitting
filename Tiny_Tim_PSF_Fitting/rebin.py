#!/usr/bin/env python2
#
# rebins a subsampled psf TimyTim to some fractional pixel 
# and applies charge diffusion kernel
#
# Created by Ole Margraff 2014
# Edited by Bryan Gillis 2014

import numpy as np
from astropy.io import fits as pyf
from fits_functions import read_fits, write_fits, add_comment


def rebin_psf(fits_struct, xy_offset, binsize, stamp_size):
    """ Shift subsampled psf by sub-pixel amount, then rebin it, then cut a postage stamp of a given
        size in the middle of the shifted, binned psf.
    
        Requires: fits_file <HDUList>
                  xy_offset (<int>,<int>) (shift length in x and y)
                  binsize <int> (subsampling factor)
                  stamp_size <int> (size of postage stamp to cut)
                  
        Returns: fits_rebinned <HDUList>
        
        Side-effects: (none)
        
        Except: Will raise an exception if:
                    binsize is larger than the size of the subsampled grid
    
    
        TO-DO: Currently this does no summing over partial pixels. This may need to be implemented
               in the future.
    """
    
    # Check for too-large binsize
    unbinned_x_size = fits_struct[0].header['NAXIS1']
    unbinned_y_size = fits_struct[0].header['NAXIS2']
    if ((binsize > unbinned_x_size) or (binsize > unbinned_y_size)):
        raise Exception("binsize passed to rebin_psf is larger than image size.")

    pad_size_x = np.abs(xy_offset[0])
    pad_size_y = np.abs(xy_offset[1])

    # Zero-pad the data frame to make the shifting a bit easier
    padded_frame = np.zeros([unbinned_x_size+pad_size_x,
                             unbinned_y_size+pad_size_y])

    # Assign values to the padded frame, shifting as we do so
    for i in range(unbinned_x_size):
        for j in range(unbinned_y_size):
            padded_frame[i+(pad_size_x+xy_offset[0])/2,j+(pad_size_y+xy_offset[1])/2] = \
                fits_struct[0].data[i,j]
            
    # Rebin the psf now, by summing up binsize x binsize squares of pixels
    
    new_x_size = int(unbinned_x_size/binsize)
    new_y_size = int(unbinned_y_size/binsize)

    binned_frame = np.zeros([new_x_size,new_y_size])
    for i in range(new_x_size):
        for j in range(new_y_size):
            binned_frame[i,j] = np.sum(padded_frame[i*binsize+(pad_size_x-xy_offset[0])/2:
                                                        (i+1)*binsize+(pad_size_x-xy_offset[0])/2, 
                                                    j*binsize+(pad_size_y-xy_offset[1])/2:
                                                        (j+1)*binsize+(pad_size_y-xy_offset[1])/2])
    
    # Check if we have an even number of pixels.
    # If so, cut off the first column
    if(2*(new_x_size/2) == new_x_size):
        binned_frame = binned_frame[1:,:]
        new_x_size -= 1
    if(2*(new_y_size/2) == new_y_size):
        binned_frame = binned_frame[:,1:]
        new_y_size -= 1
        
    # Cut a postage stamp now
    if(new_x_size > stamp_size):
        size_diff = ( new_x_size - stamp_size ) / 2
        binned_frame = binned_frame[size_diff:new_x_size-size_diff,:]
        new_x_size -= 2*size_diff
    if(new_y_size > stamp_size):
        size_diff = ( new_y_size - stamp_size ) / 2
        binned_frame = binned_frame[:,size_diff:new_y_size-size_diff]
        new_y_size -= 2*size_diff
    
    # Turn this ndarray into a proper fits data structure and set up the header
    fits_rebinned = pyf.HDUList(pyf.PrimaryHDU(binned_frame))
    fits_rebinned[0].header = fits_struct[0].header
    fits_rebinned[0].header['NAXIS1'] = new_x_size
    fits_rebinned[0].header['NAXIS2'] = new_y_size

    add_comment(fits_rebinned,
        "Rebinned with binning factor %i and shifted by effective (%.1f,%.1f)." % 
        (binsize, xy_offset[0]/float(binsize), xy_offset[1]/float(binsize))
        )

    return fits_rebinned

def find_substring(l, substring):
    """ Searches through a list of strings for a substring within one of them, and returns the
        line number where it's first found, or -1 if it isn't found.
        
        Requires: l <list of strings>
                  substring <string> (substring to search for)
                  
        Returns: i <int> (line number in which substring is first found, or -1 if it isn't found)
        
        Side-effects: (none)
    """
    
    for i, s in enumerate(l):
        if substring in s:
            return i
    return -1


def read_kernel_from_fits(fits_struct):
    """ Extract the 3x3 charge diffusion kernel from the header of a TinyTim PSF FITS file.

        We require the kernel to be stored in a 3-line COMMENT with 3 float numbers per line,
        as created by TinyTim.
        
        Requires: fits_struct <HDUList> (fits data structure to search for the kernel in
        
        Returns: kernel_n <ndarray (2D, float)> (charge-diffusion kernel)
        
        Side-effects: (none)
        
        Except: Will raise an exception if:
                    The charge diffusion kernel, as output by tinytim, is not found in the comments of
                        the passed fits data structure.
    """
    
    kernel = []

    fits_comments = fits_struct[0].header['COMMENT']

    i = find_substring(fits_comments, 'following kernel')
    if(i==-1):
        raise Exception("Cannot find charge-diffusion kernel in fits file passed to " +
                        "read_kernel_from_fits")
        
    # kernel parameters are located in the three lines following to that index
    ### the kernel is symmetrical on both axes
    ### not sure whether we have the xy-orientations of FITS and kernel correct, maybe in FITS xy are swapped?
    ### this may cause problem for non-symmetrical kernels
    for j in fits_comments[i+1:i+4]:
        kernel.append([float(x) for x in j.split()])
        
    # Convert to an ndarray        
    kernel_n = np.asarray(kernel)

    return kernel_n


def convolve_with_kernel(frame, kernel):
    """ Convolve charge diffusion kernel with rebinned FITS file.
    
        Requires: frame <HDUList> (rebinned psf data structure)
                  kernel <ndarray (3 x 3, float)>
                  
        Returns: psf_convolved <HDUList> (rebinned, convolved psf. Ready for use!)
        
        Side-effects: (none)
        
        Except: Will raise an exception if:
                    The passed kernel is not 3x3
    
    """

    # Given the small size and sharp edge of the kernel, it's probably faster to convolve in real-space
    frame_convolved = np.empty_like(frame[0].data)
    nx = frame_convolved.shape[0]
    ny = frame_convolved.shape[1]
    
    # Check shape of kernel
    if(np.shape(kernel) != (3,3)):
        raise Exception("convolve_with_kernel is only designed to work with 3x3 kernels. Ensure this\n" +
                        "size of kernel is being passed to it.")
    
    for ixc in xrange(0,nx):
        for iyc in xrange(0,ny):
            frame_convolved[ixc,iyc] = np.sum(np.multiply(frame[0].data[max(0,ixc-1):min(nx,ixc+2),
                                                                        max(0,iyc-1):min(ny,iyc+2)],
                                                          kernel[(1 if ixc==0 else 0):(2 if ixc==nx-1 else 3),
                                                                 (1 if iyc==0 else 0):(2 if iyc==ny-1 else 3)]))

    psf_convolved = pyf.HDUList(pyf.PrimaryHDU(frame_convolved.real))
    psf_convolved[0].header = frame[0].header

    add_comment(psf_convolved, "Convolved with above charge diffusion kernel.")

    return psf_convolved


def get_adjustment(f_infile, subsampling_factor, target_pos):
    """ Uses the recorded barycentre of a subsampled psf file to determine the adjustment needed to
        centre it.
        
        Requires: f_infile <string> (FITS file for subsampled psf)
                  subsampling_factor <int>
                  target_pos (<int>,<int>) (sub-pixel shift for target star)
                  
        Returns: xy_shift (<int>,<int>) (Necessary shift for this psf)
        
        Side-effects: (none)
        
        Except: Will raise an exception if:
                    f_infile cannot be read or doesn't have the necessary header values
        
    """
    
    psf = read_fits(f_infile)
    
    try:
        xp_win_psf = float(psf[0].header['XP_WIN'])-1
        yp_win_psf = float(psf[0].header['YP_WIN'])-1
        xp_cen_psf = int(np.floor((psf[0].header['NAXIS1'])/2))
        yp_cen_psf = int(np.floor((psf[0].header['NAXIS2'])/2))
    except KeyError, e:
        raise Exception("FITS file passed to get_adjustment does not have required header values:\n" +
                        str(e))

    x_shift_psf = (xp_cen_psf - xp_win_psf)
    y_shift_psf = (yp_cen_psf - yp_win_psf)
    
    # Get offset due to extra borders from tinytim in subsampled psfs
    # Might want to investigate actual cause, but for now, this scaling works
    border_shift = int(round(subsampling_factor*0.25))
    
    xy_shift = (int(round((target_pos[0]+x_shift_psf)*subsampling_factor-border_shift)),
                int(round((target_pos[1]+y_shift_psf)*subsampling_factor-border_shift)))
    
    return xy_shift

def rebin_and_shift_psf(f_infile, f_outfile, subsampling_factor, target_pos, stamp_size):
    """ Function to rebin and shift a psf, then cut it to match the size of the target star -
        this should be the function imported into other files.
    
        Requires: f_infile <string> (subsampled psf file)
                  f_outfile <string> (shifted and rebinned psf file to be output)
                  subsampling_factor <int> 
                  target_pos (<int>,<int>) (sub-pixel shift for target star)
                  stamp_size <int> (size of postage stamp to be cut)
                  
        Returns: psf_rad, (Kron radius of psf, measured by sextractor)
                 psf_e1, (e1 of psf, measured by sextractor)
                 psf_e2 (e2 of psf, measured by sextractor)
                 
        Side-effects: Overwrites f_outfile with new rebinned psf file
        
        Except: Will raise an exception if:
                    f_infile cannot be read or isn't properly formatted
                    f_outfile cannot be written to
                    subsampling_factor is larger than the size of the image in f_infile
    """

    psf_orig = read_fits(f_infile)

    xy_shift = get_adjustment(f_infile, subsampling_factor, target_pos)
    print 'Shifting by', xy_shift, '(in subsampled scale)'
    psf_binned = rebin_psf(psf_orig, xy_shift, subsampling_factor, stamp_size)
    psf_convolved = convolve_with_kernel(psf_binned, read_kernel_from_fits(psf_orig))

    write_fits(psf_convolved, f_outfile)
    
    psf_rad = psf_orig[0].header['KRON_RAD']
    psf_e1 = psf_orig[0].header['E1_IMAGE']
    psf_e2 = psf_orig[0].header['E2_IMAGE']
    
    return psf_rad, psf_e1, psf_e2
