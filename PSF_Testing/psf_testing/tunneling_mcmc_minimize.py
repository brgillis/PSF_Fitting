""" @file tunneling_mcmc_minimize.py

    Created 1 Aug 2016

    Function to globally minimize a function through a tunneling MCMC
    implementation.

    ---------------------------------------------------------------------

    Copyright (C) 2016 Bryan R. Gillis

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import numpy as np
import numpy.random as random

def tunneling_mcmc_minimize(func,
                            init_params,
                            param_steps=None,
                            param_mins=None,
                            param_maxes=None,
                            nsteps=100,
                            gamma=0.001,
                            seed=None,
                            *args, **kwargs):
    
    ndim = len(init_params)
    
    if param_steps is None:
        param_steps = np.ones(ndim)
        
    if seed is not None:
        random.seed(seed)
        
    best_in = init_params
    best_out = func(init_params,*args,**kwargs)
    
    current_in = best_in
    current_out = best_out
    current_fstun = 1 - np.exp(gamma*(best_out-current_out))
    
    for _i in range(nsteps):
        
        test_in = current_in + param_steps*random.randn(ndim)
        
        # Check test_in against mins and maxes
        if param_mins is not None:
            test_in = np.where(test_in>param_mins,test_in,param_mins)
        if param_maxes is not None:
            test_in = np.where(test_in<param_maxes,test_in,param_maxes)

        test_out = func(test_in)
        test_fstun = 1 - np.exp(gamma*(best_out-test_out))
        current_fstun = 1 - np.exp(gamma*(best_out-current_out))
        
        delta_fstun = test_fstun - current_fstun 
        
        # Check if we move
        if delta_fstun > 0:
            if random.rand() > np.exp(-0.5*delta_fstun):
                continue # Don't move
        else:
            # Check if this is the best so far
            if test_out < best_out:
                best_in = test_in
                best_out = test_out
        
        current_in = test_in
        current_out = test_out
            
    return best_in, best_out
        