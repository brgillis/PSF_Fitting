#!/usr/bin/env python

""" @file plot_stacks.py

    Created 16 Jul 2015

    This script is used to plot star/PSF/residual stacks which are
    generated as a result of running PSFit.

    ---------------------------------------------------------------------

    Copyright (C) 2015  Bryan R. Gillis

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

import click
import numpy as np
import PSF_fit_plotting.magic_values as mv

@click.command()
@click.option("--input-folder", "input_folder", default=".",
              help="Folder where stacks are located.")
@click.option("--output-folder", "output_folder", default="./plots",
              help="Folder to place generated plots in.")
@click.option("--type", "stack_type", default="residual", help="What type of stacks to plot. Allowed are " +
              "'residual', 'stamp', and 'psf'.")

@click.option("--nstack", "number_to_plot", default=9, help="Number of stacks to plot.")
@click.option("--choose", "how_to_choose", default="slice",
              help="How to choose stacks to plot. Allowed are " +
              "'first', 'random', 'slice', 'list'.")
@click.option("--size", "size_in_pixels", default=41,
              help="Size in pixels of central portion to be plotted.")
@click.option("--zmin", default=None, help="Minimum pixel value for colormap.")
@click.option("--zmax", default=None, help="Maximum pixel value for colormap.")
@click.option("--log/--lin", "logscale_z", default=True,
              help="Log- or lin-scale pixel values.")
@click.option("--colormap", "colormap_to_use", default="bw", help="Colormap for plotting.")
@click.option("--slice_offset", default=0, help="If using --choose=slice, initial offset for slice.")
@click.option("--list", "list_filename", default="stacks_to_plot.dat",
              help="If using --choose=list, filename for" +
              "list of stacks to plot.")
def main(**kwargs):
    """ @TODO main docstring
    """
    
    # Determine how to display the stacks
    num_stacks_per_col = int(np.sqrt(number_to_plot))
    num_stacks_per_row = number_to_plot / num_stacks_per_col
    if(num_stacks_per_row*num_stacks_per_col < number_to_plot):
        num_stacks_per_row += 1
        
        # Assert this is enough
        assert num_stacks_per_row*num_stacks_per_col >= number_to_plot
    
    pass

if __name__ == "__main__":
    main()
