#!/usr/bin/env python

""" @file merge_model_focus_tables.py

    Created 13 Aug 2019

    @TODO: File docstring

    ---------------------------------------------------------------------

    Copyright (C) 2019 Bryan R. Gillis

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

import os
import sys
import time

from astropy.table import Table, vstack

import numpy as np


model_focus_table_head = "Focus"
model_focus_table_tail = ".txt"

model_focus_year_start = 2003
model_focus_year_end = 2016

workdir = "/home/brg/Data/HST_Fields"

output_filename = "ModelFocus.fits"

uvis_focus_offset = -0.24

chip1_offset = 0.487
chip2_offset = -0.648


def main(argv):
    """ @TODO main docstring
    """

    num_years = model_focus_year_end - model_focus_year_start + 1
    years = np.linspace(model_focus_year_start, model_focus_year_end, num_years, endpoint=True, dtype=int)

    tables = []

    for year in years:
        filename = os.path.join(workdir, model_focus_table_head + str(year) + model_focus_table_tail)
        table = Table.read(filename, format="ascii.no_header")

        num_rows = len(table["col1"])

        times = np.zeros(num_rows)
        focuses = np.zeros(num_rows)

        for (i, row) in enumerate(table):
            time_string = str(row[1]) + " " + str(row[2]) + " " + str(row[4]) + " " + str(row[3])
            t = time.mktime(time.strptime(time_string, "%b %d %H:%M:%S %Y"))

            times[i] = t
            focuses[i] = row[5] + uvis_focus_offset

        tables.append(Table((times, focuses, focuses + chip1_offset, focuses + chip2_offset),
                            names=('t_sec', 'focus', 'focus_c1', 'focus_c2')))

    merged_table = vstack(tables)

    merged_table.write(os.path.join(workdir, output_filename), overwrite=True)

    return


if __name__ == "__main__":
    main(sys.argv)
