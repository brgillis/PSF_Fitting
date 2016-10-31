""" @file combine_results.py

    Created 31 Oct 2016

    Operations related to file system paths.

    ---------------------------------------------------------------------

    Copyright (C) 2015 Bryan R. Gillis

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

from os.path import join, isfile

def find_file_in_path(filename, path):
    """
        Searches through a colon-separated path for a file and returns the qualified name of it if found,
        None otherwise.
    """

    colon_separated_path = path.split(":")

    qualified_filename = None

    for test_path in colon_separated_path:

        test_filename = join(test_path, filename)

        if isfile(test_filename):
            qualified_filename = test_filename
            break

    return qualified_filename

def first_in_path(path):
    """
        Gets the first directory listed in the path.
    """

    return path.split(":")[0]
