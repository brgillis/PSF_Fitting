""" @file star_selection.py

    Created 16 Sep 2015

    Module containing functions to get appropriate stars from a catalog
    generated by sextractor.

    ---------------------------------------------------------------------

    Copyright (C) 2015 user_name

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

from psf_testing.star_selection.sky_object import SkyObj
from psf_testing.star_selection.star import Star

def get_objects_from_cat(sex_cat_filename):
    objects = []

    # Loop through the catalog file and load each line as a sky_object
    with open(sex_cat_filename) as fin:
        # Read in the file, except for comment lines
        for object_line in fin:
            object_line = object_line.strip()
            if len(object_line) == 0:
                continue
            if (object_line[0] != '#') and (object_line[0] != '['):
                objects.append(SkyObj(object_line))

    return objects

def get_stars(objects, min_class_star, min_star_mag, max_star_mag):
    stars = []

    for obj in objects:
        # Test it against conditions, starting with the most likely to fail
        if obj.class_star < min_class_star:
            continue
        if obj.mag > max_star_mag:
            continue
        if obj.mag < min_star_mag:
            continue

        # If we get here, it passed all tests
        stars.append(Star(obj))

    return stars

def get_isolated_stars(stars, all_objects, min_lowest_separation):

    # If there is no min lowest separation, simply return the full stars list
    if min_lowest_separation is None:
        return stars

    # Create a variable for the objects tree we'll share between stars
    object_tree = None

    # Create a list for isolated stars
    isolated_stars = []

    for star in stars:
        # Get the lowest separation for each star
        if star.SkyObj.lowest_separation is None:
            _, object_tree = star.SkyObj.get_lowest_separation(all_objects, object_tree)

        # Select those stars which have a high enough lowest separation
        if star.SkyObj.lowest_separation > min_lowest_separation:
            isolated_stars.append(star)

    return isolated_stars
