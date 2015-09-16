""" @file test_sex_cat_reading.py

    Created 16 Sep 2015

    Tests to see if we can read in a sextractor object catalog

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
import pytest
import os

from psf_testing.sextractor_utility import get_objects_from_cat
from tests.conftest import make_cat_line

@pytest.fixture(scope="module")
def mock_sex_catalog_name():
    return ".test_mock_sex_catalog"

@pytest.fixture(scope="module")
def mock_sex_catalog(request, mock_sex_catalog_name, get_default_sky_objs_list):
    
    with open(mock_sex_catalog_name,"w") as fout:
        fout.write("# Comment line\n")
        fout.write("[1] Header detail line\n")
        fout.write("\n")
        for sky_obj in get_default_sky_objs_list:
            fout.write(make_cat_line(xp=sky_obj.x_pix,
                                     yp=sky_obj.y_pix,
                                     mag=sky_obj.mag,
                                     class_star=sky_obj.class_star) + "\n")

    def tear_down():
        os.remove(mock_sex_catalog_name)
        
    request.addfinalizer(tear_down)
    
    return mock_sex_catalog_name

def test_get_objects_from_cat(mock_sex_catalog, get_default_sky_objs_list):
    objects = get_objects_from_cat(mock_sex_catalog)
    
    assert(len(objects)==len(get_default_sky_objs_list))


