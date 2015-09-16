""" test_star_selection.py

    Created 16 Sep 2015

    Unit test module for methods and classes related to selecting stars from images.

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

from psf_testing import magic_values as mv
from psf_testing.sky_object import sky_object
from psf_testing.star import star
from psf_testing.star_selection import get_isolated_stars

@pytest.fixture()
def test_x_pix():
    return 35.5

@pytest.fixture()
def test_y_pix():
    return 105.1

@pytest.fixture()
def test_x_size():
    return 18

@pytest.fixture()
def test_y_size():
    return 22

@pytest.fixture()
def test_xp_max(test_x_pix, test_x_size):
    return test_x_pix + test_x_size

@pytest.fixture()
def test_xp_min(test_x_pix, test_x_size):
    return test_x_pix - test_x_size

@pytest.fixture()
def test_yp_max(test_y_pix, test_y_size):
    return test_y_pix + test_y_size

@pytest.fixture()
def test_yp_min(test_y_pix, test_y_size):
    return test_y_pix - test_y_size

@pytest.fixture()
def test_mag():
    return 22.2

@pytest.fixture()
def test_class_star():
    return 0.997

@pytest.fixture()
def make_cat_line(test_x_pix, test_y_pix, test_mag, test_class_star,
                  test_xp_max, test_xp_min, test_yp_max,
                  test_yp_min,
                  xp=None, yp=None, mag=None, class_star=None,
                  xp_max=None, xp_min=None, yp_max=None, yp_min=None):
    
    if(xp is None):
        xp = test_x_pix
    if(yp is None):
        yp = test_y_pix
    if(mag is None):
        mag = test_mag
    if(class_star is None):
        class_star = test_class_star
    if(xp_max is None):
        xp_max = test_xp_max
    if(xp_min is None):
        xp_min = test_xp_min
    if(yp_max is None):
        yp_max = test_yp_max
    if(yp_min is None):
        yp_min = test_yp_min
    
    
    items = []
    for _ in range(mv.sex_cat_num_cols):
        items.append("0")
        
    items[mv.sex_cat_xp_col] = str(xp)
    items[mv.sex_cat_yp_col] = str(yp)
    
    items[mv.sex_cat_mag_col] = str(mag)
    items[mv.sex_cat_class_star_col] = str(class_star)
    
    items[mv.sex_cat_xp_max_col] = str(xp_max)
    items[mv.sex_cat_xp_min_col] = str(xp_min)
    items[mv.sex_cat_yp_max_col] = str(yp_max)
    items[mv.sex_cat_yp_min_col] = str(yp_min)
    
    line = ""
    
    for item in items:
        line = line + item + " "
        
    return line

def test_initialize_sky_obj(make_cat_line,test_x_pix,test_y_pix,test_mag,test_class_star,
                            test_x_size,test_y_size):
    sky_obj = sky_object(make_cat_line)
    
    assert(sky_obj.x_pix==test_x_pix)
    assert(sky_obj.y_pix==test_y_pix)
    assert(sky_obj.mag==test_mag)
    assert(sky_obj.class_star==test_class_star)
    
    assert(sky_obj.stamp_size==max((test_x_size,test_y_size))+1)

@pytest.fixture()
def make_sky_obj(pos=(0.,0.),mag=25.,class_star=0.995):
    obj=sky_object()
    
    obj.x_pix = pos[0]
    obj.y_pix = pos[1]
    
    obj.mag = mag
    obj.class_star = class_star
    
    return obj
    
@pytest.fixture(scope="module")
def make_sky_objs(positions,mags,class_stars):
    objs = []
    for pos, mag, class_star in zip(positions,mags,class_stars):
        objs.append(make_sky_obj(pos, mag, class_star))
    return objs

@pytest.fixture(scope="module")
def make_positions_list():
    positions_list = [(0.,0.),
                      (1.,0.),
                      (2.,0.),
                      (0.5,1),
                      (1.5,1),
                      (-1,-2),
                      (-1,-1),
                      (-2,-1)]
    
    return positions_list
    
@pytest.fixture(scope="module")
def make_mags_list():
    mags_list = [21.5, 22.5, 23.5, 24.5, 25.5, 26.5, 27.5, 20.5]
    
    return mags_list
    
@pytest.fixture(scope="module")
def make_class_stars_list():
    class_stars_list = [0.95, 0.96, 0.97, 0.98, 0.9, 0.945, 0.995, 0.995]
    
    return class_stars_list

@pytest.fixture(scope="module")
def get_default_sky_objs_list(make_positions_list,make_mags_list,make_class_stars_list):
    sky_objs_list = make_sky_objs(make_positions_list,
                                  make_mags_list,
                                  make_class_stars_list)
    return sky_objs_list

def test_convert_sky_obj_to_star(make_sky_obj):
    
    my_star = star(make_sky_obj)
    
    assert(my_star.sky_object.x_pix == make_sky_obj.x_pix)
    assert(my_star.sky_object.y_pix == make_sky_obj.y_pix)
    assert(my_star.sky_object.mag == make_sky_obj.mag)
    assert(my_star.sky_object.class_star == make_sky_obj.class_star)
    
def test_get_lowest_separations(get_default_sky_objs_list):
    
    sky_objs = get_default_sky_objs_list
    
    objects_tree = None
    
    for obj in sky_objs:
        obj.get_lowest_separation(sky_objs,objects_tree)
        
        assert(obj.lowest_separation==1.*mv.pixel_scale)
        
    pass

def test_get_isolated_stars(get_default_sky_objs_list):
    
    stars = []
    
    for obj in get_default_sky_objs_list:
        stars.append(star(obj))
        
    actual_min_sep = 1.*mv.pixel_scale
    
    high_min_sep =  actual_min_sep + 0.01
    low_min_sep =  actual_min_sep - 0.01
    
    high_test_isolated_stars = get_isolated_stars(stars,get_default_sky_objs_list,high_min_sep)
    assert(len(high_test_isolated_stars)==0)
    
    low_test_isolated_stars = get_isolated_stars(stars,get_default_sky_objs_list,low_min_sep)
    assert(len(low_test_isolated_stars)==len(stars))
