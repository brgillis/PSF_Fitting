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

from psf_testing import magic_values as mv
from psf_testing.sky_object import sky_object
from psf_testing.star import star
from psf_testing.star_selection import get_isolated_stars

def test_initialize_sky_obj(test_cat_line,test_x_pix,test_y_pix,test_mag,test_class_star,
                            test_x_size,test_y_size):
    sky_obj = sky_object(test_cat_line)
    
    assert(sky_obj.x_pix==test_x_pix)
    assert(sky_obj.y_pix==test_y_pix)
    assert(sky_obj.mag==test_mag)
    assert(sky_obj.class_star==test_class_star)
    
    assert(sky_obj.stamp_size==max((test_x_size,test_y_size))+1)

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
