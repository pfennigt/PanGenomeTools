"""
Tests for general coordinate calculation functionality.
"""

import pytest
from pangenometools.utils import clip_coordinates, calculate_coordinate_boundaries

# def test_coordinate_calculation(setup_test_files):
#     """Test coordinate calculation for extraction."""
#     # Test positive strand
#     left_a, left_b, right_a, right_b, additional_padding = calculate_coordinate_boundaries(
#         100, 500, "+", 50, 100, 100, 50, False, False
#     )
#     assert left_a == 50  # 100 - 50
#     assert left_b == 199  # 100 + 100 - 1
#     assert right_a == 451  # 500 - 50 + 1
#     assert right_b == 600  # 500 + 100
    
#     # Test negative strand
#     left_a, left_b, right_a, right_b, additional_padding = calculate_coordinate_boundaries(
#         100, 500, "-", 50, 100, 100, 50, False, False
#     )
#     assert left_a == 401  # 500 - 100 + 1 (reversed for negative strand)
#     assert left_b == 600  # 500 + 50
#     assert right_a == 50   # 100 - 100
#     assert right_b == 199  # 100 + 50 - 1


# def test_clip_coordinates():
#     """Test coordinate clipping."""
    
#     # Test normal clipping
#     ll, lh, rl, rh = clip_coordinates(50, 200, 400, 600, 1000)
#     assert ll == 50
#     assert lh == 200
#     assert rl == 400
#     assert rh == 600
    
#     # Test clipping to chromosome bounds
#     ll, lh, rl, rh = clip_coordinates(-50, 200, 400, 1500, 1000)
#     assert ll == 1
#     assert lh == 200
#     assert rl == 400
#     assert rh == 1000