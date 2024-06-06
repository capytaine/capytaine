# thanks to https://www.geeksforgeeks.org/how-to-check-if-a-given-point-lies-inside-a-polygon/
### RAY CASTING FUNCTIONS AND CLASS

import numpy as np

class Point2D:
    def __init__(self, x, y):
        self.x = x
        self.y = y

class Line2D:
    def __init__(self, p1, p2):
        self.p1 = p1
        self.p2 = p2

def on_line(l1, p):
    # Check whether p is on the line or not
    if (
        p.x <= max(l1.p1.x, l1.p2.x)
        and p.x >= min(l1.p1.x, l1.p2.x)
        and (p.y <= max(l1.p1.y, l1.p2.y) and p.y >= min(l1.p1.y, l1.p2.y))
    ):
        return True
    return False

def direction(a, b, c):
    val = (b.y - a.y) * (c.x - b.x) - (b.x - a.x) * (c.y - b.y)
    if val == 0:
        # Collinear
        return 0
    elif val < 0:
        # Anti-clockwise direction
        return 2
    # Clockwise direction
    return 1

def is_intersect(l1, l2):
    # Four direction for two lines and points of other line
    dir1 = direction(l1.p1, l1.p2, l2.p1)
    dir2 = direction(l1.p1, l1.p2, l2.p2)
    dir3 = direction(l2.p1, l2.p2, l1.p1)
    dir4 = direction(l2.p1, l2.p2, l1.p2)

    # When intersecting
    if dir1 != dir2 and dir3 != dir4:
        return True

    # When p2 of line2 are on the line1
    if dir1 == 0 and on_line(l1, l2.p1):
        return True

    # When p1 of line2 are on the line1
    if dir2 == 0 and on_line(l1, l2.p2):
        return True

    # When p2 of line1 are on the line2
    if dir3 == 0 and on_line(l2, l1.p1):
        return True

    # When p1 of line1 are on the line2
    if dir4 == 0 and on_line(l2, l1.p2):
        return True

    return False

def check_inside(poly, n, p):
    """Check if the point p is inside the n-gon "poly".
    Concave polygon are supported.
    """
    # When polygon has less than 3 edge, it is not polygon
    if n < 3:
        return False

    # Create a point at infinity, y is same as point p
    exline = Line2D(p, Point2D(9999, p.y))
    count = 0
    i = 0
    while True:
        # Forming a line from two consecutive points of poly
        side = Line2D(poly[i], poly[(i + 1) % n])
        if is_intersect(side, exline):
            # If side is intersects ex
            if (direction(side.p1, p, side.p2) == 0):
                return on_line(side, p)
            count += 1

        i = (i + 1) % n
        if i == 0:
            break

    # When count is odd
    return count & 1

####

def split_horizontal(boundary, x_split):
    left_points = boundary[boundary[:,0] < x_split, :]
    right_points = boundary[boundary[:,0] >= x_split, :]
    # try:
    #     onX = boundary[boundary[:,0] == x_split,:]
    return left_points, right_points

def find_center(points):
    Xcm = np.mean(points[:,0])
    Ycm = np.mean(points[:,1])
    return [Xcm, Ycm]

def order_counter_clockwise(points, pCm):
    # pCm = find_center(points)
    list_points = np.zeros((len(points),3))
    list_points[:,0:2] = points
    list_points[:,2] = np.arctan((points[:,1]-pCm[1])/(points[:,0]-pCm[0]))

    list_points = list_points[list_points[:, 2].argsort()]

    return list_points[:,0:2]

def counter_clockwise_boundary(points):
    pCm = find_center(points)
    point_left, point_right = split_horizontal(points, pCm[0])

    new1 = order_counter_clockwise(point_left, pCm)
    new2 = order_counter_clockwise(point_right, pCm)

    return np.vstack((new1,new2))

def is_instance_inside(_check_point, boundary_points):
    polygon = []
    for ii in range(len(boundary_points)):
        polygon.append(Point2D(boundary_points[ii, 0], boundary_points[ii, 1]))
    p = Point2D(_check_point[0], _check_point[1])
    n = len(polygon)

    return check_inside(polygon, n, p)
