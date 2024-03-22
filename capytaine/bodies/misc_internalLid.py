# thanks to https://www.geeksforgeeks.org/how-to-check-if-a-given-point-lies-inside-a-polygon/
### RAY CASTING FUNCTIONS AND CLASS
import numpy as np 

class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y
 
class line:
    def __init__(self, p1, p2):
        self.p1 = p1
        self.p2 = p2
 
def onLine(l1, p):
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
 
def isIntersect(l1, l2):
    # Four direction for two lines and points of other line
    dir1 = direction(l1.p1, l1.p2, l2.p1)
    dir2 = direction(l1.p1, l1.p2, l2.p2)
    dir3 = direction(l2.p1, l2.p2, l1.p1)
    dir4 = direction(l2.p1, l2.p2, l1.p2)
 
    # When intersecting
    if dir1 != dir2 and dir3 != dir4:
        return True
 
    # When p2 of line2 are on the line1
    if dir1 == 0 and onLine(l1, l2.p1):
        return True
 
    # When p1 of line2 are on the line1
    if dir2 == 0 and onLine(l1, l2.p2):
        return True
 
    # When p2 of line1 are on the line2
    if dir3 == 0 and onLine(l2, l1.p1):
        return True
 
    # When p1 of line1 are on the line2
    if dir4 == 0 and onLine(l2, l1.p2):
        return True
 
    return False
 
def checkInside(poly, n, p):
    # When polygon has less than 3 edge, it is not polygon
    if n < 3:
        return False
 
    # Create a point at infinity, y is same as point p
    exline = line(p, Point(9999, p.y))
    count = 0
    i = 0
    while True:
        # Forming a line from two consecutive points of poly
        side = line(poly[i], poly[(i + 1) % n])
        if isIntersect(side, exline):
            # If side is intersects ex
            if (direction(side.p1, p, side.p2) == 0):
                return onLine(side, p)
            count += 1
         
        i = (i + 1) % n
        if i == 0:
            break
 
    # When count is odd
    return count & 1
 
#### 

def splitHorizontal(boundary,xSplit): 
    leftPoints = boundary[boundary[:,0] < xSplit,:]
    rightPoints = boundary[boundary[:,0] >= xSplit,:]
    # try:
    #     onX = boundary[boundary[:,0] == xSplit,:]
    return leftPoints, rightPoints
    
def findCenter(points): 
    Xcm = np.mean(points[:,0])
    Ycm = np.mean(points[:,1])
    return [Xcm, Ycm]

def orderCounterClockWise(points,pCm): 
    # pCm = findCenter(points)
    listPoints = np.zeros((len(points),3)) 
    listPoints[:,0:2] = points
    listPoints[:,2] = np.arctan((points[:,1]-pCm[1])/(points[:,0]-pCm[0]))
    
    listPoints = listPoints[listPoints[:, 2].argsort()]

    return listPoints[:,0:2]

def counterClockwiseBoundary(points): 
    pCm = findCenter(points)
    pointLeft,pointRight = splitHorizontal(points,pCm[0])

    new1 = orderCounterClockWise(pointLeft,pCm)
    new2 = orderCounterClockWise(pointRight,pCm)

    return np.vstack((new1,new2))

def isInstanceInside(checkPoint,boundaryPoints):
    polygon = []
    for ii in range(len(boundaryPoints)): 
        polygon.append(Point(boundaryPoints[ii,0],boundaryPoints[ii,1]))
    p = Point( checkPoint[0], checkPoint[1] )
    n = len(polygon)

    return checkInside(polygon, n, p)


