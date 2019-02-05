#!/usr/bin/python

import numpy,sys

# magnitude of a vector
def magvect(v):
    magnitude = numpy.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
    return magnitude


# angle between vectors
def angle(v1,v2):
    ang = numpy.arccos(numpy.dot(v1,v2)/(magvect(v1)*magvect(v2)))
    ang = ang*180/numpy.pi
    return ang


# distance between vectors
def distance(start1,v1,start2,v2):
    cross_prod = numpy.cross(v1,v2)
    mx = magvect(cross_prod)
    norm = cross_prod/mx
    diff = start1 - start2
    dist = numpy.fabs(numpy.dot(norm,diff))
    return dist


# least-squares fit of two arrays, x and y, and return first and last points from fit
def lsq(x,y):
    slope,intcpt =  numpy.polyfit(x,y,1,full=False)
    first_point = slope*x[0]+intcpt
    last_point = slope*x[-1]+intcpt
    return first_point,last_point

def calc_angles(sseCalp1,sseCalp2):
    start1 = []
    end1 = []
    vect1 = []
    start2 = []
    end2 = []
    vect2 = []
    sseCalp1 = numpy.array(sseCalp1)
    sseCalp2 = numpy.array(sseCalp2)

    count1 = numpy.arange(numpy.shape(sseCalp1)[0])
    count2= numpy.arange(numpy.shape(sseCalp2)[0])

#    print count1,sseCalp1 

    (s,e) = lsq(count1,sseCalp1[:,0])
    start1.append(s)
    end1.append(e)
    (s,e) = lsq(count1,sseCalp1[:,1])
    start1.append(s)
    end1.append(e)
    (s,e) = lsq(count1,sseCalp1[:,2])
    start1.append(s)
    end1.append(e)

    vect1.append(end1[0] - start1[0])
    vect1.append(end1[1] - start1[1])
    vect1.append(end1[2] - start1[2])

    (s,e) = lsq(count2,sseCalp2[:,0])
    start2.append(s)
    end2.append(e)
    (s,e) = lsq(count2,sseCalp2[:,1])
    start2.append(s)
    end2.append(e)
    (s,e) = lsq(count2,sseCalp2[:,2])
    start2.append(s)
    end2.append(e)

    vect2.append(end2[0] - start2[0])
    vect2.append(end2[1] - start2[1])
    vect2.append(end2[2] - start2[2])

    start1 = numpy.asarray(start1)
    start2 = numpy.asarray(start2)
    vect1 = numpy.asarray(vect1)
    vect2 = numpy.asarray(vect2)

    ang1_2 = angle(vect1,vect2)
    # dist1_2 = distance(start1,vect1,start2,vect2)
    # return (ang1_2,dist1_2)
    return ang1_2

class angleCalc():
    def __init__(self):
        self.sseCalp1 = []
        self.sseCalp2 = []

    def getAngles(self,sseCalp1,sseCalp2):
        # (angle,dist) = calc_angles(sseCalp1,sseCalp2)
        angle = calc_angles(sseCalp1,sseCalp2)
        return angle