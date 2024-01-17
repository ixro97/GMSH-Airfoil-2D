import math

import numpy as np
from gmshairfoil2d.geometry_def import Airfoil, SplineInterpolation
from gmshairfoil2d.basicGoemetry import Point, Spline
import gmsh

class SurfaceDomain:
    """
    description
    ...

    Attributes
    ----------
    extDomain : float
        External domain around the airfoil
    airfoil : float
        Airfoil
    struct : int
        Type of structured Region
        1: limited length which must be given by extension length
        2: wake refinement extends to the external domain boundary
    """
    def __init__(self, airfoil, extDomain, struct=None):
        self.airfoil = airfoil
        self.extDomain = extDomain
        self.struct = struct
        
        print(self.airfoil)
        print(self.extDomain)
    
    def createStruct(self, offsetValue, extensionLength, transitionLength, conicWakeAngle):
        self.offsetValue = offsetValue
        self.airfoil.offset(offsetValue)
        if self.struct == 2:
            extensionLength
        self.airfoil.extendTE(extensionLength, transitionLength)
        for spline in self.airfoil.splines:
            spline.offsetSpline.generate()
        for spline in self.airfoil.extensionSplines:
            spline.generate()
            
        for spline in self.airfoil.extensionSplines:
            spline.offsetSplines = []       

        # Lines
        angle = conicWakeAngle * (math.pi / 180)
        direction = 1
        refStartPoint = self.airfoil.extensionSplines[-1].points[0]
        refEndPoint = self.airfoil.extensionSplines[-1].points[1]

        for idx in range(2):
            # Calculate offset for transition part
            startPoint = self.airfoil.offsetPoints[self.airfoil.idxTE + idx]
            endPoint = Point(refStartPoint.x, refStartPoint.y + direction*self.offsetValue, refStartPoint.z)
            lastPoint = self.airfoil.offsetPoints[self.airfoil.idxTE + idx - direction]

            tangVector = np.array([startPoint.x - lastPoint.x, startPoint.y - lastPoint.y, 0])
            tangVector = tangVector/np.linalg.norm(tangVector)

            controlPoint = SplineInterpolation.generateControlPoint(startPoint, tangVector, endPoint)
            interpolationPoints = SplineInterpolation.quadBezierInterpolation(startPoint, controlPoint, endPoint)
            spline = Spline(interpolationPoints)
            self.airfoil.extensionSplines[0].offsetSplines.append(spline)
            spline.generate()
            
            # Calculate offset for linear extension part
            startPoint = endPoint
            # Prepare external domain for type
            match self.struct:
                case 1:
                    print(f"struct = 1")
                    endPoint = Point(refEndPoint.x, refEndPoint.y + direction*(self.offsetValue + (refEndPoint.x - refStartPoint.x)*math.tan(angle)), refEndPoint.z)
                case 2:
                    angle1 = math.acos(((math.sqrt(self.extDomain.radius**2 - (math.cos(angle)*startPoint.y - math.cos(angle)*self.extDomain.yc - direction*math.sin(angle)*startPoint.x + direction*math.sin(angle)*self.extDomain.xc)**2)*math.cos(angle) + math.cos(angle)**2*self.extDomain.xc - direction*math.sin(angle)*(startPoint.y - self.extDomain.yc)*math.cos(angle) + direction*math.sin(angle)**2*startPoint.x) - self.extDomain.xc)/self.extDomain.radius)
                    angle2 = math.asin()
                    print(angle1)
                    # self.extDomain.split()
            gmsh.model.occ.synchronize()
            gmsh.fltk.run()
            spline = Spline([startPoint, endPoint])
            self.airfoil.extensionSplines[1].offsetSplines.append(spline)
            spline.generate()
            
            direction = -1
        
        