import numpy as np
from numpy import pi, sqrt, cos, sin, tan, sign
from numpy import arctan2 as atan2
from numpy import arctan as atan
from numpy import arcsinh as asinh
from numpy import log as ln

import scipy.optimize as opt
from timeit import default_timer as timer

import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["mathtext.fontset"] = "stix"

class airfoil(object):
    def __init__(self, naca='0012', nPts=400, flapType=0, xf=1.0, yf=-100.0, delta=0.0):
        super(airfoil, self).__init__()

        self.naca = naca
        self.nPts = nPts
        self.flapType = flapType
        self.xf = xf
        self.yf = yf
        self.delta = delta

        self.setup()

    def setup(self):
        delta = self.delta*pi/180.0
        n = self.nPts

        # Set up theta evenly spaced (cosine clustering in x)
        theta = np.linspace(-pi, pi, n)

        # Compute x values along camber line
        s = 0.5 * (1.0 - cos(theta))

        # Compute nodal x and y coordinates for a symmetric airfoil
        t = float(self.naca[-2:]) * 0.01
        yt = 5.0 * t * (0.2969 * sqrt(s) - 0.1260 * s - 0.3516 * s**2 + 0.2843 * s**3 - 0.1015 * s**4) * sign(theta)

        # Get the max camber and location of max camber
        m = 0.01 * float(self.naca[0])  # Maximum camber (% chord)
        p = 0.1 * float(self.naca[1])  # Location of maximum chamber (% chord)

        # Compute x and y coordinates for a cambered airfoil
        ycamber = np.zeros(n)
        dydx = np.zeros(n)
        xcamber = np.copy(s)

        xf = self.xf
        yf = self.yf
        #Find vertical hinge point if not given
        if(yf==-100.0):
            if xf < p:
                yf = m / p**2 * (2.0 * p * xf - xf**2)
            else:
                yf = m / (1.0 - p)**2 * (1.0 - 2.0 * p + 2.0 * p * xf - xf**2)

        self.yf = yf
        #print("xf,yf = ",xf,yf)
        #Calculate camber line and slope
        for i in range(n):
            # Leading-edge Geometry
            if s[i] < p:
                ycamber[i] = m / p**2 * (2.0 * p * s[i] - s[i]**2)
                dydx[i] = 2.0 * m / p**2 * (p - s[i])
            # Trailing-edge Geometry
            else:
                ycamber[i] = m / (1.0 - p)**2 * (1.0 - 2.0 * p + 2.0 * p * s[i] - s[i]**2)
                dydx[i] = 2.0 * m / (1.0 - p)**2 * (p - s[i])

            #Flap Settings offset yc and dydx
            if((s[i]>xf) and (self.flapType>0) and (delta != 0.0)):
                if(self.flapType==1): #Traditional Flap
                    r = sqrt((ycamber[i]-yf)**2 + (s[i]-xf)**2)
                    psi = atan((ycamber[i]-yf)/(s[i]-xf))
                    xcamber[i] = xf + r*cos(delta-psi)
                    ycamber[i] = yf - r*sin(delta-psi)
                    dydx[i] = (dydx[i] - tan(delta))/(1+dydx[i]*tan(delta))
                if(self.flapType==2): #Parabolic Flap
                    length = sqrt(yf**2+(1-xf)**2)
                    ghi = -atan(yf/(1-xf))
                    R = sqrt(4*tan(delta)**2 + 1) + asinh(2*tan(delta))/2.0/tan(delta)
#                    if(delta<0.01):
#                        R = 1.0+sqrt(4*delta**2+1.0)
                    xite = 2*length/R
                    etate = -2*length/R*tan(delta)
                    xio = (xcamber[i]-xf)*length/(1-xf)
                    xip = opt.newton(self.f_eq_28,xite*xio/length,fprime=None,args=(xio,length,R,delta),tol=1.0e-12,maxiter=50,fprime2=None)
                    etap = -xip**2/xite*tan(delta)
                    detadxi = -2*xip/xite*tan(delta)
                    xp = xf + xip*cos(ghi) - etap*sin(ghi)
                    yp = yf + xip*sin(ghi) + etap*cos(ghi)
                    yo = yf*(1-xio/length)
                    dyc = ycamber[i] - yo
                    xcamber[i] = xp + dyc*sin(atan(2*xip/xite*tan(delta)))
                    ycamber[i] = yp + dyc*cos(atan(2*xip/xite*tan(delta)))
                    dydx[i] = (dydx[i] - 2*xip*tan(delta)/xite)/(1 + 2*xip*tan(delta)/xite*dydx[i])

        # Add thickness offset to camber location for airfoil surface
        angle = atan(dydx)
        self.x = xcamber - yt * sin(angle)
        self.y = ycamber + yt * cos(angle)
        self.theta = theta
        self.xcamber = xcamber
        self.ycamber = ycamber
        self.dydx = dydx
#        for i in range(n):
#            print(i,self.x[i],self.y[i])

    def f_eq_28(self,x,xio,length,R,delta):
        return -xio + x/2.0*sqrt((x/length*R*tan(delta))**2 + 1) + length*asinh(x*R*tan(delta)/length)/(2*R*tan(delta))

    def panel_geom(self):
        x = self.x
        y = self.y

        # Calculate the length of each panel
        length = sqrt((x[1:] - x[:-1])**2 + (y[1:] - y[:-1])**2)

        # Calculate the normal vector for each panel
        xnorm = -(y[1:] - y[:-1]) / length
        ynorm = (x[1:] - x[:-1]) / length

        # Calculate the centroid (control point) for each panel
        xc = (x[:-1] + x[1:]) / 2.0
        yc = (y[:-1] + y[1:]) / 2.0

        return length, xc, yc, xnorm, ynorm
