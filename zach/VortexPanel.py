import numpy as np
from numpy import pi, sqrt, cos, sin, tan, sign
from numpy import arctan2 as atan2
from numpy import arctan as atan
from numpy import arcsinh as asinh
from numpy import log as ln
import csv
import os
import scipy.optimize as opt

import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["mathtext.fontset"] = "stix"

import airfoil as af

class driver(object):
    def __init__(self):
        super(driver,self).__init__()

    def add_airfoil(self,airfoil):
        self.myaf = airfoil

    def airfoil_coefficient_matrix(self):
        n = self.myaf.nPts
        # Calculate control points and length of each panel
        l, xc, yc, xnorm, ynorm = self.myaf.panel_geom()

        # Compute the "P matrix"
        P = self.p_matrix()

        # Compute the "A matrix"
        A = np.zeros([n, n])
        for i in range(n - 1):
            A[i, :-1] += ynorm[i] * P[i, :, 1, 0] + xnorm[i] * P[i, :, 0, 0]
            A[i, 1:] += ynorm[i] * P[i, :, 1, 1] + xnorm[i] * P[i, :, 0, 1]

        # Last row of A matrix enforces kutta condition
        A[n - 1, 0] = 1.0
        A[n - 1, n - 1] = 1.0

        return A

    def p_matrix(self,offset = 0.0):
        x = self.myaf.x
        y = self.myaf.y
        n = self.myaf.nPts
        # Calculate control points and length of each panel
        l, xc, yc, xnorm, ynorm = self.myaf.panel_geom()
        xc += xnorm * offset
        yc += ynorm * offset

        # Calculate panel coordinates for each control point
        xi = np.zeros([n - 1, n - 1])
        eta = np.zeros([n - 1, n - 1])
        Phi = np.zeros([n - 1, n - 1])
        Psi = np.zeros([n - 1, n - 1])
        for i in range(n - 1):
            xi[i] = (xc[i] - x[:-1]) * ynorm + (yc[i] - y[:-1]) * (-xnorm)
            eta[i] = -(xc[i] - x[:-1]) * (-xnorm) + (yc[i] - y[:-1]) * ynorm
            Phi[i] = atan2(eta[i] * l, eta[i]**2 + xi[i]**2 - xi[i] * l)
            Psi[i] = 0.5 * ln((xi[i]**2 + eta[i]**2) / ((xi[i] - l)**2 + eta[i]**2))

        # Calculate the panel coefficient matrix
        P = np.zeros([n - 1, n - 1, 2, 2])
        af_array = np.asarray([[(x[1:] - x[:-1]), -(y[1:] - y[:-1])],
                            [(y[1:] - y[:-1]),  (x[1:] - x[:-1])]])
        for i in range(n - 1):
            p_array = np.asarray([[((l - xi[i]) * Phi[i] + eta[i] * Psi[i]), (xi[i] * Phi[i] - eta[i] * Psi[i])],
                                [(eta[i] * Phi[i] - (l - xi[i]) * Psi[i] - l), (-eta[i] * Phi[i] - xi[i] * Psi[i] + l)]])
            for j in range(n - 1):
                P[i, j] = 0.5 / (pi * l[j]**2) * np.dot(af_array[:,:,j], p_array[:,:,j])

        return P

    def aero_coefficients(self, m, M, inc, A):
        x = self.myaf.x
        y = self.myaf.y
        n = self.myaf.nPts

        # Set up an array of angles of attack to evaluate
        alphas = np.linspace(m, M, int((M - m) / inc) + 1)
        aero_coeff = []
        cp_profiles = []

        # Calculate control points and length of each panel
        l, xc, yc, xnorm, ynorm = self.myaf.panel_geom()

        # Get a P matrix with slight offset from the surface (used for Cp calculations)
        offset = 1.0e-3
        P = self.p_matrix(offset)

        # Invert the airfoil coefficient matrix to solve for multiple operating conditions
        A_inv = np.linalg.inv(A)

        vec = np.zeros(n)
        for alpha in alphas:
            rad = alpha * pi / 180.0
            vec[:-1] = ((y[1:] - y[:-1]) * cos(rad) - (x[1:] - x[:-1]) * sin(rad)) / l

            # Calculate vortex strengths at each panel endpoint
            gamma = np.dot(A_inv, vec)

            # Compute aerodynamic coefficients
            cl = 0.0
            cmLE = 0.0
            for i in range(n - 1):
                cl += l[i] * (gamma[i] + gamma[i + 1])
                cmLE += l[i] * ((2.0 * x[i] * gamma[i] + x[i] * gamma[i + 1] + x[i + 1] * gamma[i]
                                 + 2.0 * x[i + 1] * gamma[i + 1]) * cos(rad) +
                                (2.0 * y[i] * gamma[i] + y[i] * gamma[i + 1] + y[i + 1] * gamma[i]
                                 + 2.0 * y[i + 1] * gamma[i + 1]) * sin(rad))

            cmLE /= -3.0
            cmQC = cmLE + 0.25 * cl * cos(rad)

            # Compute the pressure coefficients along the surface of the airfoil
            cp_profile = []
            for i in range(n - 1):
                v = [cos(rad), sin(rad)]
                for j in range(n - 1):
                    v += np.dot(P[i, j], gamma[j:j+2])
                cp_profile.append([xc[i], yc[i], np.linalg.norm(v), 1.0 - np.dot(v, v)])

            aero_coeff.append([alpha, cl, cmLE, cmQC])
            cp_profiles.append([alpha, cp_profile])

        return list(map(list, zip(*aero_coeff))), cp_profiles

    def view(self, filename=None):
        # Calculate geometry parameters
#        theta, x, y, xcamber, ycamber, dydx = geometry(airfoil, nnodes, flap_type, hinge, delta)

        x = self.myaf.x
        y = self.myaf.y
        xcamber = self.myaf.xcamber
        ycamber = self.myaf.ycamber
        dydx = self.myaf.dydx
        # Plot the airfoil geometry
        fig = plt.figure(figsize=(10,5))
        ax = fig.add_subplot(111)
    #    ax.set_title('Geometry and Camber Line of a NACA {} Airfoil ({} nodes)'.format(airfoil, nnodes))
        ax.set_xlabel('$x/c$', size=20)
        ax.set_ylabel('$y/c$', size=20)
        ax.plot(x, y, '-', label = 'Airfoil Geometry', lw=2, color = 'black', fillstyle = 'none')
        ax.plot(xcamber, ycamber, '-', lw=1, label = 'Camber Line', color = 'black')
        ax.set_aspect('equal', 'datalim')
    #    plt.legend()
        ax2=ax.twinx()
        ax2.plot(xcamber, dydx, '.', markersize = 2, label = 'Camber-Line Slope', color = 'black')
        ax2.set_ylabel('$dy_c/dx$', style = 'italic', size=20)
        ax2.set_ylim([-1.0,0.6])
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1+h2, l1+l2, loc=1, prop={'size':16})
    #    plt.legend()
        plt.show()
        if(filename != None):
            filename1 = filename + '.pdf'
            print('Saving Figure: ' + filename)
            fig.savefig(filename1,format='pdf')
            cwd = os.getcwd()
            mystring = '/Applications/Inkscape.app/Contents/Resources/bin/inkscape'
            mystring += ' --file '+cwd+'/'+filename1+' --export-wmf '+cwd+'/'+filename+'.wmf'
            #after you run this, you have to open in Inkscape and save as emf
            os.system(mystring)

    def write_xfoil(self, filename = 'coordinates.dat'):
        x = self.myaf.x
        y = self.myaf.y
        n = self.myaf.nPts

        a = open(filename,'w')
        for i in range(n):
            print("{:25.16E}{:25.16E}".format(x[i],y[i]),file=a)
        a.close()
        print("Xfoil file written: ",filename)


    def fit(self, minAlpha, maxAlpha, incAlpha,iprint=False):
        # Calculate geometry parameters
    #    theta, x, y, xcamber, ycamber, dydx = geometry(airfoil, nnodes, flap_type, hinge, delta)

        # Calculate the inverse of the Airfoil Coefficient Matrix
        A = self.airfoil_coefficient_matrix()

        # Calculate aerodynamic coefficients for each desired angle of attack
        aero_coeff, cp_profiles = self.aero_coefficients(minAlpha, maxAlpha, incAlpha, A)
        aero_coeff = np.asarray(aero_coeff)
        alphas = aero_coeff[0][:]
        CL = aero_coeff[1]
        CmLE = aero_coeff[2]
        rads = alphas[:]*np.pi/180.0

        if(iprint):
            print("   Alpha[deg]               CL                      CmLE")
            for i in range(len(alphas)-1):
                print("{:25.16E}{:25.16E}{:25.16E}".format(alphas[i],CL[i],CmLE[i]))

        aL0 = np.sum(CL*cos(rads))*np.sum(sin(rads)**2) - np.sum(CL*sin(rads))*np.sum(sin(rads)*cos(rads))
        aL0 = atan(aL0/(np.sum(CL*cos(rads))*np.sum(sin(rads)*cos(rads)) - np.sum(CL*sin(rads))*np.sum(cos(rads)**2)))

        CL0a = np.sum(CL*cos(rads))/(np.sum(sin(rads)*cos(rads)) - tan(aL0)*np.sum(cos(rads)**2))

        Cmat = np.zeros([3, 3])
        Cmat[0,0] =   np.sum(sin(2.0*rads)**2)
        Cmat[0,1] =   np.sum(CL*sin(2.0*rads)*cos(rads))
        Cmat[0,2] = - np.sum(CL*sin(2.0*rads)*sin(rads))
        Cmat[1,0] =   Cmat[0,1]
        Cmat[1,1] =   np.sum(CL**2*cos(rads)**2)
        Cmat[1,2] = - np.sum(CL**2*cos(rads)*sin(rads))
        Cmat[2,0] = - Cmat[0,2]
        Cmat[2,1] = - Cmat[1,2]
        Cmat[2,2] = - np.sum(CL**2*sin(rads)**2)

        RHS = np.zeros(3)
        RHS[0] = np.sum(CmLE*sin(2.0*rads))
        RHS[1] = np.sum(CmLE*CL*cos(rads))
        RHS[2] = np.sum(CmLE*CL*sin(rads))

        #Check for symmetric airfoil. If so, Cmat is ill-conditioned, and cannot be solved
        if(abs(aL0)<1.0e-12):
            Cans = np.zeros(3)
            Cans[1] = RHS[1]/Cmat[1,1]
        else:
            C_inv = np.linalg.inv(Cmat)
            Cans = np.dot(C_inv, RHS)

        Cm0a = Cans[0]
        CmN = Cans[1]
        CmA = Cans[2]

        if(iprint):
            print("   aL0                      CL0a                     Cm0a                     CmN                      CmA")
            print("{:25.16E}{:25.16E}{:25.16E}{:25.16E}{:25.16E}".format(aL0,CL0a,Cm0a,CmN,CmA))

        return aL0, CL0a, Cm0a, CmN, CmA

    def eq_deflection(self,naca,xf=0.7,flapType=2,alpha=0.0,targetCL=0.5,iprint=False):

        ans = opt.newton(self.f_one_delta,0.0,fprime=None,args=(naca,xf,flapType,alpha,targetCL),tol=1.0e-6,maxiter=50,fprime2=None)
        if(iprint): print("Equivalent Deflection = ",ans)
        return ans

    def f_one_delta(self,x,naca,xf,flapType,alpha,targetCL):
        self.myaf = af.airfoil(naca=naca,xf=xf,flapType=flapType,delta=x)
        CL,CmLE,dummy = self.single(alpha)
        return (CL-targetCL)

    def single(self,alpha):
        theta = self.myaf.theta
        nnodes = self.myaf.nPts

        # Calculate the inverse of the Airfoil Coefficient Matrix
        A = self.airfoil_coefficient_matrix()

        # Calculate aerodynamic coefficients for each desired angle of attack
        aero_coeff, cp_profiles = self.aero_coefficients(alpha, alpha, 1, A)
        alphas = aero_coeff[0][:]
        CL = aero_coeff[1]
        CmLE = aero_coeff[2]
        CmQC = aero_coeff[3]
        return CL[0], CmLE[0], CmQC[0]

    def sweep(self,minAlpha, maxAlpha, incAlpha):
        theta = self.myaf.theta
        nnodes = self.myaf.nPts
        # Calculate geometry parameters
#        theta, x, y, xcamber, ycamber, dydx = geometry(airfoil, nnodes, flap_type, hinge, delta)

        # Calculate the inverse of the Airfoil Coefficient Matrix
        A = self.airfoil_coefficient_matrix()

        # Calculate aerodynamic coefficients for each desired angle of attack
        aero_coeff, cp_profiles = self.aero_coefficients(minAlpha, maxAlpha, incAlpha, A)
        alphas = aero_coeff[0][:]
        cl = aero_coeff[1]
        cmLE = aero_coeff[2]
        cmQC = aero_coeff[3]

        # Plot the cl and cm(LE) as functions of alpha
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title('Lift and Moment Coefficients for a NACA {}\nas a function of Angle of Attack'.format(airfoil))
        ax.set_xlabel(r'$\alpha$ (deg)')
        ax.set_ylabel(r'$C_L, C_m$')
        ax.plot(alphas, cl, 'o', label = r'$C_L$ (Vortex Panel Method)', color = 'black', fillstyle = 'none')
        ax.plot(alphas, cmLE, 's', label = r'$C_{m_{LE}}$ (Vortex Panel Method)', color = 'black')
        plt.legend(loc = 'upper left')
        plt.show()

        # Generate a CSV file of cl and cm data
        with open('CL_CmLE_{}.csv'.format(self.myaf.naca), 'w') as cl_cmle_file:
            writer = csv.writer(cl_cmle_file, lineterminator="\n")
            for i in range(len(alphas)):
                writer.writerow([alphas[i], cl[i], cmLE[i]])

        output_alphas = [-10, 0, 10]
        for alpha in alphas:
            ind = int(alpha) - minAlpha
            cp = list(map(list, zip(*cp_profiles[ind][1])))

            if alpha in output_alphas:
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.set_title(r'Pressure distribution over a NACA {} at $\alpha$ = {} deg'
                    .format(airfoil, alpha))
                ax.set_xlabel('$x/c$')
                ax.set_ylabel('$C_p$')
                ax.plot(cp[0][:int(nnodes/2)], cp[3][:int(nnodes/2)], '--', label = 'Bottom Surface', color = 'black')
                ax.plot(cp[0][int(nnodes/2):], cp[3][int(nnodes/2):], '-', label = 'Top Surface', color = 'black')
                lims = ax.get_ylim()
                ax.set_ylim([lims[1], lims[0]])
                plt.legend(loc='lower right')
                plt.show()

            # Generate a CSV file of cp data
            with open('Cp_{}_{}.csv'.format(self.myaf.naca, alpha), 'w') as cp_file:
                writer = csv.writer(cp_file, lineterminator="\n")
                for i in range(len(cp[0])):
                    writer.writerow([theta[i], cp[0][i], cp[1][i], cp[2][i], cp[3][i]])
