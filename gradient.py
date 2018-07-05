# use: python3 gradient.py n k [anim]

import numpy as np
import math
from sympy import cos, sin, acos, gegenbauer, var
from mayavi import mlab
from sympy.solvers import solve
from sympy import Symbol
from sympy import *
import sys
import time

class Gradient:
    def __init__(self, n, k):
        self.points = []
        self.n = int(n)
        self.k = int(k)
        self.plt = None
        self.current_k = 0
        

    def pol(self, t):
        c = ((self.n + self.k) * (self.n + self.k - 1)) / (2 * (2 * self.k - 1))
        d = (self.k + 0.5)
        return c * gegenbauer(self.n - self.k, self.k - 0.5, t) + d * (1 - t * t) * gegenbauer(self.n - self.k - 2,
                                                                                               self.k + 1.5, t)

    def getPoints(self):
        self.points = []
        # polos
        polo_norte = (math.pi, 0)
        polo_sur = (0, 0)
        self.points.append(polo_norte)
        self.points.append(polo_sur)
        if self.k > 0 :
            # soluciones de cos kfi= 0
            phis = []
            for i in range(0, 2 * self.k):
                phis.append(-math.pi / (2 * self.k) + (i * math.pi) / self.k)

            # soluciones del gegenbauer
            print("calculating")
            t = Symbol('t')
            solved = solve(self.pol(t), t,quintics=False,quartics=False,cubics=False)
            # solucion de polinomio n-k = 0
            thetas = [acos(point) for point in solved]

            # generar los puntos cruzando los thetas y los phis
            for phi in phis:
                for theta in thetas:
                    self.points.append([theta, phi])

    def pintarEsfera(self):
        theta, phi = var('theta phi')
        theta, phi = np.linspace(0, 2 * np.pi, 50), np.linspace(0, np.pi, 25)

        THETA, PHI = np.meshgrid(theta, phi)

        X = np.sin(PHI) * np.cos(THETA)
        Y = np.sin(PHI) * np.sin(THETA)
        Z = np.cos(PHI)

        mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(600, 600))
        mlab.clf()
        mlab.mesh(X, Y, Z, color=(0.9, 0.9, 0.9))

    def printPoints(self, color=(1, 0, 0.3)):
        curva_pt = np.array([(sin(theta) * sin(phi), sin(theta) * cos(phi), cos(theta)) for theta, phi in self.points])
        xx = np.array([np.float(pt[0]) for pt in curva_pt])
        yy = np.array([np.float(pt[1]) for pt in curva_pt])
        zz = np.array([np.float(pt[2]) for pt in curva_pt])
        self.plt = mlab.points3d(xx, yy, zz, scale_factor=0.05, color=color)

    def pintarPuntosGradiente(self, color=(1, 0, 0.3),animacion=False):
        self.getPoints();
        self.pintarEsfera()
        self.printPoints(color)

        name = 'points' + str(n) + str(k) + '.png'
        # mlab.savefig(name) (bug)
        if animacion:
            self.anim()
        mlab.show()
    
    #funcion para generar la animacion
    @mlab.animate(delay=100)
    def anim(self):
        msplt = self.plt.mlab_source
        k = self.k
        start_time = time.time()  
        for i in range(self.current_k, self.n+1):
            print(str(i) + "/" + str(self.n))
            # paramos el proceso 2 segundos para mejorar la visualización
            time.sleep(2)
            self.k = i
            self.getPoints()
            print("printing")
            curva_pt = np.array(
                [(sin(theta) * sin(phi), sin(theta) * cos(phi), cos(theta)) for theta, phi in self.points])
            xx = np.array([np.float(pt[0]) for pt in curva_pt])
            yy = np.array([np.float(pt[1]) for pt in curva_pt])
            zz = np.array([np.float(pt[2]) for pt in curva_pt])
            msplt.reset(x=xx, y=yy, z=zz)
            print("--- %s seconds ---" % (time.time() -start_time))
            self.current_k = i
            yield
            
        #mlab.close()


n = sys.argv[1]
k = sys.argv[2]
animacion = False
if (len(sys.argv) > 3):
    animacion = sys.argv[3]
    k = 0

gradiente = Gradient(n, k)
gradiente.pintarPuntosGradiente(animacion=animacion)
