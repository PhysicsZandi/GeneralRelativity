import sympy as sp
from einsteinpy.symbolic import MetricTensor, ChristoffelSymbols, RiemannCurvatureTensor, RicciTensor, RicciScalar, EinsteinTensor, WeylTensor
from IPython.display import display, Math
sp.init_printing()

import time

class Manifold:
    def __init__(self, metric, coordinates):
        self.metric = MetricTensor(metric.tolist(), coordinates)
        start = time.time()
        self.christoffel = ChristoffelSymbols.from_metric(self.metric)
        self.christoffel.simplify()
        end = time.time()
        print("christoffel: ", end-start)
        start = time.time()
        self.riemann = RiemannCurvatureTensor.from_christoffels(self.christoffel)
        self.riemann.simplify()
        end = time.time()
        print("riemann: ", end-start)
        start = time.time()
        self.ricci = RicciTensor.from_riemann(self.riemann)
        self.ricci.simplify()
        end = time.time()
        print("ricci: ", end-start)
        start = time.time()
        self.scalar = RicciScalar.from_riccitensor(self.ricci)
        self.scalar.simplify()
        end = time.time()
        print("scalar: ", end-start)
        start = time.time()
        self.einstein = EinsteinTensor.from_metric(self.metric)
        self.einstein.simplify()
        end = time.time()
        print("einstein: ", end-start)
        start = time.time()
        self.weyl = WeylTensor.from_metric(self.metric)
        self.weyl.simplify()
        end = time.time()
        print("weyl: ", end-start)

    def print_metric(self):
        for i in range(4):
            for j in range(4):
                if self.metric[i, j] != 0:
                    display(Math(f"g_{{{i}{j}}} = " + sp.latex(self.metric[i, j])))

    def print_christoffel(self):
        for i in range(4):
            for j in range(4):
                for k in range(4):
                    if self.christoffel[i, j, k] != 0:
                        display(Math(f"\Gamma^{{{i}}}_{{{j}{k}}} = " + sp.latex(self.christoffel[i, j, k])))

    def print_riemann(self):
        for i in range(4):
            for j in range(4):
                for k in range(4):
                    for l in range(4):
                        if self.riemann[i, j, k, l] != 0:
                            display(Math(f"R^{{{i}}}_{{{j}{k}{l}}} = " + sp.latex(self.riemann[i, j, k, l])))

    def print_ricci(self):
        for i in range(4):
            for j in range(4):
                if self.ricci[i, j] != 0:
                    display(Math(f"R_{{{i}{j}}} = " + sp.latex(self.ricci[i, j])))

    def print_scalar(self):
        display(Math(f"R = " + sp.latex(self.scalar.simplify())))

    def print_einstein(self):
        for i in range(4):
            for j in range(4):
                if self.einstein[i, j] != 0:
                    display(Math(f"G_{{{i}{j}}} = " + sp.latex(self.einstein[i, j])))

    def print_weyl(self):
        for i in range(4):
            for j in range(4):
                for k in range(4):
                    for l in range(4):
                        if self.weyl[i, j, k, l] != 0:
                            display(Math(f"W^{{{i}}}_{{{j}{k}{l}}} = " + sp.latex(self.weyl[i, j, k, l])))
