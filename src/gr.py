from IPython.display import display, Math
import sympy as sp
import time

sp.init_printing()

class Manifold: 
    def __init__(self, metric, coordinates):
        self.coordinates = coordinates
        self.metric = metric
        self.metric = sp.simplify(self.metric)
        self.inverse_metric = metric.inv()
        self.inverse_metric = sp.simplify(self.inverse_metric)

        start = time.time()
        self.compute_christoffel()
        end = time.time()
        print("christoffel: ", end-start)
        start = time.time()
        self.compute_riemann()
        end = time.time()
        print("riemann: ", end-start)
        start = time.time()
        self.compute_ricci()
        end = time.time()
        print("ricci: ", end-start)
        start = time.time()
        self.compute_scalar()
        end = time.time()
        print("scalar: ", end-start)
        start = time.time()
        self.compute_einstein()
        end = time.time()
        print("einstein: ", end-start)
        start = time.time()
        self.compute_kretschmann()
        end = time.time()
        print("kretschmann: ", end-start)

    def compute_christoffel(self):
        self.christoffel = [sp.zeros(4, 4) for _ in range(4)]
        for i in range(4):
            for j in range(4):
                for k in range(4):
                    for l in range(4):
                        self.christoffel[k][i,j] += 1/2 * self.inverse_metric[k,l] * (sp.diff(self.metric[i,l], self.coordinates[j]) + sp.diff(self.metric[j,l], self.coordinates[i]) - sp.diff(self.metric[i,j], self.coordinates[l]))
                    self.christoffel[k][i,j] = sp.simplify(self.christoffel[k][i,j]) 

    def compute_riemann(self):
        self.riemann = [[sp.zeros(4, 4),  sp.zeros(4, 4), sp.zeros(4, 4), sp.zeros(4, 4)], [sp.zeros(4, 4),  sp.zeros(4, 4), sp.zeros(4, 4), sp.zeros(4, 4)], [sp.zeros(4, 4),  sp.zeros(4, 4), sp.zeros(4, 4), sp.zeros(4, 4)], [sp.zeros(4, 4),  sp.zeros(4, 4), sp.zeros(4, 4), sp.zeros(4, 4)]]
        for i in range(4):
            for j in range(4):
                for k in range(4):
                    for l in range(4):
                        self.riemann[i][j][k,l] += sp.diff(self.christoffel[i][j,l], self.coordinates[k]) - sp.diff(self.christoffel[i][j,k], self.coordinates[l])
                        for a in range(4):
                            self.riemann[i][j][k,l] += self.christoffel[a][j,l] * self.christoffel[i][a,k] - self.christoffel[a][j,k] * self.christoffel[i][a,l] 
                        self.riemann[i][j][k,l] = sp.simplify(self.riemann[i][j][k,l])

    def compute_ricci(self):
        self.ricci = sp.zeros(4, 4)
        for i in range(4):
            for j in range(4):
                for k in range(4):
                    self.ricci[i,j] += sp.diff(self.christoffel[k][i,j], self.coordinates[k]) - sp.diff(self.christoffel[k][i,k], self.coordinates[j])
                    for a in range(4):
                        self.ricci[i,j] += self.christoffel[a][i,j] * self.christoffel[k][a,k] - self.christoffel[a][i,k] * self.christoffel[k][a,j]
                self.ricci[i,j] = sp.simplify(self.ricci[i,j])
            
    def compute_scalar(self):
        self.scalar = 0
        for i in range(4):
            self.scalar += self.inverse_metric[i,i] * self.ricci[i,i]
        self.scalar = sp.simplify(self.scalar)

    def compute_einstein(self):
        self.einstein = sp.zeros(4, 4)
        for i in range(4):
            for j in range(4):
                self.einstein[i,j] += self.ricci[i,j] - 0.5 * self.metric[i,j] * self.scalar
                self.einstein[i,j] = sp.simplify(self.einstein[i,j])

    def compute_kretschmann(self):
        self.riemann_down = [[sp.zeros(4, 4),  sp.zeros(4, 4), sp.zeros(4, 4), sp.zeros(4, 4)], [sp.zeros(4, 4),  sp.zeros(4, 4), sp.zeros(4, 4), sp.zeros(4, 4)], [sp.zeros(4, 4),  sp.zeros(4, 4), sp.zeros(4, 4), sp.zeros(4, 4)], [sp.zeros(4, 4),  sp.zeros(4, 4), sp.zeros(4, 4), sp.zeros(4, 4)]]
        self.riemann_up = [[sp.zeros(4, 4),  sp.zeros(4, 4), sp.zeros(4, 4), sp.zeros(4, 4)], [sp.zeros(4, 4),  sp.zeros(4, 4), sp.zeros(4, 4), sp.zeros(4, 4)], [sp.zeros(4, 4),  sp.zeros(4, 4), sp.zeros(4, 4), sp.zeros(4, 4)], [sp.zeros(4, 4),  sp.zeros(4, 4), sp.zeros(4, 4), sp.zeros(4, 4)]]
        for i in range(4):
            for j in range(4):
                for k in range(4):
                    for l in range(4):
                        for a in range(4):
                            self.riemann_down[i][j][k,l] += self.metric[i,a] * self.riemann[a][j][k,l]
                        self.riemann_down[i][j][k,l] = sp.simplify(self.riemann_down[i][j][k,l])
        for i in range(4):
            for j in range(4):
                for k in range(4):
                    for l in range(4):
                        for a in range(4):
                            for b in range(4):
                                for c in range(4):
                                    self.riemann_up[i][j][k,l] += self.inverse_metric[j,a] * self.inverse_metric[k,b] * self.inverse_metric[l,c] * self.riemann[i][a][k,l]
                        self.riemann_up[i][j][k,l] = sp.simplify(self.riemann_up[i][j][k,l])
        self.kretschmann = 0
        for i in range(4):
            for j in range(4):
                for k in range(4):
                    for l in range(4):
                        self.kretschmann += self.riemann_up[i][j][k,l] * self.riemann_down[i][j][k,l]
        self.kretschmann = sp.simplify(self.kretschmann)          

    def print_metric(self):
        for i in range(4):
            for j in range(4):
                if self.metric[i, j] != 0:
                    display(Math(f"g_{{{i}{j}}} = " + sp.latex(self.metric[i, j])))

    def print_inverse_metric(self):
        for i in range(4):
            for j in range(4):
                if self.metric[i, j] != 0:
                    display(Math(f"g^{{{i}{j}}} = " + sp.latex(self.inverse_metric[i, j])))

    def print_christoffel(self):
        for i in range(4):
            for j in range(4):
                for k in range(4):
                    if self.christoffel[i][j,k] != 0:
                        display(Math(f"\Gamma^{{{i}}}_{{{j}{k}}} = " + sp.latex(self.christoffel[i][j,k])))

    def print_riemann(self):
        for i in range(4):
            for j in range(4):
                for k in range(4):
                    for l in range(4):
                        if self.riemann[i][j][k,l] != 0:
                            display(Math(f"R^{{{i}}}_{{{j}{k}{l}}} = " + sp.latex(self.riemann[i][j][k,l])))

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

    def print_kretschmann(self):
        display(Math(f"k = " + sp.latex(self.kretschmann.simplify())))
