import sympy as sp
from IPython.display import display, Math
sp.init_printing()
    
class Manifold:
    def __init__(self, metric, coordinates):
        self.coordinates = coordinates
        self.metric = metric 
        self.inverse_metric = self.metric.inv()
        self.compute_christoffel()
        self.compute_riemann()
        self.compute_ricci()
        self.compute_scalar()
        self.compute_einstein()
        self.compute_kretschmann()

        self.check_metric()
        self.check_ricci()
        self.check_einstein()

    def print_metric(self):
        latex_text = f"g = {sp.latex(self.metric)}"
        display(Math(latex_text))

    def print_inverse_metric(self):
        latex_text = f"g^{{-1}} = {sp.latex(self.inverse_metric)}"
        display(Math(latex_text))

    def compute_christoffel(self):
        self.christoffel = [sp.zeros(4, 4) for _ in range(4)]
        for i in range(4):
            for j in range(4):
                for k in range(4):
                    for l in range(4):
                        self.christoffel[k][i,j] += 1/2 * self.inverse_metric[k,l] * (sp.diff(self.metric[i,l], self.coordinates[j]) + sp.diff(self.metric[j,l], self.coordinates[i])-sp.diff(self.metric[i,j], self.coordinates[l]))
                    self.christoffel[k][i,j] = sp.simplify(self.christoffel[k][i,j]) 
        
    def print_christoffel(self):
        for i in range(4):
            for j in range(4):
                for k in range(4):
                    latex_text = f"\\Gamma^{{{k}}}_{{{i}{j}}} = {sp.latex(self.christoffel[k][i, j])}"
                    display(Math(latex_text))

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

    def print_riemann(self):
        for i in range(4):
            for j in range(4):
                for k in range(4):
                    for l in range(4):
                        latex_text = f"R^{{{i}}}_{{{j}{k}{l}}} = {sp.latex(self.riemann[i][j][k, l])}"
                        display(Math(latex_text))

    def compute_ricci(self):
        self.ricci = sp.zeros(4, 4)
        for i in range(4):
            for j in range(4):
                for k in range(4):
                    self.ricci[i,j] += sp.diff(self.christoffel[k][i,j], self.coordinates[k]) - sp.diff(self.christoffel[k][i,k], self.coordinates[j])
                    for a in range(4):
                        self.ricci[i,j] += self.christoffel[a][i,j] * self.christoffel[k][a,k] - self.christoffel[a][i,k] * self.christoffel[k][a,j]
                self.ricci[i,j] = sp.simplify(self.ricci[i,j])

    def print_ricci(self):
        for i in range(4):
            for j in range(4):
                latex_text = f"R_{{{i}{j}}} = {sp.latex(self.ricci[i,j])}"
                display(Math(latex_text))
            
    def compute_scalar(self):
        self.scalar = 0
        for i in range(4):
            self.scalar += self.inverse_metric[i,i] * self.ricci[i,i]
        self.scalar = sp.simplify(self.scalar)

    def print_scalar(self):
        latex_text = f"R = {sp.latex(self.scalar)}"
        display(Math(latex_text))

    def compute_einstein(self):
        self.einstein = sp.zeros(4, 4)
        for i in range(4):
            for j in range(4):
                self.einstein[i,j] += self.ricci[i,j] - 0.5 * self.metric[i,j] * self.scalar
                self.einstein[i,j] = sp.simplify(self.einstein[i,j])

    def print_einstein(self):
        for i in range(4):
            for j in range(4):
                latex_text = f"G_{{{i}{j}}} = {sp.latex(self.einstein[i,j])}"
                display(Math(latex_text))

    def compute_kretschmann(self):
        self.riemann_down = [[sp.zeros(4, 4),  sp.zeros(4, 4), sp.zeros(4, 4), sp.zeros(4, 4)], [sp.zeros(4, 4),  sp.zeros(4, 4), sp.zeros(4, 4), sp.zeros(4, 4)], [sp.zeros(4, 4),  sp.zeros(4, 4), sp.zeros(4, 4), sp.zeros(4, 4)], [sp.zeros(4, 4),  sp.zeros(4, 4), sp.zeros(4, 4), sp.zeros(4, 4)]]
        self.riemann_up = [[sp.zeros(4, 4),  sp.zeros(4, 4), sp.zeros(4, 4), sp.zeros(4, 4)], [sp.zeros(4, 4),  sp.zeros(4, 4), sp.zeros(4, 4), sp.zeros(4, 4)], [sp.zeros(4, 4),  sp.zeros(4, 4), sp.zeros(4, 4), sp.zeros(4, 4)], [sp.zeros(4, 4),  sp.zeros(4, 4), sp.zeros(4, 4), sp.zeros(4, 4)]]
        for i in range(4):
            for j in range(4):
                for k in range(4):
                    for l in range(4):
                        for a in range(4):
                            self.riemann_down[i][j][k,l] += self.metric[i,a] * self.riemann[a][j][k,l]
        for i in range(4):
            for j in range(4):
                for k in range(4):
                    for l in range(4):
                        for a in range(4):
                            for b in range(4):
                                for c in range(4):
                                    self.riemann_up[i][j][k,l] += self.inverse_metric[j,a] * self.inverse_metric[k,b] * self.inverse_metric[l,c] * self.riemann[i][a][k,l]
        self.kretschmann = 0
        for i in range(4):
            for j in range(4):
                for k in range(4):
                    for l in range(4):
                        self.kretschmann += self.riemann_up[i][j][k,l] * self.riemann_down[i][j][k,l]
        self.kretschmann = sp.simplify(self.kretschmann)

    def print_kretschmann(self):
        latex_text = f"k = {sp.latex(self.kretschmann)}"
        display(Math(latex_text))
    
    def check_metric(self):
        metric_T = self.metric.T
        for i in range(4):
            for j in range(4):
                if self.metric[i, j] != metric_T[i, j]:
                    raise Exception("Metric is not symmetric")

    def check_ricci(self):
        ricci_T = self.ricci.T
        for i in range(4):
            for j in range(4):
                if self.ricci[i, j] != ricci_T[i, j]:
                    raise Exception("Ricci tensor is not symmetric")

    def check_einstein(self):
        einstein_T = self.einstein.T
        for i in range(4):
            for j in range(4):
                if self.einstein[i, j] != einstein_T[i, j]:
                    raise Exception("Einstein tensor is not symmetric")
             