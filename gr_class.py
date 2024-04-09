import sympy as sp

sp.init_printing()

class Manifold: 
    def __init__(self, coordinates, metric):
        #Implement coordinates and metric
        self.coordinates = coordinates
        self.metric = metric
        self.inverse_metric = metric.inv()

        self.compute_christoffel()
        self.compute_riemann()
        self.compute_kretschmann()
        self.compute_ricci()
        self.compute_scalar()
        self.compute_einstein()

        self.check_metric()
        self.check_ricci()
        self.check_einstein()

    # Compute geometrical quantities
    def compute_christoffel(self):
        self.christoffel = [sp.zeros(4, 4),  sp.zeros(4, 4), sp.zeros(4, 4), sp.zeros(4, 4)]
        for i in range(4):
            for j in range(4):
                for k in range(4):
                    for l in range(4):
                        self.christoffel[k][i,j] += 1/2 * self.inverse_metric[k,l] * (sp.diff(self.metric[i,l], self.coordinates[j]) + sp.diff(self.metric[j,l], self.coordinates[i])-sp.diff(self.metric[i,j], self.coordinates[l]))
                self.christoffel[i][j,k] = sp.simplify(self.christoffel[i][j,k])

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
            self.scalar= sp.simplify(self.scalar)

    def compute_einstein(self): 
        self.einstein = sp.zeros(4, 4)
        for i in range(4):
            for j in range(4):
                self.einstein[i,j] += self.ricci[i,j] - 0.5 * self.metric[i,j] * self.scalar
                self.einstein[i,j] = sp.simplify(self.einstein[i,j])

    # Check symmetries
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
             
    # Get geometrical quantities
    def get_metric(self):
        return self.metric

    def get_inverse_metric(self):
        return self.inverse_metric
    
    def get_christoffel(self):
        return self.christoffel

    def get_riemann(self):
        return self.riemann

    def get_ricci(self):
        return self.ricci

    def get_scalar(self):
        return self.scalar

    def get_einstein(self):
        return self.einstein

    def get_kretschmann(self):
        return self.kretschmann
    
    # Print geometrical quantities in LaTeX
    def print_metric(self):
        for i in range(4):
            for j in range(4):
                if self.metric[i, j] != 0:
                    print(r"\begin{equation*}")
                    print(r"g_{", i, j,"} = ", sp.latex(self.metric[i, j]))
                    print(r"\end{equation*}")

    def print_inverse_metric(self):
        for i in range(4):
            for j in range(4):
                if self.inverse_metric[i, j] != 0:
                    print(r"\begin{equation*}")
                    print(r"g^{", i, j,"} = ", sp.latex(self.inverse_metric[i, j]))
                    print(r"\end{equation*}")

    def print_christoffel(self):
        for i in range(4):
            for j in range(4):
                for k in range(4):
                    if self.christoffel[i][j, k] != 0:
                        print(r"\begin{equation*}")
                        print(r"\Gamma^{", i, "}_{", j, k, "} = ", sp.latex(self.christoffel[i][j, k]))
                        print(r"\end{equation*}")

    def print_riemann(self):
        for i in range(4):
            for j in range(4):
                for k in range(4):
                    for l in range(4):
                        if self.riemann[i][j][k, l] != 0:
                            print(r"\begin{equation*}")
                            print(r"R^{", i, r"}_{ \phantom", i, j, k, l, "} = ", sp.latex(self.riemann[i][j][k, l]))
                            print(r"\end{equation*}")

        
    def print_ricci(self):
        for i in range(4):
            for j in range(4):
                if self.ricci[i, j] != 0:
                        print(r"\begin{equation*}")
                        print(r"R_{", i, j, "} = ", sp.latex(self.ricci[i,j]))
                        print(r"\end{equation*}")

    def print_scalar(self):
        print(r"\begin{equation*}")
        print(r"R = ", sp.latex(self.scalar))
        print(r"\end{equation*}")

    def print_einstein(self):
        for i in range(4):
            for j in range(4):
                if self.einstein[i, j] != 0:
                    print(r"\begin{equation*}")
                    print(r"G_{", i, j, "} = ", sp.latex(self.einstein[i, j]))
                    print(r"\end{equation*}")

    def print_kretschmann(self):
        print(r"\begin{equation*}")
        print(r"k = ", sp.latex(self.kretschmann))
        print(r"\end{equation*}")
        
        
