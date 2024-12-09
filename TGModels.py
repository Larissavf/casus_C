import math

class TGmodels:
    def __init__(self):
        self.solvingsystem = False

    def formule():
        return None
    def params():
        return None
    def setparams():
        return None

    def solve_rungekutta(self, t):
        self.solvingsystem = "RungeKutta"

        steps = max(1, int(abs(t) / 0.01))
        dt = t / steps
        y = self.y
        for step in range(steps):
            # Tijdelijke stappen:
            dydt1 = self.formule(y)                # Differentiaalvergelijking stap 1
            y1 = y + dydt1 * 0.5 * dt
            dydt2 = self.formule(y1)                # Differentiaalvergelijking stap 2
            y2 = y + dydt2 * 0.5 * dt
            dydt3 = self.formule(y2)                # Differentiaalvergelijking stap 3
            y3 = y + dydt3 * dt
            dydt4 = self.formule(y3)                # Differentiaalvergelijking stap 4
            # Definitieve stap:
            dydt = (dydt1 + 2.0 * dydt2 + 2.0 * dydt3 + dydt4) / 6.0
            y += dydt * dt
        return y
    
    def mean_squared_error(self, x, y, solvingsystem = None):
        solvingsystems = {
            "RungeKutta": self.solve_rungekutta
        }


        # issue
        solvingsystem = "RungeKutta"

            # if self.solvingsystem != False:
            #     solvingsystem = self.solvingsystem
            # if solvingsystem not in solvingsystems:
            #     return f"Ongeldig solvesystem '{solvingsystem}'. Kies uit: {', '.join(solvingsystems.values())}"

        N = len(x)
        total = 0.0
        for i in range(N):
            error = y[i] - solvingsystems[solvingsystem](x[i])
            total += error * error
        return total / N
    

    def fit(self, x, y):
        params = self.params()

        deltas = {key: 1.0 for key in params} #adjusting?
        counter = 1
        # Herhaaldelijke aanpassing
        mse = self.mean_squared_error(x, y)
        while max(abs(delta) for delta in deltas.values()) > 0.0001: #anders?
            for key in params:
                new_params = params.copy()
                # Probeer de betreffende parameter the verhogen
                new_params[key] = params[key] + deltas[key]
                self.setparams(**new_params)
                new_mse = self.mean_squared_error(x, y)
                if new_mse < mse:
                    params = new_params
                    mse = new_mse
                    deltas[key] *= 1.2 #adjust?
                    continue
                # Probeer de betreffende parameter the verlagen
                new_params[key] = params[key] - deltas[key]
                self.setparams(**new_params)
                new_mse = self.mean_squared_error(x, y)
                if new_mse < mse:
                    params = new_params
                    mse = new_mse
                    deltas[key] *= -1.0 #adjust?
                    continue
                # Verklein de stapgrootte
                deltas[key] *= 0.5 #adjust?
            counter += 1
            print(f"loop{counter}")
        self.setparams(**params)
        


class Lineair_growth(TGmodels):
    def __init__(self,y , c =  0.0):
        self.c = c
        self.y = y

    def formule(self, y):
        return self.c
    
    def params(self):
        return {"y": self.y, "c": self.c}
    
    def setparams(self, y, c):
        self.c = c
        self.y = y
        
class Mendelsohn_growth(TGmodels):
    def __init__(self, y ,c = 0.0, d = 0.0):
        self.c = c
        self.d = d
        self.y = y

    def formule(self, y):
        if y <= 0:
            self.y = 1e-8
            y = 1e-8
        return self.c*y**self.d

    def params(self):
        return {"c": self.c, "d": self.d, "y": self.y}
    
    def setparams(self, c, d, y):
        self.c = c
        self.d = d
        self.y = y
    
class Logistic_growth(TGmodels):
    def __init__(self, y, vmin, c = 0.0):
        self.c = c
        self.vmin = vmin
        self.y = y
    
    def formule(self, y):
        if y <= 0:
            y = 1e-8
            self.y = 1e-8
        return self.c*y *(self.vmin - y)
    
    def params(self):
        return {"c": self.c, "y": self.y, "vmin": self.vmin}
    
    def setparams(self, c, y, vmin):
        self.c = c
        self.vmin = vmin
        self.y = y


class Allee_effect_growth(TGmodels):
    def __init__(self,y , c, vmin, vmax ):
        self.c = c
        self.vmin = vmin
        self.vmax = vmax
        self.y = y
    
    def formule(self, y):
        if y <= 0:
            y = 1e-8
            self.y = 1e-8
        return self.c*(y - self.vmin) * (self.vmax - y)

    def params(self):
        return {"c": self.c, "y": self.y, "vmin": self.vmin, "vmax": self.vmax}
    
    def setparams(self, c, y, vmin, vmax):
        self.c = c
        self.vmin = vmin
        self.vmax = vmax
        self.y = y

    
class surface_limited_growth(TGmodels):
    def __init__(self,y , c = 0.0, d = 0.0):
        self.c = c
        self.d = d
        self.y = y
    
    def formule(self, y):
        if y <= 0:
            y = 1e-8
            self.y = 1e-8
        return self.c* (y/(y+self.d)**(1/3))
    
    def params(self):
        return {"c": self.c, "d": self.d, "y": self.y}
    
    def setparams(self, c, d, y):
        self.c = c
        self.d = d
        self.y = y


class Gompertz_growth(TGmodels):
    def __init__(self, y,  vmax, c = 0.0):
        self.c = c
        self.vmax = vmax
        self.y = y
    
    def formule(self, y):
        if y <= 0:
            y = 1e-8
            self.y = 1e-8
        return self.c*y *math.log(self.vmax / y) 
    
    def params(self):
        return {"y": self.y, "c": self.c, "vmax": self.vmax}
    
    def setparams(self, y, c, vmax):
        self.c = c
        self.vmax = vmax
        self.y = y

