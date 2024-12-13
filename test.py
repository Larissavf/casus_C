import math
from scipy.optimize import curve_fit
from scipy.integrate import solve_ivp


"""
A class with different models to possibly fit a tumor growth.
Different methods to either solve the differentiate formule, fit data to the
model to determine the best fitting variables and to test wich model fits the best.
"""
class TGmodels:

    """
    Contains the formule for the model 
    """
    def formule():
        return None
    """
    Containts the params in a dict format
    """
    def params():
        return None
    
    """
    Set the params to new values
    """
    def setparams():
        return None

    """
    To solve the differantiate formule using steps to determine 
    using the Runge Kutta method
    input: int with the point
    ouput: int with the outcome  
    """
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
    
    """
    To solve the differantiate formule using steps to determine 
    using the Euler's method
    input: int with the point
    ouput: int with the outcome 
    """
    def solve_euler(self, t):
        self.solvingsystem = "Euler"
        steps = max(1, int(abs(t) / 0.01))
        dt = t / steps
        y = self.y()
        for step in range(steps):
            dydt = self.formule(y)                 # Differentiaalvergelijking
            y += dydt * dt
        return y
    
    """ 
    Optimalisation of model to known points.
    input: x data, y data, string with type of solvingsystem
    output: float with the MSE
    """
    def mean_squared_error(self, x, y, solvingsystem = None):
        solvingsystems = {
            "RungeKutta": self.solve_rungekutta,
            "Euler": self.solve_euler
        }
        # if chosen for a different solvingsystem
        if solvingsystem not in solvingsystems:
            if solvingsystem == None:
                solvingsystem = "RungeKutta"
            else: return f"Ongeldig solvesystem '{solvingsystem}'. Kies uit: {', '.join(solvingsystems.values())}"

        # determine MSE
        N = len(x)
        total = 0.0
        for i in range(N):
            error = y[i] - solvingsystems[solvingsystem](x[i])
            total += error * error
        return total / N
    

    """
    Fit the model with the correct params to the data using 
    the Hooke & Jeeves, direct search.
    input: x data and y data
    """
    def fit(self, x, y):
        params = self.params()

        deltas = {key: 1.0 for key in params} 
        # Herhaaldelijke aanpassing
        mse = self.mean_squared_error(x, y)
        while max(abs(delta) for delta in deltas.values()) > 1e-8: 
            for key in params:
                new_params = params.copy()
                # Probeer de betreffende parameter the verhogen
                new_params[key] = params[key] + deltas[key]
                if new_params[key] < 0:
                    new_params[key] = 1e-8
                self.setparams(**new_params)
                new_mse = self.mean_squared_error(x, y)
                if new_mse < mse:
                    params = new_params
                    mse = new_mse
                    deltas[key] *= 1.2 
                    continue
                # Probeer de betreffende parameter the verlagen
                new_params[key] = params[key] - deltas[key]
                if new_params[key] < 0:
                    new_params[key] = 1e-8
                self.setparams(**new_params)
                new_mse = self.mean_squared_error(x, y)
                if new_mse < mse:
                    params = new_params
                    mse = new_mse
                    deltas[key] *= -1.0 
                    continue
                # Verklein de stapgrootte
                deltas[key] *= 0.5 
        self.setparams(**params)

    """
    Determine the BIC value to compare the fit of models.
    """
    def BIC(self, x, y ):
        return len(x) * math.log(self.mean_squared_error(x, y)) + (math.log(len(x)) * len(self.params()))
        
"""
The lineair growth model
"""
class Lineair_growth(TGmodels):
    def __init__(self,y , c =  1.0):
        self.c = c
        self.y = y

    def formule(self, y):
        return self.c
    
    def params(self):
        return {"y": self.y, "c": self.c}
    
    def setparams(self, y, c):
        self.c = c
        self.y = y

"""
The exponantial increasing growth model
"""
class Exponentieel_toenemende_groei(TGmodels):
    def __init__(self,y , c =  1.0):
        self.c = c
        self.y = y

    def formule(self, y):
        return self.c *y
    
    def params(self):
        return {"y": self.y, "c": self.c}
    
    def setparams(self, y, c):
        self.c = c
        self.y = y

"""
The exponantiol decrease growth model
"""
class Exponentieel_afvlakkende_groei(TGmodels):
    def __init__(self,y , vmax, c =  1.0):
        self.c = c
        self.y = y
        self.vmax = vmax

    def formule(self, y):
        if y >= self.vmax:
            y = self.vmax
        return self.c *(self.vmax -y)
    
    def params(self):
        return {"y": self.y, "c": self.c, "vmax": self.vmax}
    
    def setparams(self, y, c, vmax):
        self.c = c
        self.y = y
        self.vmax = vmax

"""
The mendelsohn growth model
"""             
class Mendelsohn_growth(TGmodels):
    def __init__(self, y ,c = 1.0, d = 1.0):
        self.c = c
        self.d = d
        self.y = y

    def formule(self, y):
        return self.c*y**self.d

    def params(self):
        return {"c": self.c, "d": self.d, "y": self.y}
    
    def setparams(self, c, d, y):
        self.c = c
        self.d = d
        self.y = y


"""
The logistic growth model
"""    
class Logistic_growth(TGmodels):
    def __init__(self, y, vmin, c = 1.0):
        self.c = c
        self.vmin = vmin
        self.y = y
    
    def formule(self, y):
        if y <= self.vmin:
            y= self.vmin
        return self.c*y *(self.vmin - y)
    
    def params(self):
        return {"c": self.c, "y": self.y, "vmin": self.vmin}
    
    def setparams(self, c, y, vmin):
        self.c = c
        self.vmin = vmin
        self.y = y

"""
The allee effect growth model
"""
class Allee_effect_growth(TGmodels):
    def __init__(self,y , vmin, vmax, c = 1.0 ):
        self.c = c
        self.vmin = vmin
        self.vmax = vmax
        self.y = y
    
    def formule(self, y):
        if y >= self.vmax:
            y = self.vmax
        if y <= self.vmin:
            y= self.vmin
        return self.c*(y - self.vmin) * (self.vmax - y)

    def params(self):
        return {"c": self.c, "y": self.y, "vmin": self.vmin, "vmax": self.vmax}
    
    def setparams(self, c, y, vmin, vmax):
        self.c = c
        self.vmin = vmin
        self.vmax = vmax
        self.y = y

"""
The lineair limited growth model
"""
class Lineair_gelimiteerde_groei(TGmodels):
    def __init__(self,y , c = 1.0, d = 1.0):
        self.c = c
        self.d = d
        self.y = y
    
    def formule(self, y ):
        return self.c* (y/(y+self.d))
    
    def params(self):
        return {"c": self.c, "d": self.d, "y": self.y}
    
    def setparams(self, c, d, y):
        self.c = c
        self.d = d
        self.y = y

"""
The surface limited growth model
"""    
class surface_limited_growth(TGmodels):
    def __init__(self,y , c = 1.0, d = 1.0):
        self.c = c
        self.d = d
        self.y = y
    
    def formule(self, y ):
        return self.c* (y/(y+self.d)**(1/3))
    
    def params(self):
        return {"c": self.c, "d": self.d, "y": self.y}
    
    def setparams(self, c, d, y):
        self.c = c
        self.d = d
        self.y = y

"""
The von bertalanffy limited growth model
"""  
class Von_Bertalanffy_groei(TGmodels):
    def __init__(self,y , c = 1.0, d = 1.0):
        self.c = c
        self.d = d
        self.y = y
    
    def formule(self, y ):
        return self.c* y** (2/3) - (y*self.d)
    
    def params(self):
        return {"c": self.c, "d": self.d, "y": self.y}
    
    def setparams(self, c, d, y):
        self.c = c
        self.d = d
        self.y = y

"""
The gompertz growth model
"""  
class Gompertz_growth(TGmodels):
    def __init__(self, y,  vmax, c = 1.0):
        self.c = c
        self.vmax = vmax
        self.y = y
    
    def formule(self, y):
        if y >= self.vmax:
            y = self.vmax
        return self.c*y *math.log(self.vmax / y) 
    
    def params(self):
        return {"y": self.y, "c": self.c, "vmax": self.vmax}
    
    def setparams(self, y, c, vmax):
        self.c = c
        self.vmax = vmax
        self.y = y

