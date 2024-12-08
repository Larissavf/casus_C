import math

class TGmodels:
    def formule(y):
        return True

    def solve_rungekutta(self, y, t, steps=10):
        """Lost de differentiaalvergelijking
            dy/dt = y + 1
        met beginconditie y(0) = 0 op met Runge-Kutta's methode,
        en retourneert de waarde y(t) voor een gevraagde t.
        """
        dt = t / steps
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


class Lineair_growth(TGmodels):
    def __init__(self, c):
        self.c = c

    def formule(self, y):
        return self.c
    
class Mendelsohn_growth(TGmodels):
    def __init__(self, c, d):
        self.c = c
        self.d = d
    
    def formule(self, y):
        return self.c*y**self.d
    
class Logistic_growth(TGmodels):
    def __init__(self, c, vmin):
        self.c = c
        self.vmin = vmin
    
    def formule(self, y):
        return self.c*y *(self.vmin - y)

class Allee_effect_growth(TGmodels):
    def __init__(self, c, vmin, vmax):
        self.c = c
        self.vmin = vmin
        self.vmax = vmax
    
    def formule(self, y):
        return self.c*(y - self.vmin) * (self.vmax - y)
    
class surface_limited_growth(TGmodels):
    def __init__(self, c, d):
        self.c = c
        self.d = d
    
    def formule(self, y):
        return self.c* (y/(y+self.d)**(1/3))

class Gompertz_growth(TGmodels):
    def __init__(self, c, vmax):
        self.c = c
        self.vmax = vmax
    
    def formule(self, y):
        return self.c*y *math.log(self.vmax / y) 

