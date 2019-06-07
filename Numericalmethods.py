import math
class Numericalmethods:
    #Function that returns the numerical derivative of f at point x using an interval h
    def df(self, f, x, h):
        return (f(x+h)-f(x-h))/2/h

    #Returns the zero of function f using x0 as the starting point with a maximum error e and using intervals h for the diferentiation
    def newton(self, f, x0, e, h):
        x=x0
        error=abs(f(x))
        while(error>e):
            x=x-f(x)/Numericalmethods().df(f, x, h)
            error=abs(f(x))
        return x
