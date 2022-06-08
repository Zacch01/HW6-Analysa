import math
import sympy as sp
from sympy.utilities.lambdify import lambdify

ACCURACY = 0.0001

def simpsonMethod(f, startPoint, endPoint, parts):
    """
    :param f: original function
    :param startPoint: start of range
    :param endPoint: end of range
    :param parts: amount of segments
    :return: approximate area of the integral
    """
    if parts % 2 == 1:  # if there is not even numbers of parts
        print("Amount of parts must be even")
        return None
    x = sp.symbols('x')
    func = lambdify(x, f)
    gap = abs(endPoint - startPoint) / parts  # calculate h
    string = "Integral(" + str(startPoint) + ", " + str(endPoint) + ") = 1/3 * " + str(gap) + "[f(" + str(
        startPoint) + ")"
    appr = func(startPoint)  # placing the start point in the function
    for i in range(1, parts):  # run over the parts
        if i % 2 == 0:  # if is the even place
            string += " + 2 * f(" + str((i * gap) + startPoint) + ")"
            appr += 2 * func((i * gap) + startPoint)
        else:  # if is not the even place
            string += " + 4 * f(" + str((i * gap) + startPoint) + ")"
            appr += 4 * func((i * gap) + startPoint)
        if i % 4 == 0:  # for the printing
            string += "\n"
    string += " * f(" + str(endPoint) + ")]\n"
    print(string)  # print the equation
    appr += func(endPoint)
    appr *= 1 / 3 * gap
    return appr


def errorCalculation(f, parts, endPoint, startPoint):
    """
    :param f: original function
    :param parts: amount of segments
    :param endPoint: end of range
    :param startPoint: start of range
    :return: value of error
    """
    x = sp.symbols('x')
    gap = abs(endPoint - startPoint) / parts
    func = calcDerived(calcDerived(calcDerived(calcDerived(f))))
    ksi = findUpperBound(func)
    if ksi is None:
        print("Couldn't calculate the error")
        return
    fourthDerivative = lambdify(x, func)
    return 1 / 180 * math.pow(gap, 4) * (endPoint - startPoint) * fourthDerivative(ksi)


def findUpperBound(f):
    """
    :return: estimated upper bound
    """
    x = sp.symbols('x')
    fTag = calcDerived(f)
    s = sp.solve(fTag)
    func = lambdify(x, f)
    if len(s) > 0:
        maximum = func(s[0])
        for i in range(1, len(s)):
            if func(s[i]) > maximum:
                maximum = func(s[i])
        return maximum
    return None


def calcDerived(f):
    """
    :param f: original func
    :return: the derived without lambdify
    """
    # calc the derivative from func -> lambdify
    x = sp.symbols('x')
    f_prime = f.diff(x)
    return f_prime

def trapezoidMethod(f, a, b, n):
    """
    :param f: Original function
    :param a: start of the range
    :param b: end of the range
    :param n: the number of the segments
    :return: The area in the range
    """
    """if not TrapezoidError(f, a, b, n):
        print("Not good")"""
    x = sp.symbols('x')
    f = lambdify(x, f)
    h = (b - a) / n
    sum = 0
    save = a
    count = 0
    while a < b:
        sum += 0.5 * ((a + h) - a) * (f(a) + f(a + h))
        count += 1
        if a is not save:
            print(" + ", end="")
        if count == 3:
            print("\n       ", end="")
            count = 0
        print("1/2 * (" + str(b) + " - " + str(a) + ") * (f(" + str(a) + " + f(" + str(a + h) + "))", end="")
        a += h
    return sum


def MaxFunctionValue(f, startAt, endAt):
    """
    Method for finding the function Roots

    :param f: Our function
    :param startAt: Left domain of the function
    :param endAt: Right domain of the function
    """
    # Variable to store the derivative function
    g = f.diff(x)

    # Activating the functions to be able to get an X
    f = lambdify(x, f)
    g = lambdify(x, g)

    # Variable to store the maximum iteration in order to find the function roots
    maxIteration = 1000

    # Variable to store the max value of the function
    maxValue = max(f(startAt), f(endAt))

    # Divide our function domain range into multiply domains with 0.1 range, then search for each one of them for a root
    while startAt < endAt:

        # In case the derivative function changes its sign (Mean there's a possibility for an extreme point)
        if g(startAt) * g(startAt + 0.1) < 0:

            # Getting a possibility for an extreme point (Might be a Root or an Extreme point)
            possiblePoint = Secant(g, startAt, startAt + 0.1, maxIteration)

            # In case we found an extreme point
            if f(possiblePoint) > maxValue:
                maxValue = f(possiblePoint)

        # Update our domain for the next iteration
        startAt = startAt + 0.1

    return maxValue


def Secant(f, previewX, currentX, maxIteration):
    """
    Finding the function root in the domain range [left To right]

    :param f: Our function
    :param previewX: Left domain of the function
    :param currentX: Right domain of the function
    :param maxIteration: The maximum iteration for finding the root
    :return: The root of the function if existed, else according failed message
    """
    # Search the root within the maximum allowed iteration
    for i in range(maxIteration):

        # Variable to store the next X
        nextX = (previewX * f(currentX) - currentX * f(previewX)) / (f(currentX) - f(previewX))

        # In case we found our root, Return the root and the iteration number
        if abs(f(nextX)) < ACCURACY:
            return int(nextX * 10 * 5) / 10 * 5

        # Update the previewX to be the currentX
        previewX = currentX

        # Update the currentX to be new one
        currentX = nextX

    # In case we didn't find the root within the allowed amount iteration, Print fail message and shut down the program
    print("Failed to find the root, Secant Method isn't suitable")




x = sp.symbols('x')

# define function
f = (sp.sin(x ** 2 + 5 * x + 6)) / (2 * sp.exp(-x))


# define range
startPoint = 0
endPoint = 1
print(simpsonMethod(f, startPoint, endPoint, 4))
#print(errorCalculation(f, 4, endPoint, startPoint))
print()

#result = trapezoidMethod(f, startPoint, endPoint, 4)
#print("\n" + str(result)+ "\n"+ "\n"+ "\n")




