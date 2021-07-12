import turtle, time
import random
import numpy as np

window = turtle.Screen()
window.setup(700, 700)
window.bgcolor('black')

sun = turtle.Turtle()
sun.shape('circle')
sun.turtlesize(3)
sun.color('yellow')

planets = []
# planets.append([planet, a, b, c, T])

earth = turtle.Turtle()
a = 100
b = 100
T = 500
c = np.sqrt(a**2 - b**2)
earth.shape('circle')
earth.color('green')
earth.turtlesize(0.8)
earth.speed(0)
earth.up()
earth.goto(a-c, 0)
earth.down()
planets.append([earth, a, b, c, T])

merc = turtle.Turtle()
a = 39
b = 39
T = 100
c = np.sqrt(a**2 - b**2)
merc.shape('circle')
merc.turtlesize(0.5)
merc.color('brown')
merc.speed(0)
merc.up()
merc.goto(a-c, 0)
merc.down()
planets.append([merc, a, b, c, T])

comet = turtle.Turtle()
a = 259
b = 150
T = 700
c = np.sqrt(a**2 - b**2)
comet.shape('circle')
comet.turtlesize(0.5)
comet.color('brown')
comet.speed(0)
comet.up()
comet.goto(a-c, 0)
comet.down()
planets.append([comet, a, b, c, T])


for t in range(0, 3000, 1):
    for j in range(len(planets)):
        planets[j][0].goto(planets[j][1]*np.cos(2*np.pi/planets[j][4]*t) - planets[j][3],
                           planets[j][2]*np.sin(2*np.pi/planets[j][4]*t))
