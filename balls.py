import turtle, time
import random
import numpy as np

window = turtle.Screen()
window.setup(400, 400)

qwe = turtle.Turtle()
qwe.speed(0)
qwe.up()
qwe.hideturtle()
qwe.pensize(5)
qwe.goto(200, 200)
qwe.down()
qwe.goto(200, -200)
qwe.goto(-200, -200)
qwe.goto(-200, 200) 
qwe.goto(200, 200)

planets = []
# planets.append([planet, x, v, R, m])

earth = turtle.Turtle()
earth.turtlesize(0.8)
earth.speed(1)
earth.up()
earth.goto(-100, 0)
earth.down()
planets.append([earth, -100, 1, 50, 10])

merc = turtle.Turtle()
merc.shape('circle')
merc.turtlesize(0.5)
merc.color('brown')
merc.speed(1)
merc.up()
merc.goto(100, 0)
merc.down()
planets.append([merc, 100, -1, 50, 5])
koof = 1

coords = [0, 0]
count = 5

    
for t in range(0, 3000, 1):
  for j in range(len(planets)):
    x, y = planets[j][0].position()
    if x + planets[j][1] >= 200 or x +planets[j][1] <= -200:
      planets[j][1] = - planets[j][1]
    planets[j][0].goto(x+ planets[j][1], 0)
    planets[j][0].goto(planets[j][1]+planets[j][2]*t, 0)
    for k in range(len(planets)):
      coords[k], y = planets[j][0].position()
    if np.absolute(coords[0]-coords[1]) <= planets[0][3] + planets[1][3]:
      planets[j][2] = planets[j][2] + (1 + koof * planets[j][4] * (planets[j][2]-planets[j][2])/(planets[j][4]+planets[j][4]))
    
