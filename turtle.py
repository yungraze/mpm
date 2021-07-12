import turtle, time
import random

window = turtle.Screen()
window.setup(700, 700)

board = turtle.Turtle()
board.speed(0)
board.up()
board.hideturtle()
board.pensize(5)
board.goto(300, 300)
board.down()
board.goto(300, -300)
board.goto(-300, -300)
board.goto(-300, 300)
board.goto(300, 300)

balls = []
count = 5

for i in range(count):
    ball = turtle.Turtle()
    ball.shape('circle')

    randx = random.randint(-250, 250)
    randy = random.randint(-250, 250)
    ball.speed(0)
    ball.up()
    ball.goto(randx, randy)
    ball.down()

    dx = random.randint(-5, 5)
    dy = random.randint(-5, 5)

    balls.append([ball, dx, dy])


for i in range(0, 1000, 1):
    for j in range(count):
        x, y = balls[j][0].position()
        if x+balls[j][1] >= 300 or x+balls[j][1] <= -300:
            balls[j][1] = - balls[j][1]
        if y+balls[j][2] >= 300 or y+balls[j][2] <= -300:
            balls[j][2] = - balls[j][2]
        balls[j][0].goto(x+balls[j][1], y+balls[j][2])
