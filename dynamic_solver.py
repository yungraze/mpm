import fenics as fen


def dynamic_solver_func(ncells=10,  # количество узлов на заданном итервале
                        init_time=0,  # начальный момент времени
                        end_time=10,  # конечный момент времени
                        dxdphi=1,  # производная от потенциала по х
                        dydphi=1,  # производня от потенциала по у
                        x0=0,  # начальное положение по оси х
                        vx0=1,  # проекция начальной скорости на ось х
                        y0=0,  # начальное положение по оси у
                        vy0=1):  # проекция начальной скорости на ось у
    """ Функция на вход которой подается производная от потенциала
        гравитационного поля, возвращающая координаты смещенных
        материальных точек (частиц).
    """
    # генерация сетки на заданном интервале времени
    mesh = fen.IntervalMesh(ncells, 0, end_time-init_time)

    welm = fen.MixedElement([fen.FiniteElement('Lagrange', fen.interval, 2),
                             fen.FiniteElement('Lagrange', fen.interval, 2),
                             fen.FiniteElement('Lagrange', fen.interval, 2),
                             fen.FiniteElement('Lagrange', fen.interval, 2)])

    # генерация функционального рростаанства
    W = fen.FunctionSpace(mesh, welm)

    # постановка начальных условий задачи
    bcsys = [fen.DirichletBC(W.sub(0), fen.Constant(x0), 'near(x[0], 0)'),
             fen.DirichletBC(W.sub(1), fen.Constant(vx0), 'near(x[0], 0)'),
             fen.DirichletBC(W.sub(2), fen.Constant(y0), 'near(x[0], 0)'),
             fen.DirichletBC(W.sub(3), fen.Constant(vy0), 'near(x[0], 0)')]

    # опееделение тестовых функций для решения задачи
    up = fen.Function(W)
    x_cor, v_x, y_cor, v_y = fen.split(up)
    v1, v2, v3, v4 = fen.split(fen.TestFunction(W))

    # постановка задачи drdt = v; dvdt = - grad(phi) в проекциях на оси системы координат
    weak_form = (x_cor.dx(0) - v_x) * v1 * fen.dx + (v_x.dx(0) + dxdphi) * v2 * fen.dx \
                + (y_cor.dx(0) - v_y) * v3 * fen.dx + (v_y.dx(0) + dydphi) * v4 * fen.dx

    # решние поставленной задачи
    fen.solve(weak_form == 0, up, bcs=bcsys)

    # определение момента времени
    time = fen.Point(end_time - init_time)

    # расчет координат и скоростей
    x_end_time = up(time.x())[0]
    vx_end_time = up(time.x())[1]
    y_end_time = up(time.x())[2]
    vy_end_time = up(time.x())[3]

    return x_end_time, y_end_time, vx_end_time, vy_end_time


print(dynamic_solver_func())
