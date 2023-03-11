def getHeaderTable(x):
    answer = "    t      delta           x:"
    for xi in x:
        answer += "%9.3f" % xi
    answer += "\n"
    return answer

def table_for_shooting_method(history):
    answer = " Itr  z(0)      y(1)      Delta\n"
    for i in history:
        answer += "{:2d}  {:7.5f}  {:8.6f}  {:18.14e}\n".format(i[0],i[1],i[2],i[3])
    return answer

def table_for_runge_kutta(history):
    answer = "   x      y(x)     Ypr(x)     z(x)       Dela\n"
    for i in history:
        answer += "{:7.5f}  {:7.5f}  {:7.5f}  {:7.5f} {:18.14e}\n".format(i[0],i[1],i[2],i[3],i[4])
    return answer