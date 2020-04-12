from math import *
from graphics import *
from graphics import Point as Poing
from copy import deepcopy as dc
import cProfile
inf = float('inf')
PROFILE = False
PRINT = False
CRASH = False
DEBUG = False
if not PRINT:
    CRASH = False
    DEBUG = False
    _print = print
    print('Printing is disabled.')
    def print(*args):
        pass

_sqrt = sqrt
def sqrt(x):
    if equal(x, 0): return 0
    return _sqrt(x)

global n
n = None
global g
g = -9.8
global dirs
dirs = ('up', 'down', 'none', 'left', 'right')
global tau
tau = 2*pi

STOP = False
INVALID = False

def between(num, lower, upper):
    return lower < num and upper > num or lower > num and upper < num

def mod(a, q):
    return ((a % q) + q) % q;

def rightq(a):
    if a > pi/2 or a < -pi/2:
        return pi - a
    else:
        return a

def otherq(a):
    return pi - a

def norm(a):
    return mod(a + pi, 2*pi) - pi

def norm360(a):
    return mod(a + 180, 360) - 180

def equal(var1, var2):
    if var1 in dirs and var2 in dirs:
        return var1 == var2
    try:
        var1 = float(var1)
        var2 = float(var2)
        return abs(var1 - var2) < 0.0000000001 or abs((var1 - var2) / var1) < 0.00001
    except (TypeError, ValueError):
        return False
    except ZeroDivisionError:
        return abs(var1 - var2) < 0.00001

def zero(var):
    return equal(var, 0)

def fround(num, places=0, func=round):
    return func(num * 10 ** places) / float(10 ** places)

def sigfigs(num, figs):
    if equal(num, 0):
        return 0
    if num >= 0: factor = num
    else: factor = -num
    power = figs - int(ceil(log(factor, 10)))
    magnitude = 10 ** power
    shifted = round(num * magnitude)
    return shifted / magnitude

def isfloat(n):
    try:
        float(n)
        return True
    except Exception:
        return False

def enableField(field):
    i = fields.index(field)
    checks[i].checkon(win)
    fields[i].entry.configure(state='normal')

def checkUpdate(var, i):
    global STOP
    if STOP: return
    invalidValues.undraw()
    nums = []
    for field in fields:
        if not checks[fields.index(field)].checked:
            STOP = True
            field.setText('')
            STOP = False
        try:
            t = field.getText()
            if t in dirs:
                nums.append(t)
            else:
                nums.append(float(t))
        except (TypeError, ValueError):
            if not between(fields.index(field), len(fields)-5, len(fields)) and t != '-' and not (len(t) >= 1 and t[-1] == 'e' and isfloat(t[:-1])) and not (len(t) >= 2 and t[-2:] == 'e-' and isfloat(t[:-2])) and not (len(t) >= 2 and t[-2:] == 'e+' and isfloat(t[:-2])):
                STOP = True
                fields[fields.index(field)].setText('')
                STOP = False
            nums.append(n)
    if nums[6] != n: nums[6] *= pi / 180; nums[6] = norm(nums[6])
    if nums[13] != n: nums[13] *= pi / 180; nums[13] = norm(nums[13])
    if nums[20] != n: nums[20] *= pi / 180; nums[20] = norm(nums[20])
    (ps.i.x, ps.i.y, ps.i.v, ps.i.vx, ps.i.vy, ps.i.t, ps.i.a,
     ps.f.x, ps.f.y, ps.f.v, ps.f.vx, ps.f.vy, ps.f.t, ps.f.a,
     ps.c.x, ps.c.y, ps.c.v, ps.c.vx, ps.c.vy, ps.c.t, ps.c.a,
     ps.i.d, ps.f.d, ps.c.d) = nums
    ps.solve()
    nums = [ps.i.x, ps.i.y, ps.i.v, ps.i.vx, ps.i.vy, ps.i.t, ps.i.a,
            ps.f.x, ps.f.y, ps.f.v, ps.f.vx, ps.f.vy, ps.f.t, ps.f.a,
            ps.c.x, ps.c.y, ps.c.v, ps.c.vx, ps.c.vy, ps.c.t, ps.c.a,
            ps.i.d, ps.f.d, ps.c.d]
    if nums[6] != n: nums[6] *= 180 / pi; nums[6] = norm360(nums[6])
    if nums[13] != n: nums[13] *= 180 / pi; nums[13] = norm360(nums[13])
    if nums[20] != n: nums[20] *= 180 / pi; nums[20] = norm360(nums[20])
    global INVALID
    if not INVALID:
        for i in range(len(nums)):
            try:
                # Only bother checking to disable the field if
                # it's enabled.
                if checks[i].checked:
                    # Only overwrite the field if it has a valid value
                    text = fields[i].getText()
                    if i in (6, 13, 20):
                        try:
                            text = str(norm360(float(nums[i])))
                        except TypeError:
                            pass
                    if text or nums[i] == n: # ... or if it's blank
                        try:
                            float(text) # If it's a number, then proceed.
                        except ValueError:
                            assert text in dirs
                            # If it's not a number, then only proceed if it's
                            # a valid word.
                    # Only overwrite the field if the value is wrong
                    if not equal(text, nums[i]):
                        checks[i].checkoff(win)
                        fields[i].entry.configure(state='readonly')
            except AssertionError:
                pass
            finally:
                if not checks[i].checked:
                    # If the value is consistent, only change if not entered.
                    # Otherwise, uncheck the field, then set the value.
                    if nums[i] == n:
                        STOP = True
                        fields[i].setText('')
                        STOP = False
                        checks[i].checkoff(win)
                    elif nums[i] in dirs:
                        STOP = True
                        fields[i].setText(nums[i])
                        STOP = False
                        checks[i].checkcheck(win)
                    else:
                        STOP = True
                        fields[i].setText(str(float(sigfigs(nums[i], 5))))
                        STOP = False
                        checks[i].checkcheck(win)
    for line in ps.lines:
        line.undraw()
    if INVALID:
        for field in fields:
            if not checks[fields.index(field)].checked:
                STOP = True
                field.setText('')
                STOP = False
        try:
            invalidValues.draw(win)
        except GraphicsError:
            pass
        INVALID = False
    else:
        ps.plot(win)
    win.flush()

def ex(self, form):
    args = form.split(' ')
    tests = []
    for arg in args:
        tests.append(eval('self.' + arg[0] + '.' + arg[1:] + ' != None'))
    for test in tests:
        if not test:
            return False
    return True

def transform(oldcoords1, oldcoords2, newcoords1, newcoords2, point):
    ratio = ((point[0] - oldcoords1[0]) / (oldcoords2[0] - oldcoords1[0]),
             (point[1] - oldcoords1[1]) / (oldcoords2[1] - oldcoords1[1]))
    relative = (ratio[0] * (newcoords2[0] - newcoords1[0]),
                ratio[1] * (newcoords2[1] - newcoords1[1]))
    final = (newcoords1[0] + relative[0],
             newcoords1[1] + relative[1])
    return final

class Pointset():
    def __init__(self, i=n, f=n, c=n, bother=True):
        self.i = i
        self.f = f
        self.c = c
        self.bother = bother
        self.lines = []
    def deduce(self):
        # No negative velocities!
        if ex(self, 'iv') and self.i.v < 0 and not zero(self.i.v) or ex(self, 'fv') and self.f.v < 0 and not zero(self.f.v):
            self.invalid()
        self.set('cvx', 0)
        # Change calculations
        for var in ('x', 'y', 'v', 'vx', 'vy', 't', 'a'):
            if ex(self, 'i'+var+' f'+var):
                self.set('c'+var, eval('self.f.'+var+' - self.i.'+var))
            if ex(self, 'i'+var+' c'+var):
                self.set('f'+var, eval('self.i.'+var+' + self.c.'+var))
            if ex(self, 'f'+var+' c'+var):
                self.set('i'+var, eval('self.f.'+var+' - self.c.'+var))
        # Same point detection
        if ex(self, 'cx') or ex(self, 'cvy') or ex(self, 'ct') or ex(self, 'ca'):
            if (zero(self.c.vy) or
                zero(self.c.t) or
                zero(self.c.x) and self.i.vx != n and not zero(self.i.vx)):
                # Last: zero ∂x AND not simply vertical
                self.set('cx', 0)
                self.set('cy', 0)
                self.set('cv', 0)
                self.set('cvx', 0)
                self.set('cvy', 0)
                self.set('ct', 0)
                self.set('ca', 0)
        if ex(self, 'cy') or ex(self, 'cv'):
            if equal(self.c.y, 0): self.set('cv', 0)
            if equal(self.c.v, 0): self.set('cy', 0)
        if ex(self, 'cx id fd'):
            if zero(self.c.x) and self.i.d != self.f.d:
                self.set('ivx', 0)
                self.set('fvx', 0)
        # Direction calculations - horizontal
        if ex(self, 'cd'):
            if self.c.d == 'none':
                self.set('ivx', 0)
        if ex(self, 'ivx'):
            if zero(self.i.vx):
                self.set('cd', 'none')
            elif self.i.vx > 0:
                self.set('cd', 'right')
            elif self.i.vx < 0:
                self.set('cd', 'left')
        if ex(self, 'ia'):
            m = norm(self.i.a-pi/2)
            if equal(self.i.v, 0) or equal(mod(self.i.a-pi/2, pi), 0):
                self.set('cd', 'none')
            elif m > 0 and self.i.v != n:
                self.set('cd', 'left')
            elif m < 0 and self.i.v != n:
                self.set('cd', 'right')
        if ex(self, 'fa'):
            m = norm(self.f.a-pi/2)
            if equal(self.f.v, 0) or equal(mod(self.f.a-pi/2, pi), 0):
                self.set('cd', 'none')
            elif m > 0 and self.f.v != n:
                self.set('cd', 'left')
            elif m < 0 and self.f.v != n:
                self.set('cd', 'right')
        # Direction calculations - vertical
        if ex(self, 'id'):
            if self.i.d == 'none':
                self.set('ivy', 0)
        if ex(self, 'fd'):
            if self.f.d == 'none':
                self.set('fvy', 0)
        if ex(self, 'ivy'):
            if zero(self.i.vy):
                self.set('id', 'none')
            elif self.i.vy > 0:
                self.set('id', 'up')
            elif self.i.vy < 0:
                self.set('id', 'down')
        if ex(self, 'fvy'):
            if zero(self.f.vy):
                self.set('fd', 'none')
            elif self.f.vy > 0:
                self.set('fd', 'up')
            elif self.f.vy < 0:
                self.set('fd', 'down')
        if ex(self, 'ia'):
            if equal(mod(self.i.a, pi), 0):
                self.set('id', 'none')
            else:
                m = norm(self.i.a)
                if m > 0:
                    self.set('id', 'up')
                if m < 0:
                    self.set('id', 'down')
        if ex(self, 'fa'):
            if equal(mod(self.f.a, pi), 0):
                self.set('fd', 'none')
            else:
                m = norm(self.f.a)
                if m > 0:
                    self.set('fd', 'up')
                if m < 0:
                    self.set('fd', 'down')
        # Vertical motion only
        if ex(self, 'cx ct'):
            if zero(self.c.x) and not zero(self.c.t):
                self.set('ivx', 0)
                self.set('fvx', 0)
        # Direction none to velocity
        if ex(self, 'cd id'):
            if self.i.d == 'none' and self.c.d == 'right':
                self.set('ia', 0)
            if self.i.d == 'none' and self.c.d == 'left':
                self.set('ia', pi)
            if self.i.d == 'none' and self.c.d == 'none':
                self.set('ia', 0)
        if ex(self, 'cd fd'):
            if self.f.d == 'none' and self.c.d == 'right':
                self.set('fa', 0)
            if self.f.d == 'none' and self.c.d == 'left':
                self.set('fa', pi)
            if self.f.d == 'none' and self.c.d == 'none':
                self.set('fa', 0)
        if ex(self, 'cd id'):
            if self.c.d == 'none' and self.i.d == 'up':
                self.set('ia', pi/2)
            if self.c.d == 'none' and self.i.d == 'down':
                self.set('ia', -pi/2)
            if self.c.d == 'none' and self.i.d == 'none':
                self.set('ia', 0)
        if ex(self, 'cd fd'):
            if self.c.d == 'none' and self.f.d == 'up':
                self.set('fa', pi/2)
            if self.c.d == 'none' and self.f.d == 'down':
                self.set('fa', -pi/2)
            if self.c.d == 'none' and self.f.d == 'none':
                self.set('fa', 0)
        # Angle calculations
        if ex(self, 'iv'):
            if zero(self.i.v):
                self.set('ivx', 0)
                self.set('ivy', 0)
        if ex(self, 'fv'):
            if zero(self.f.v):
                self.set('fvx', 0)
                self.set('fvy', 0)
        if ex(self, 'iv ia'):
            self.set('ivx', self.i.v * cos(self.i.a))
            self.set('ivy', self.i.v * sin(self.i.a))
        if ex(self, 'ivx ia'):
            if not zero(cos(self.i.a)):
                self.set('iv', self.i.vx / cos(self.i.a))
        if ex(self, 'ivy ia'):
            if not zero(sin(self.i.a)):
                self.set('iv', self.i.vy / sin(self.i.a))
        if ex(self, 'ivx ivy'):
            self.set('iv', sqrt(self.i.vx ** 2 + self.i.vy ** 2))
            self.set('ia', atan2(self.i.vy, self.i.vx))
        if ex(self, 'iv ivx'):
            if not zero(self.i.v):
                try:
                    olda = self.i.a
                    newa = abs(acos(self.i.vx / self.i.v))
                    if zero(newa):
                        self.set('ia', 0)
                    elif self.i.d == 'up':
                        self.set('ia', newa)
                    elif self.i.d == 'down':
                        self.set('ia', -newa)
                except ValueError:
                    self.invalid()
        if ex(self, 'iv ivy'):
            if not zero(self.i.v):
                try:
                    newa = rightq(asin(self.i.vy / self.i.v))
                    if self.c.d in ('none', 'right'):
                        self.set('ia', newa)
                    elif self.c.d == 'left':
                        self.set('ia', otherq(newa))
                except ValueError:
                    self.invalid()
        if ex(self, 'fv fa'):
            self.set('fvx', self.f.v * cos(self.f.a))
            self.set('fvy', self.f.v * sin(self.f.a))
        if ex(self, 'fvx fa'):
            if not zero(cos(self.f.a)):
                self.set('fv', self.f.vx / cos(self.f.a))
        if ex(self, 'fvy fa'):
            if not zero(sin(self.f.a)):
                self.set('fv', self.f.vy / sin(self.f.a))
        if ex(self, 'fvx fvy'):
            self.set('fv', sqrt(self.f.vx ** 2 + self.f.vy ** 2))
            self.set('fa', atan2(self.f.vy, self.f.vx))
        if ex(self, 'fv fvx'):
            if not zero(self.f.v):
                try:
                    newa = abs(acos(self.f.vx / self.f.v))
                    if zero(newa):
                        self.set('fa', 0)
                    elif self.f.d == 'up':
                        self.set('fa', newa)
                    elif self.f.d == 'down':
                        self.set('fa', -newa)
                except ValueError:
                    self.invalid()
        if ex(self, 'fv fvy'):
            if not zero(self.f.v):
                try:
                    newa = rightq(asin(self.f.vy / self.f.v))
                    if self.c.d in ('none', 'right'):
                        self.set('fa', newa)
                    elif self.c.d == 'left':
                        self.set('fa', otherq(newa))
                except ValueError:
                    self.invalid()
        # X position is x = vt
        if ex(self, 'cx ct'):
            if not zero(self.c.t):
                vx = self.c.x / self.c.t
                self.set('ivx', vx)
                self.set('fvx', vx)
        if ex(self, 'ivx cx'):
            if not zero(self.i.vx):
                self.set('ct', self.c.x / self.i.vx)
        if ex(self, 'ivx ct'):
            self.set('cx', self.i.vx * self.c.t)
        # Y acceleration is -g
        if ex(self, 'ct'):
            self.set('cvy', self.c.t * g)
        if ex(self, 'cvy'):
            self.set('ct', self.c.vy / g)
        # y = vt - 0.5gt^2
        if ex(self, 'ivy ct'):
            self.set('cy', self.i.vy * self.c.t + 0.5 * g * self.c.t ** 2)
        # v^2 = v0^2 + 2gy
        if ex(self, 'ivy cy'):
            try:
                newvy = sqrt(self.i.vy ** 2 + 2 * g * self.c.y)
                if zero(newvy):
                    self.set('fvy', 0)
                elif self.f.d == 'up':
                    self.set('fvy', newvy)
                elif self.f.d == 'down':
                    self.set('fvy', -newvy)
            except ValueError:
                self.invalid()
        # complex angles
        if ex(self, 'cx cy ia') and not equal(self.c.x, 0):
            try:
                newt = sqrt(2*(self.c.y-self.c.x*tan(self.i.a))/g)
                if zero(newt):
                    self.set('ct', 0)
                else:
                    try:
                        self.set('ct', newt)
                    except AssertionError:
                        self.set('ct', -newt)
            except ValueError:
                self.invalid()
        if ex(self, 'cx cy fa') and not equal(self.c.x, 0):
            try:
                newt = sqrt(2*(-self.c.y-self.c.x*tan(pi-self.f.a))/g)
                if zero(newt):
                    self.set('ct', 0)
                else:
                    try:
                        self.set('ct', newt)
                    except AssertionError:
                        self.set('ct', -newt)
            except ValueError:
                self.invalid()
        if ex(self, 'cx cy ct') and not zero(self.c.x):
            try:
                newa = atan((self.c.y-0.5*g*self.c.t**2)/self.c.x)
                if self.c.d == 'right':
                    self.set('ia', newa)
                if self.c.d == 'left':
                    self.set('ia', norm(newa-pi))
            except ValueError:
                self.invalid()
    def set(self, var, value): # grapple
        try:
            name = 'self.' + var[0] + '.' + var[1:]
            if var in ('ia', 'fa'):
                value = norm(value)
            if eval(name) == n:
                exec(name + ' = value')
            else:
                exec('assert equal(' + name + ', value)')
        except AssertionError:
            if self.bother:
                if DEBUG:
                    print(name, eval(name), value, equal(eval(name), value))
                    print('Invalidating.')
                self.invalid()
    def overset(self, var, value):
        exec('self.' + var[0] + '.' + var[1:] + ' = value')
    def solve(self):
        try:
            for i in range(10):
                self.deduce()
        except GeneratorExit:
            return
    def qsolve(self):
        try:
            self.set('cvx', 0)
            for var in ('x', 'y', 'v', 'vx', 'vy', 't', 'a'):
                if ex(self, 'i'+var+' f'+var):
                    self.set('c'+var, eval('self.f.'+var+' - self.i.'+var))
                if ex(self, 'i'+var+' c'+var):
                    self.set('f'+var, eval('self.i.'+var+' + self.c.'+var))
                if ex(self, 'f'+var+' c'+var):
                    self.set('i'+var, eval('self.f.'+var+' - self.c.'+var))
            if ex(self, 'cx ivx'):
                self.set('ct', self.c.x / self.i.vx)
            if ex(self, 'cvy'):
                self.set('ct', self.c.vy / g)
            if ex(self, 'ivx ct'):
                self.set('cx', self.i.vx * self.c.t)
            if ex(self, 'ivy ct'):
                self.set('cy', self.i.vy * self.c.t + 0.5 * g * self.c.t ** 2)
            for var in ('x', 'y', 'v', 'vx', 'vy', 't', 'a'):
                if ex(self, 'i'+var+' f'+var):
                    self.set('c'+var, eval('self.f.'+var+' - self.i.'+var))
                if ex(self, 'i'+var+' c'+var):
                    self.set('f'+var, eval('self.i.'+var+' + self.c.'+var))
                if ex(self, 'f'+var+' c'+var):
                    self.set('i'+var, eval('self.f.'+var+' - self.c.'+var))
        except GeneratorExit:
            eval('lolcat. you cannot possibly comprehend.')
            return
    def invalid(self):
        if self.bother:
            global CRASH
            if CRASH: assert False
            global INVALID
            INVALID = True
            raise GeneratorExit
    def plot(self, win):
        self.lines = []
        b = 30
        win1 = (10+b, 10+b)
        win2 = (730-b, 315-b)
        twin1 = (10+b, 315-b)
        twin2 = (730-b, 10+b)
        start = [self.i.x, self.i.y]
        end = [self.f.x, self.f.y]
        if start[0] == n:
            start[0] = 0
            end[0] = self.c.x
        if start[1] == n:
            start[1] = 0
            end[1] = self.c.y
        if equal(start[0], end[0]): return
        new1 = Pointset(i=dc(self.i), f=Point(vy=0), c=Point(), bother=False)
        new1.set('ix', start[0])
        new1.set('iy', start[1])
        new1.qsolve()
        skipeak = False
        try:
            peak = (new1.i.x + new1.c.x, new1.i.y + new1.c.y)
            if not between(peak[0], start[0], end[0]):
                skipeak = True
                new1 = Pointset(i=dc(self.i), f=Point(x=start[0]+self.c.x/2), c=Point(), bother=True)
                new1.set('ix', start[0])
                new1.set('iy', start[1])
                new1.qsolve()
                peak = (new1.i.x + new1.c.x, new1.i.y + new1.c.y)
            func1 = (min(start[0], end[0], peak[0]), min(start[1], end[1], peak[1]))
            func2 = (max(start[0], end[0], peak[0]), max(start[1], end[1], peak[1]))
        except TypeError:
            return
        start, end, peak = (transform(func1, func2, twin1, twin2, start),
                            transform(func1, func2, twin1, twin2, end),
                            transform(func1, func2, twin1, twin2, peak))
        quad = Quadratic().fit(start[0], start[1], end[0], end[1], peak[0], peak[1])
        self.lines = []
        diff = 5
        for x in range(win1[0], win2[0]+diff, diff):
            if x != win1[0]:
                oldy = newy
            newy = quad.f(x)
            if x != win1[0]:
                self.lines.append(Line(Poing(x-diff, oldy), Poing(x, newy)))
                self.lines[-1]._reconfig('width', 3)
                self.lines[-1].draw(win)
        for p in (start, end, peak):
            if p == peak and skipeak: continue
            self.lines.append(Circle(Poing(p[0], p[1]), 5))
            self.lines[-1].setFill({start:'green',end:'red',peak:'blue'}[p])
            self.lines[-1].setOutline({start:'green',end:'red',peak:'blue'}[p])
            self.lines[-1].draw(win)

class Point():
    def __init__(self, x=n, y=n, v=n, vx=n, vy=n, t=n, a=n, d=n):
        self.x = x
        self.y = y
        self.v = v
        self.vx = vx
        self.vy = vy
        self.t = t
        self.a = a
        self.d = d

def mmult(a, b):
    assert len(a[0]) == len(b)
    c = []
    for i in range(len(a)):
        c.append([])
        for j in range(len(b[0])):
            s = 0
            for x in range(len(a[i])):
                s += a[i][x]*b[x][j]
            c[i].append(s)
    return c

def minv3(m1):
    a,b,c,d,e,f,g,h,i = m1[0][0],m1[0][1],m1[0][2],m1[1][0],m1[1][1],m1[1][2],m1[2][0],m1[2][1],m1[2][2]
    m2 = [[e*i-f*h,c*h-b*i,b*f-c*e],[f*g-d*i,a*i-c*g,c*d-a*f],[d*h-e*g,b*g-a*h,a*e-b*d]]
    coef = 1/(a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g))
    for x in range(len(m2)):
        for y in range(len(m2[x])):
            m2[x][y] *= coef
    return m2

class Quadratic():
    def __init__(self, a=n, b=n, c=n):
        self.a = a
        self.b = b
        self.c = c
    def fit(self, x1, y1, x2, y2, x3, y3):
        m1 = [[x1**2,x1,1],[x2**2,x2,1],[x3**2,x3,1]]
        m2 = [[y1],[y2],[y3]]
        m3 = mmult(minv3(m1),m2)
        self.a, self.b, self.c = m3[0][0], m3[1][0], m3[2][0]
        return self
    def f(self, x):
        return self.a*x**2+self.b*x+self.c
    def s(self, y):
        d = sqrt(self.b**2-4*self.a*self.c)
        return (-self.b+d)/self.a, (-self.b-d)/self.a
    def derivative(self):
        return Linear(a/2, b)

def launch():
    global win, fields, checks, labels, svars, cbutton, invalidValues, ps
    global STOP
    win = GraphWin('Projectile Simulation 2013 - Radon Rosborough', 800, 600, autoflush=False)
    fields, checks, labels, svars = [], [], [], []
    names = ['Initial x (m)', 'Initial y (m)', 'Initial v (m/s)', 'Initial vx (m/s)', 'Initial vy (m/s)', 'Initial time (s)', 'Initial angle (deg)', 'Final x (m)', 'Final y (m)', 'Final v (m/s)', 'Final vx (m/s)', 'Final vy (m/s)', 'Final time (s)', 'Final angle (deg)', 'Change in x (m)', 'Change in y (m)', 'Change in v (m/s)', 'Change in vx (m/s)', 'Change in vy (m/s)', 'Change in time (s)', 'Change in angle (deg)']
    for y in range(3):
        for x in range(7):
            fields.append(Entry(Poing(85+x*110, 425+y*75), 10))
            i = len(fields)-1
            fields[-1].text.trace('w', lambda name, index, mode, var=fields[i].text, i=i:checkUpdate(var, i))
            labels.append(Text(Poing(85+x*110, 405+y*75), names[y*7+x]))
            checks.append(Checkbox(Poing(35+x*110, 425+y*75), False))

    fields.append(Entry(Poing(85+0*110, 425+-1*75), 10))
    i = len(fields)-1
    fields[-1].text.trace('w', lambda name, index, mode, var=fields[i].text, i=i:checkUpdate(var, i))
    labels.append(Text(Poing(85+0*110, 405+-1*75), 'Initial direction (up/down/none)'))
    checks.append(Checkbox(Poing(35+0*110, 425+-1*75), False))

    fields.append(Entry(Poing(85+6*110, 425+-1*75), 10))
    i = len(fields)-1
    fields[-1].text.trace('w', lambda name, index, mode, var=fields[i].text, i=i:checkUpdate(var, i))
    labels.append(Text(Poing(85+6*110-25, 405+-1*75), 'Final direction (up/down/none)'))
    checks.append(Checkbox(Poing(35+6*110, 425+-1*75), False))

    fields.append(Entry(Poing(85+3*110, 425+-1*75), 10))
    i = len(fields)-1
    fields[-1].text.trace('w', lambda name, index, mode, var=fields[i].text, i=i:checkUpdate(var, i))
    labels.append(Text(Poing(85+3*110, 405+-1*75), 'Horizontal direction (left/right/none)'))
    checks.append(Checkbox(Poing(35+3*110, 425+-1*75), False))

    for field in fields:
        field.draw(win)
        field.entry.configure(state='readonly')
    fields[-1].entry.configure(state='normal')
    for label in labels:
        label.draw(win)
    for check in checks:
        check.draw(win)

    # Warning: if any more fields are added, this will crash!
    checks[-1].checkon(win)
    STOP = True
    fields[-1].setText('right')
    STOP = False

    for field in fields:
        field.entry.bind('<space>', lambda event, field=field: enableField(field))

    cbutton = Button(Poing(740, 10), Poing(790, 40), 'Close')
    cbutton.draw(win)

    rbutton = Button(Poing(740, 50), Poing(790, 80), 'Reset')
    rbutton.draw(win)

    invalidValues = Text(Poing(400, 150), 'Error: Contradictary values')
    invalidValues.setSize(36)

    ps = Pointset(i=Point(), f=Point(), c=Point())

    checkUpdate(n, n)

    Line(Poing(10, 10), Poing(730, 10)).draw(win)
    Line(Poing(10, 10), Poing(10, 315)).draw(win)
    Line(Poing(730, 10), Poing(730, 315)).draw(win)
    Line(Poing(10, 315), Poing(730, 315)).draw(win)

    win.flush()

    while True:
        click = win.checkMouse()
        if click:
            if cbutton.clicked(click):
                break
            if rbutton.clicked(click):
                for check in checks:
                    check.checkoff(win)
                STOP = True
                for field in fields:
                    field.setText('')
                    field.entry.configure(state='readonly')
                STOP = False
                checkUpdate(n, n)
            for check in checks:
                if (fields[checks.index(check)].getText() == '' or check.checked) and check.check(click, win):
                    if check.checked:
                        state = 'normal'
                        fields[checks.index(check)].entry.focus_set()
                    else:
                        state = 'readonly'
                        fields[checks.index(check)].setText('')
                    fields[checks.index(check)].entry.configure(state=state)
                    win.flush()
    win.close()

if __name__ == '__main__':
    if PROFILE:
        cProfile.run('launch()')
    else:
        launch()
