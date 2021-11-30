import numpy as np
from numpy import rad2deg,deg2rad
import sympy
from sympy import cos, sin, solve
from math import atan2,sqrt

theta1, theta2, theta3, theta4, theta5, theta6,u = sympy.symbols("theta1 theta2 theta3 theta4 theta5 theta6 u")
Td2w = [
        [1, 0, 0, 830],
        [0, 1, 0, 20],
        [0, 0, 1, 330],
        [0, 0, 0, 1]]

theta = deg2rad(35)

To2d = [
        [cos(theta), -sin(theta), 0, -280],
        [sin(theta), cos(theta), 0, 250],
        [0, 0, 1, 62.5],
        [0, 0, 0, 1]]

Td2w = np.array(Td2w)
To2d = np.array(To2d)

To2w = Td2w.dot(To2d)
print(To2w)

T02w = [
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, 373],
        [0, 0, 0, 1]]

To26 = [
        [0, 0, 1, 0],
        [0, -1, 0, 0],
        [1, 0, 0, 206],
        [0, 0, 0, 1]]

T02w = np.array(T02w)
To26 = np.array(To26)
T02w_inv = np.linalg.inv(T02w)
To26_inv =  np.linalg.inv(To26)

T620  = T02w_inv.dot(To2w)
T620 = T620.dot(To26_inv)

print(T620)

px = T620[0][3]
py = T620[1][3]
pz = T620[2][3]

a0 = 0
a1 = -30
a2 = 340
a3 = -40
a4 = 0
a5 = 0

d1 = 0
d2 = 0
d3 = 0
d4 = 338
d5 = 0
d6 = 0


alpha0 = deg2rad(0)
alpha1 = deg2rad(-90)
alpha2 = deg2rad(0)
alpha3 = deg2rad(-90)
alpha4 = deg2rad(90)
alpha5 = deg2rad(-90)

# c1 = cos(theta1);
c2 = cos(theta2)
c3 = cos(theta3)
# c4 = cos(theta4);
# c5 = cos(theta5);
# c6 = cos(theta6);

# s1 = sin(theta1);
s2 = sin(theta2)
s3 = sin(theta3)
# s4 = sin(theta4);
# s5 = sin(theta5);
# s6 = sin(theta6);

f1 = a3 * c3 + d4 * sin(alpha3)*s3 + a2
f2 = a3 * cos(alpha2) * s3 - d4 * sin(alpha3) * cos(alpha2) * c3 - d4 * sin(alpha2) * cos(alpha3) - d3 * sin(alpha2)
f3 = a3 * sin(alpha2) * s3 - d4 * sin(alpha3) * sin(alpha2) * c3 + d4 * cos(alpha2) * cos(alpha3) + d3 * cos(alpha2)

k1 = f1
k2 = -f2
k3 = f1**2 + f2**2 + f3**2 + a1**2 + d2**2 + 2 * d2 * f3
k4 = f3 * cos(alpha1) + d2 * cos(alpha1)

# r = (k1*c2 + k2*s2)*2 *a1 + k3
# z = (k1*s2-k2*c2) * sin(alpha1) + k4

r = px**2 + py**2 + pz**2
z = pz

print(r)
print(z)
equation1 = (r - k3) ** 2 / (4 * (a1 ** 2)) + (z - k4) ** 2 / (sin(alpha1) ** 2) - k1 ** 2 - k2 ** 2
res = (solve(equation1, theta3))
for i in res:
    print("theta3 degree:", rad2deg(float(i)))
