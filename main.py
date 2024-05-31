import math

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import figure
from flytmath.vector import Vector
from flytmath.polygon import Polygon
def plot_polygon(polygon, sweep_dir, start, sweep_direction, ls):
    # Extract x and y coordinates from vertices
    #print(polygon)

    x = [vertex.x for vertex in polygon]
    y = [vertex.y for vertex in polygon]
    path = calc_path(polygon, sweep_dir, start, sweep_direction, ls)
    X = [point[0] for point in path]
    Y = [point[1] for point in path]
    # print(X)
    # print(Y)
    plt.plot(x + [x[0]], y + [y[0]], 'r-')  # Connect last point to the first to close the polygon
    plt.plot(X + [X[0]], Y + [Y[0]], 'b-')  # Connect last point to the first to close the polygon
    # Plot polygon
    #  plt.plot(x + [x[0]], y + [y[0]], 'r-')  # Connect last point to the first to close the polygon
    # plt.arrow(35,35, sweep_dir[0],sweep_dir[1], head_width=0.1, head_length=0.1)
    # Plot start and end points
    plt.plot(start[0], start[1], 'go')  # Red circle for start point
    # Set axis limits
    #plt.xlim(min(x) - 1, max(x) + 1)
    #plt.ylim(min(y) - 1, max(y) + 1)

    # Add labels and title
    #plt.xlabel('X')
    #plt.ylabel('Y')
    #plt.title('Polygon with Start and End Points')

    # Show plot
    #plt.grid(True)
    #plt.gca().set_aspect('equal', adjustable='box')
    #plt.show()


# Example usage
def perp_dist(point, line_vector, line_start):
    line_end = Polygon.addvect(line_start , line_vector)

    # Vector from line start to the point
    point_vector = Polygon.subvect(point , line_start)

    # Calculate the projection of point_vector onto line_vector
    projection_length = np.dot(Vector.ret_arr(point_vector),Vector.ret_arr(line_vector)) / np.dot(Vector.ret_arr(line_vector), Vector.ret_arr(line_vector))

    # Calculate the projection point on the line
    projection_point = Polygon.addvect(line_start , Vector.__mul__(line_vector, projection_length))

    # Calculate the perpendicular distance from the point to the line
    opposed_vector = Polygon.subvect(point , projection_point)
    distance = np.linalg.norm(Vector.ret_arr(opposed_vector))

    return distance, opposed_vector


def perpendicular_unit_vector(vector):
    # Find a vector not parallel to the given vector
    if Vector.__getitem__(vector,0) != 0 or Vector.__getitem__(vector,1) != 0:
        perpendicular_vector = Vector(-1*Vector.__getitem__(vector,1), Vector.__getitem__(vector,0))
    else:
        perpendicular_vector = Vector(-1*Vector.__getitem__(vector,0), -1*Vector.__getitem__(vector,2))

    # Normalize the perpendicular vector to obtain a unit vector
    unit_perpendicular_vector = Vector.__truediv__(perpendicular_vector , np.linalg.norm(Vector.ret_arr(perpendicular_vector)))

    return unit_perpendicular_vector


def find_line_sweep(polygon):
    edges = Polygon.edge_polygon(polygon)
    sweep_direction=Vector(0,0)
    optimal = 0
    for j in range(Polygon.__len__(edges)):
        max_dist_edge = 0

        for k in range(Polygon.__len__(polygon)):  # for each point in the polygon,
            perp_dis, opp_vect = perp_dist(Polygon.__getitem__(polygon,k), Polygon.__getitem__(edges,j), Polygon.__getitem__(polygon,j))
            # print(perp_dis)
            # print(opp_vect)
            if (perp_dis > max_dist_edge):
                max_dist_edge = perp_dis
                sweep_direction = opp_vect

        if (max_dist_edge < optimal or j == 0):
            optimal = max_dist_edge
            line_sweep = perpendicular_unit_vector(sweep_direction)
            start_point = Polygon.__getitem__(polygon,j)

    # print(optimal)
    return line_sweep, start_point, sweep_direction, optimal



def euclidean_distance(point1, point2):
    return np.linalg.norm(Vector.ret_arr(Polygon.subvect(point1 , point2)))


def is_point_inside_polygon(point, vertices):
    x, y = point[0], point[1]
    n = len(vertices)
    inside = False
    p1x, p1y = vertices[0]
    for i in range(n + 1):
        p2x, p2y = vertices[i % n]
        if y > min(p1y, p2y):
            if y <= max(p1y, p2y):
                if x <= max(p1x, p2x):
                    if p1y != p2y:
                        xinters = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                    if p1x == p2x or x <= xinters:
                        inside = not inside
        p1x, p1y = p2x, p2y
    return inside


def line_equation(point, unit_vector):
    # Calculate line equation ax + by = c
    a = unit_vector.y
    b = -unit_vector.x
    c = a * point.x + b * point.y
    return a, b, c


def intersection_point(edge_start, edge_end, a, b, c):
    # Calculate intersection point between line and edge
    x1, y1 = edge_start.x, edge_start.y
    x2, y2 = edge_end.x, edge_end.y
    det = a * x1 + b * y1 - c
    det2 = a * x2 + b * y2 - c
    if det * det2 > 0:  # Both points are on the same side of the line
        return None
    elif det == 0 and det2 == 0:  # Edge lies on the line
        return None
    else:
        dx = x2 - x1
        dy = y2 - y1
        denom = a * dx + b * dy
        if denom == 0:  # Parallel edge and line
            return None
        else:
            t = (c - a * x1 - b * y1) / denom
            if 0 <= t <= 1:  # Intersection point lies within the edge
                return Vector(x1 + t * dx, y1 + t * dy)
            else:
                return None


def find_intersection_points(polygon, point, unit_vector):
    # Calculate line equation ax + by = c
    a, b, c = line_equation(point, unit_vector)
    intersections = []
    len = Polygon.__len__(polygon)
    for i in range(len):
        edge_start = Polygon.__getitem__(polygon,i)
        edge_end = Polygon.__getitem__(polygon,((i + 1) % len))
        intersect = intersection_point(edge_start, edge_end, a, b, c)
        if intersect is not None:
            intersections.append(intersect)
    return intersections


def calc_path(polygon, sweep_along, start, sweep_dir, ls):
    foot_x = 3
    foot_y = 5
    overlap_x = 0.2
    overlap_y = 0.4
    d = foot_x * (1 - overlap_x)
    path = []
    start = Polygon.addvect(start , Vector.__mul__( (Vector.normalize(sweep_dir)), (foot_x/2)))
    intersection_points = find_intersection_points(polygon, start, sweep_along)
    if not intersection_points:
        sweep_dir = Vector.__mul__(sweep_dir , -1)
        start = Polygon.addvect(start , Vector.__mul__(Vector.normalize(sweep_dir) ,(foot_x) ) )
    # print(start)
    #print(ls)
    #print(d)
    z = ls - (foot_x / 2)
    turns = 0
    if z % d <= foot_x / 2:
        x = 4 * math.ceil(z / d)
    elif z % d > foot_x / 2:
        x = 4 * (math.ceil(z / d) + 1)
    #print(x)

    while (True ):
        intersection_points = find_intersection_points(polygon, start, sweep_along)
        if not intersection_points:
            break
        # print(start)
        if euclidean_distance(start, intersection_points[0]) < euclidean_distance(start, intersection_points[1]):
            path.append(intersection_points[0])
            path.append(intersection_points[1])
        else:
            path.append(intersection_points[1])
            path.append(intersection_points[0])
        turns = turns + 2
        start = Polygon.addvect(path[-1] , Vector.__mul__((Vector.normalize(sweep_dir)),d))
    #print(turns)

    return path

#vector=Vector(5.3,7.5)
#print(vector)
Center=Vector(35,35)
vertices = [(35,35), (40,35), (40,40), (35,40)]
polygon1= Polygon.from_tuples(vertices)
polygon1=polygon=Polygon.regular(Center, 25, 4)

polygon=Polygon.regular(Center, 85, 7)
polygon.points[6]=Vector(65,0)
print(polygon)
#polygon1.points[0]=Vector(85,35)
#polygon1.points[12]=Vector(27.1,28.3)
print(polygon1)
lista= Polygon.convex_decompose(polygon,[polygon1])
k=1
d=0
figure(figsize=(16, 12), dpi=80)
for poly in lista :
    sweep_along, start_point, sweep_dir, ls = find_line_sweep(poly)
    plot_polygon(poly, sweep_along, start_point, sweep_dir, ls)

    #plt.xlim(min(x) - 1, max(x) + 1)
    #plt.ylim(min(y) - 1, max(y) + 1)

    # Add labels and title
    #plt.xlabel('X')
    #plt.ylabel('Y')
    #plt.title('Polygon with Start and End Points')

    # Show plot
    #plt.grid(True)
    #plt.gca().set_aspect('equal', adjustable='box')
plt.show()

#print(polygon)
#print(yo)

#sweep_along, start_point, sweep_dir, ls = find_line_sweep(polygon)
#plot_polygon(polygon, sweep_along, start_point, sweep_dir, ls)
# print(normalize_vector(sweep_dir)    plt.savefig('2nd.png'))
# print(ls)




