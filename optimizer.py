import math
from scipy.optimize import differential_evolution

NEWTON_ITERATIONS = 6
class Orbit:
    def __init__(self, e, i, node, periapsis, m_0, p):
        self.sm = None              # Semi-Major Axis
        self.e =  e                 # Eccentricity
        self.i = i                  # Inclination
        self.node = node            # Longitude of the Ascending Node
        self.periapsis = periapsis  # Argument of Periapsis
        self.m_0 = m_0              # Mean anomaly on January 1st, 2000
        self.p = p                  # Period

        self.beta = e / (1 + math.sqrt(1 - e ** 2))
        self.predicted_positions = []

    def optimize_semi_major(self, data):
        parameter_squared = 0
        resultant = 0
        for pos in data:
            t = pos[0]
            x_actual = pos[1]
            y_actual = pos[2]

            pos_predict = self.__calculate_pos_scaled(t)
            self.predicted_positions.append(pos_predict)

            parameter_squared += pos_predict[0] ** 2
            parameter_squared += pos_predict[1] ** 2

            resultant += pos_predict[0] * x_actual
            resultant += pos_predict[1] * y_actual

        self.sm = resultant / parameter_squared

    def __calculate_mean(self, t):
        return self.m_0 + 2 * math.pi / self.p * (t - 2000)

    def __calculate_eccentric(self, mean_anomaly):
        guess = mean_anomaly

        for i in range(NEWTON_ITERATIONS):
            guess = guess + (mean_anomaly - guess + self.e * math.sin(guess)) / (
                        1 - self.e * math.cos(guess))

        return guess

    def __calculate_true(self, eccentric_anomaly):
        return eccentric_anomaly + 2 * math.atan(self.beta * math.sin(eccentric_anomaly) / (1 - self.beta * math.cos(eccentric_anomaly)))

    def __calculate_pos_scaled(self, t):
        mean_anomaly = self.__calculate_mean(t)
        eccentric_anomaly = self.__calculate_eccentric(mean_anomaly)
        true_anomaly = self.__calculate_true(eccentric_anomaly)

        radius_scaled = 1 - self.e * math.cos(eccentric_anomaly)
        planar_angle = (math.cos(true_anomaly + self.periapsis), math.sin(true_anomaly + self.periapsis))
        node_angle = (math.cos(self.node - 3 * math.pi / 2), math.sin(self.node - 3 * math.pi / 2))
        inclined_angle = math.cos(self.i)

        x = radius_scaled * (planar_angle[0] * node_angle[0] - inclined_angle * planar_angle[1] * node_angle[1])
        y = radius_scaled * (inclined_angle * planar_angle[1] * node_angle[0] + planar_angle[0] * node_angle[1])

        return x,y

    def calculate_pos(self, t):
        return self.sm * self.__calculate_pos_scaled(t)

    def calculate_error(self, data):
        error = 0
        for i in range(len(data)):
            error += (data[i][1] - self.sm * self.predicted_positions[i][0]) ** 2
            error += (data[i][2] - self.sm * self.predicted_positions[i][1]) ** 2

        if self.sm < 0:
            return error * 100
        return error

    def calculate_r_squared(self, data):
        x_mean = sum([element[1] for element in data]) / len(data)
        y_mean = sum([element[2] for element in data]) / len(data)
        return 1 - (self.calculate_error(data) / sum([(x_mean - element[1]) ** 2 + (y_mean - element[2]) ** 2 for element in data]))

file_data = open("data.csv", 'r').read() # File format must be comma separated t,x,y
file_data = [[float(j) for j in line.split(",")] for line in file_data.split("\n")]

def orbit_error(x, data):
    orbit = Orbit(x[0], x[1], x[2], x[3], x[4], x[5])
    orbit.optimize_semi_major(data)
    return orbit.calculate_error(data)


bounds = [(0, 1), (0, 2 * math.pi), (0, 2 * math.pi), (0, 2 * math.pi), (0, 2 * math.pi), (10, 1000)]
result = differential_evolution(orbit_error, bounds, args=(file_data,))

best_orbit = Orbit(result.x[0], result.x[1], result.x[2], result.x[3], result.x[4], result.x[5])
best_orbit.optimize_semi_major(file_data)

print("Optimal values:")
print("a: ", best_orbit.sm)
print("e: ", best_orbit.e)
print("i: ", best_orbit.i)
print("Ω: ", best_orbit.node)
print("⍵: ", best_orbit.periapsis)
print("M: ", best_orbit.m_0)
print("p: ", best_orbit.p)
print()
print("R^2 value:")
print(best_orbit.calculate_r_squared(file_data))
