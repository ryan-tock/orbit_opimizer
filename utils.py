import math

NEWTON_ITERATIONS = 6
class Orbit:
    def __init__(self, e, i, node, periapsis, m_0, p, data):
        self.e =  e                #Eccentricity
        self.beta = e / (1 + math.sqrt(1 - e ** 2))
        self.i = i                 #Inclination
        self.node = node           #Longitude of the Ascending Node
        self.periapsis = periapsis #Argument of Periapsis
        self.m_0 = m_0             #Mean anomaly on January 1st, 2000
        self.p = p                 #Period
        self.predicted_positions = []
        sm = self.__calculate_semi_major(data)
        self.sm = sm               # Semi-Major Axis

        self.data = data

    def __calculate_semi_major(self, data):
        parameter_squared = 0
        resultant = 0
        for pos in data:
            t = pos[0]
            x_actual = pos[1]
            y_actual = pos[2]

            mean = self.calculate_mean(t)
            eccentric = self.calculate_eccentric(mean)
            true = self.calculate_true(eccentric)

            pos_predict = self.__calculate_pos_scaled(eccentric, true)
            self.predicted_positions.append(pos_predict)

            parameter_squared += pos_predict[0] ** 2
            parameter_squared += pos_predict[1] ** 2

            resultant += pos_predict[0] * x_actual
            resultant += pos_predict[1] * y_actual

        print(parameter_squared)
        return resultant / parameter_squared

    def calculate_mean(self, t):
        return self.m_0 + 2 * math.pi / self.p * (t - 2000)

    def calculate_eccentric(self, mean_anomaly):
        guess = mean_anomaly

        for i in range(NEWTON_ITERATIONS):
            guess = guess + (mean_anomaly - guess + self.e * math.sin(guess)) / (
                        1 - self.e * math.cos(guess))

        return guess

    def calculate_true(self, eccentric_anomaly):
        return eccentric_anomaly + 2 * math.atan(self.beta * math.sin(eccentric_anomaly) / (1 - self.beta * math.cos(eccentric_anomaly)))


    def __calculate_pos_scaled(self, eccentric_anomaly, true_anomaly):
        radius_scaled = 1 - self.e * math.cos(eccentric_anomaly)
        planar_angle = (math.cos(true_anomaly + self.periapsis), math.sin(true_anomaly + self.periapsis))
        node_angle = (math.cos(self.node - 3 * math.pi / 2), math.sin(self.node - 3 * math.pi / 2))
        inclined_angle = math.cos(self.i)

        x = radius_scaled * (planar_angle[0] * node_angle[0] - inclined_angle * planar_angle[1] * node_angle[1])
        y = radius_scaled * (inclined_angle * planar_angle[1] * node_angle[0] + planar_angle[0] * node_angle[1])

        return x,y

    def calculate_pos(self, eccentric_anomaly, true_anomaly):
        return self.sm * self.__calculate_pos_scaled(eccentric_anomaly, true_anomaly)

    def calculate_error(self):
        error = 0
        for i in range(len(self.data)):
            error += (self.data[i][1] - self.sm * self.predicted_positions[i][0]) ** 2
            error += (self.data[i][2] - self.sm * self.predicted_positions[i][1]) ** 2

        return error



file_data = open("data.csv", 'r').read()
file_data = [[float(j) for j in line.split(",")] for line in file_data.split("\n")]

orbit = Orbit(0.33, 0.66, 3.55, 0.4, 5.97, 40, file_data)
print(orbit.calculate_error())