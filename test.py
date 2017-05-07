import numpy as np
import matplotlib.pyplot as plt
import ujson, time
import csv
from collections import namedtuple

State = namedtuple('State', ['x', 'y', 'theta'])
ExtendedState = namedtuple('ExtendedState', ['x', 'y', 'theta', 'throttle']) # throttle: 1 - forwards, -1: backwards
Action = namedtuple('Action', ['steering_angle', 'throttle'])

def make_arrow(state, length=1.75):
	return [state.x, state.y, np.cos(state.theta)*length, np.sin(state.theta)*length]

class MotionModel(object):
 	""" Base class for describing car dynamics """
 	def __init__(self):
 		pass
 	
 	def apply(self, action, state, step):
 		raise NotImplementedError("Must be defined by child class")

 	def make_curves(self, state, action, step, disc, scale_step=False):
 		xs = []
		ys = []
		xs.append(state.x)
		ys.append(state.y)
		steps = np.linspace(step,0.0,disc,endpoint=False)[::-1]
		for i in steps:
			new_state = self.apply(action, state, i, scale_step=scale_step)
			xs.append(new_state.x)
			ys.append(new_state.y)
		return xs, ys

class EmpericalModel(MotionModel):
 	""" Emperical model based off measured car data """
 	def __init__(self, world_scale=1.0):
 		self.max_steer = 0.335
 		self.L = 0.23
 		self.ackerman_transition = [0.095, 0.105]
 		self.polynomial_coeffs = [-4.035, 5.153, -0.018]
 		self.poly = lambda a: 1.0 / (a*a*self.polynomial_coeffs[0] + a*self.polynomial_coeffs[1] + self.polynomial_coeffs[2])

 		self._EPSILON = 0.0000001
 		self._SCALE = world_scale

 	# only accepts positive, nonzero steering angles
 	def steering_arc_radius(self, angle):
 		unscaled_radius = 0.0
 		if angle <= self.ackerman_transition[0]:
 			unscaled_radius = self.L / np.tan(angle)
 		elif angle >= self.ackerman_transition[1]:
 			unscaled_radius = self.poly(min(angle, self.max_steer))
 		else:
 			width = self.ackerman_transition[1] - self.ackerman_transition[0]
 			t = (angle - self.ackerman_transition[0]) / width
 			ackL  = self.L / np.tan(angle)
 			polyL = self.poly(angle)
 			unscaled_radius = (1.0 - t) * ackL + t * polyL

 		if unscaled_radius == 0.0:
	 		print "THIS SHOULD NOT HAPPEN"
 		return unscaled_radius * self._SCALE

 	def apply(self, action, state, step, scale_step=False):
 		if scale_step:
 			step *= self._SCALE

 		if abs(action.steering_angle) < self._EPSILON:
 			dt = 0
 			dy = 0
 			dx = action.throttle * step
 		else:
 			sign = np.sign(action.steering_angle)
 			R = self.steering_arc_radius(abs(action.steering_angle))
 			dt = sign * action.throttle * step / R
			dy = sign * (1.0-np.cos(dt))*R
			dx = np.sin(dt)*R

 		c, s = np.cos(state.theta), np.sin(state.theta)
		dx_global = c*dx - s*dy
		dy_global = s*dx + c*dy

		return State(dx_global + state.x, dy_global + state.y, (dt + state.theta)%(2.0*np.pi))

class AckermannModel(MotionModel):
 	""" Standard kinematic bicycle motion model """
 	def __init__(self, world_scale=1.0):
 		self._EPSILON = 0.0000001
 		self._SCALE = world_scale

 		self.L = 0.23
 		self.L *= world_scale

 	def apply(self, action, state, step, scale_step=False):
 		if scale_step:
 			step *= self._SCALE

 		if abs(action.steering_angle) < self._EPSILON:
 			dt = 0
 			dy = 0
 			dx = action.throttle * step
 		else:
 			tanangle = np.tan(action.steering_angle)
 			dt = action.throttle * step * tanangle / self.L
			dy = (1.0-np.cos(dt))*self.L / tanangle
			dx = self.L * np.sin(dt) / tanangle

 		c, s = np.cos(state.theta), np.sin(state.theta)
		dx_global = c*dx - s*dy
		dy_global = s*dx + c*dy

		return State(dx_global + state.x, dy_global + state.y, (dt + state.theta)%(2.0*np.pi))

def make_curves(state, action, step, disc):
	xs = []
	ys = []
	xs.append(state.x)
	ys.append(state.y)
	steps = np.linspace(step,0.0,disc,endpoint=False)[::-1]
	for i in steps:
		new_state = motion_model(action, state,step_size=i)
		xs.append(new_state.x)
		ys.append(new_state.y)
	return xs, ys

wheelbase = 4.7619
STEP = 2.0

def unscale_t(t,disc):
	return disc * t / (2.0*np.pi)

def rescale_t(ind,disc):
	return (2.0*np.pi) * ind / float(disc)

def motion_model(action, state, step_size=1.0):
	tanangle = np.tan(action.steering_angle)
	if tanangle == 0:
		dt = 0
		dy = 0
		dx = action.throttle * step_size
	else:
		dt = action.throttle * step_size * tanangle / wheelbase
		dy = (1.0-np.cos(dt))*wheelbase / tanangle
		dx = wheelbase * np.sin(dt) / tanangle

	c, s = np.cos(state.theta), np.sin(state.theta)
	dx_global = c*dx - s*dy
	dy_global = s*dx + c*dy

	return State(dx_global + state.x, dy_global + state.y, (dt + state.theta)%(2.0*np.pi))

class NHPolicy(object):
	"""docstring for NHPolicy"""
	def __init__(self, path):
		print "Loading Policy:", path
		self.path = path
		policy_file = open(path, 'r')
		print "..loading json"
		policy_raw = ujson.load(policy_file)

		if not "policy" in policy_raw:
			print "Incorrectly formatted data, exiting."
			return

		policy_raw = policy_raw["policy"]
		print "..parsing"

		self.d0 = int(policy_raw["d0"])
		self.d1 = int(policy_raw["d1"])
		self.d2 = int(policy_raw["d2"])

		print self.d0, self.d1, self.d2

		self.actions = map(lambda x: Action(float(x["a"]), int(x["t"])), policy_raw["data"])

	def at(self, state, interpolate=False):
		x = int(np.round(state.x))
		y = int(np.round(state.y))
		z = int(np.round(state.theta) % self.d2)
		return self.actions[x+self.d0*(y+self.d1*z)]

class NHPolicy2(object):
	"""docstring for NHPolicy"""
	def __init__(self, path):
		print "Loading Policy:", path
		self.path = path
		summary_file = open(path + "/summary.json", 'r')
		print "..loading json"
		summary_raw = ujson.load(summary_file)

		if not "policy" in summary_raw:
			print "Incorrectly formatted data, exiting."
			return

		summary_raw = summary_raw["policy"]
		print "..parsing summary"

		self.d0 = int(summary_raw["d0"])
		self.d1 = int(summary_raw["d1"])
		self.d2 = int(summary_raw["d2"])

		print self.d0, self.d1, self.d2

		self.backwards_actions = np.zeros((self.d0 * self.d1 * self.d2, 2))
		self.forwards_actions = np.zeros((self.d0 * self.d1 * self.d2, 2))
		
		print "..loading forwards policy"
		f_file = open(path + "/policy_forwards.object", 'r')
	
		i = 0
		for row in csv.reader(f_file, delimiter=' ', quotechar='|'):
			self.forwards_actions[i, :] = row
			i += 1

		print "..loading backwards policy"
		b_file = open(path + "/policy_forwards.object", 'r')
		i = 0
		for row in csv.reader(b_file, delimiter=' ', quotechar='|'):
			self.backwards_actions[i, :] = row
			i += 1

		f_file.close()
		b_file.close()
		summary_file.close()

	def at(self, state, interpolate=False):
		x = int(np.round(state.x))
		y = int(np.round(state.y))
		z = int(np.round(state.theta) % self.d2)
		if state.throttle > 0:
			action = self.forwards_actions[x+self.d0*(y+self.d1*z)]
		else:
			action = self.backwards_actions[x+self.d0*(y+self.d1*z)]
		return Action(action[0], action[1])

class SimulatorExtendedState(object):
	""" Simulates using the given optimal policy """
	def __init__(self, Policy, Cost=None):
		self.Policy = Policy
		self.Cost = Cost
		self.debug = True
		self.has_image = False

	def load_image(self, path):
		self.im = plt.imread(path)
		self.has_image = True

	def n_simulation_steps(self, n, start_state, make_history_array=False):
		state = start_state
		if make_history_array:
			state_hist = [state]
			action_hist = []
		
		for i in xrange(n):
			state, action = self.simulate(state, give_action=True)
			if make_history_array:
				state_hist.append(state)
				action_hist.append(action)
		
		if make_history_array:
			return state, state_hist, action_hist
		else:
			return state

	# run a single simulation step based on the transition function and optimal policy
	def simulate(self, start_state, interpolate=False, give_action=False):
		if start_state.throttle == 0:
			# print "STOP"
			optimal_control = self.Policy.at(start_state, interpolate=interpolate)
			new_state = start_state
		else:
			optimal_control = self.Policy.at(start_state, interpolate=interpolate)
			new_state = self.transition(start_state, optimal_control)

		if not optimal_control.throttle == self.last_throttle:
			print "REVERSAL"
			self.last_throttle = optimal_control.throttle

		if give_action:
			return new_state, optimal_control
		else:
			return new_state

	def transition(self, state, action):
		# print
		# print state.theta
		unscaled_state = State(state.x,state.y, rescale_t(state.theta, self.Policy.d2))
		# print unscaled_state.theta
		ns = motion_model(action, unscaled_state, step_size=STEP)
		# print ns.theta
		scaled_new_state = ExtendedState(ns.x,ns.y, unscale_t(ns.theta, self.Policy.d2), action.throttle)
		# print scaled_new_state.theta
		return scaled_new_state

	def on_keydown(self, evt):
		self.down_spot = evt.xdata, evt.ydata

	def on_keyup(self, evt):
		up_spot = evt.xdata, evt.ydata

		start_x = self.down_spot[0]
		start_y = self.down_spot[1]

		deltaY = up_spot[1] - start_y
		deltaX = up_spot[0] - start_x

		angle = np.arctan2(deltaY, deltaX)
		print self.down_spot, "-->", up_spot, angle, unscale_t(angle, self.Policy.d2)
		self.do_the_thing(State(start_x, start_y, unscale_t(angle, self.Policy.d2)))
		
	def do_the_thing(self, state):
		self.last_throttle = 1
		start = ExtendedState(state.x,state.y,state.theta, 1)
		end, state_hist, action_hist = self.n_simulation_steps(200,start, make_history_array=True)
		# print state_hist
		self.ax.cla()

		if self.has_image:
			implot = self.ax.imshow(self.im)

		arrows = []
		line_xs = []
		line_ys = []
		for i in xrange(len(state_hist)):
			state = state_hist[i]
			s = State(state.x, state.y, rescale_t(state.theta, self.Policy.d2))
			arrows.append(make_arrow(s))
			
			if i < len(action_hist):
				lx, ly = make_curves(s, action_hist[i], STEP, 10)
				line_xs += lx
				line_ys += ly
		
		soa = np.array(arrows)
		# soa = np.array([make_arrow(State(1,1,0)),make_arrow(State(1,1,np.pi/2.0)),make_arrow(State(1,1,np.pi)),make_arrow(State(1,1,3.0*np.pi/2.0)),make_arrow(State(1,1,0))])
		X, Y, U, V = zip(*soa)
		plt.plot(line_xs, line_ys)
		self.ax.quiver(X, Y, U, V, angles='xy', scale_units='xy', scale=1)
		self.ax.set_xlim([0, self.Policy.d0])
		self.ax.set_ylim([0, self.Policy.d1])
		# # ax.set_xticks(np.arange(-3, 3, 0.1))
		# # ax.set_yticks(np.arange(-3, 3., 0.1))
		# plt.grid()
		plt.draw()

	def plot(self, state_hist, action_hist):
		arrows = []
		line_xs = []
		line_ys = []
		for i in xrange(len(state_hist)):
			state = state_hist[i]
			s = State(state.x, state.y, rescale_t(state.theta, self.Policy.d2))
			arrows.append(make_arrow(s))
			
			if i < len(action_hist):
				lx, ly = make_curves(s, action_hist[i], STEP, 10)
				line_xs += lx
				line_ys += ly
		
		soa = np.array(arrows)
		# soa = np.array([make_arrow(State(1,1,0)),make_arrow(State(1,1,np.pi/2.0)),make_arrow(State(1,1,np.pi)),make_arrow(State(1,1,3.0*np.pi/2.0)),make_arrow(State(1,1,0))])
		X, Y, U, V = zip(*soa)
		plt.figure()
		plt.plot(line_xs, line_ys)
		ax = plt.gca()


		if self.has_image:
			implot = ax.imshow(self.im)

		ax.quiver(X, Y, U, V, angles='xy', scale_units='xy', scale=1)
		ax.set_xlim([0, self.Policy.d0])
		ax.set_ylim([0, self.Policy.d1])
		# ax.set_xticks(np.arange(-3, 3, 0.1))
		# ax.set_yticks(np.arange(-3, 3., 0.1))
		plt.grid()
		plt.draw()

	def interactive_plot(self):
		fig = plt.figure()
		fig.canvas.mpl_connect('button_press_event', self.on_keydown)
		fig.canvas.mpl_connect('button_release_event', self.on_keyup)

		self.ax = plt.gca()
		if self.has_image:
			implot = self.ax.imshow(self.im)
		self.ax.set_xlim([0, self.Policy.d0])
		self.ax.set_ylim([0, self.Policy.d1])
		plt.show()

# if __name__ == '__main__':
# 	t = time.time()
# 	# policy = NHPolicy2("./build/serialized/parallel_park")
# 	policy = NHPolicy2("./build/serialized/parallel_park_basement_left")
# 	print "Loaded in:", time.time() - t

# 	s = ExtendedState(0,0,0,0)
# 	print policy.at(s)
# 	sim = SimulatorExtendedState(policy)

# 	# sim.load_image("70x70_map.png")
# 	# sim.load_image("./maps/210x210_kitchen.png")
# 	sim.load_image("./maps/parking_goal_map_left.png")
# 	# sim.load_image("160x160_parallel2.png")
# 	sim.interactive_plot()
# 	print "DONE"

if __name__ == '__main__':
	world_scale = 1.0
	max_angle = 0.35
	disc = 19
	steering_angles = np.linspace(max_angle, -max_angle, num=disc)
	speeds = [-1, 1]
	actions = []

	emp = EmpericalModel(world_scale=world_scale)
	ack = AckermannModel(world_scale=world_scale)

	# sa = np.linspace(0.05, 0.15, num=18)
	# for a in sa:
	# 	print a, "-->", model.steering_arc_radius(a)

	for a in steering_angles:
		actions.append(Action(a, speeds[0]))
		actions.append(Action(a, speeds[1]))

	state1 = State(-1.25,0,np.pi/2.0)
	state2 = State(1.25,0,np.pi/2.0)

	step = 1.75

	for a in actions:
		print a

		lx, ly = ack.make_curves(state1, a, step, 10, scale_step=True)
		plt.plot(lx, ly, color="red")

		lx, ly = emp.make_curves(state2, a, step, 10, scale_step=True)
		plt.plot(lx, ly, color="black")

	mv = 2.0 * world_scale
	plt.xlim([-1.5*mv, 1.5*mv])
	plt.ylim([-mv, mv])
	plt.show()
