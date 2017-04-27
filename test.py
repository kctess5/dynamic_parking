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

wheelbase = 10.4
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

if __name__ == '__main__':
	t = time.time()
	policy = NHPolicy2("./build/serialized/parallel_park")
	print "Loaded in:", time.time() - t

	s = ExtendedState(0,0,0,0)
	print policy.at(s)
	sim = SimulatorExtendedState(policy)

	# sim.load_image("70x70_map.png")
	sim.load_image("./maps/210x210_kitchen.png")
	# sim.load_image("160x160_parallel2.png")
	sim.interactive_plot()
	print "DONE"