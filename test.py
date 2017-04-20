import numpy as np
import matplotlib.pyplot as plt
import ujson

from collections import namedtuple

State = namedtuple('State', ['x', 'y', 'theta'])
ExtendedState = namedtuple('State', ['x', 'y', 'theta', 'throttle']) # throttle: 1 - forwards, -1: backwards
Action = namedtuple('Action', ['steering_angle', 'throttle'])

def make_arrow(state, length=1.75):
	# print [state.x, state.y, state.x+np.cos(state.theta)*length, state.y+np.sin(state.theta)*length], state.theta
	# return [state.x, state.y, state.x+np.cos(state.theta)*length, state.y+np.sin(state.theta)*length]
	return [state.x, state.y, np.cos(state.theta)*length, np.sin(state.theta)*length]

def make_curves(state, action, step, disc):
	xs = []
	ys = []

	xs.append(state.x)
	ys.append(state.y)
	steps = np.linspace(4.0,0.0,4,endpoint=False)[::-1]
	for i in steps:
		new_state = motion_model(action, state,step_size=i)
		xs.append(new_state.x)
		ys.append(new_state.y)
	return xs, ys

wheelbase = 1.0
STEP = 2.5

def unscale_t(t,disc):
	return disc * t / (2.0*np.pi)

def rescale_t(ind,disc):
	return (2.0*np.pi) * ind / float(disc)


#define RESCALE_T(ind, disc) (M_2PI * (state_T)(ind) / (state_T)(disc))
#define UNSCALE_T(val, disc) ((disc) * (val) / M_2PI)

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
		self.actions = map(lambda x: Action(float(x["a"]), int(x["t"])), policy_raw["data"])

	def at(self, state, interpolate=False):
		x = int(np.round(state.x))
		y = int(np.round(state.y))
		z = int(np.round(state.theta) % self.d2)
		return self.actions[x+self.d0*(y+self.d1*z)]

class Cost(object):
	"""docstring for Cost"""
	def __init__(self, path):
		self.d0 = None
		self.d1 = None
		self.d2 = None

class Simulator(object):
	""" Simulates using the given optimal policy """
	def __init__(self, Policy, Cost=None):
		self.Policy = Policy
		self.Cost = Cost
		self.debug = True

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
		if self.debug: print "   at state: ", start_state
		optimal_control = self.Policy.at(start_state, interpolate=interpolate)
		new_state = self.transition(start_state, optimal_control)
		if self.debug: print "     -- control  : ", optimal_control
		if self.debug: print "     -- new state: ", new_state
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
		scaled_new_state = State(ns.x,ns.y, unscale_t(ns.theta, self.Policy.d2))
		# print scaled_new_state.theta
		return scaled_new_state

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
		ax.quiver(X, Y, U, V, angles='xy', scale_units='xy', scale=1)
		ax.set_xlim([0, self.Policy.d0])
		ax.set_ylim([0, self.Policy.d1])
		# ax.set_xticks(np.arange(-3, 3, 0.1))
		# ax.set_yticks(np.arange(-3, 3., 0.1))
		plt.grid()
		plt.draw()
		plt.show()


# class InteractivePlot(object):
# 	def __init__(self):
# 		# self.fig, (self.ax1,self.ax2) = plt.subplots(2, 1)
# 		self.fig, self.ax = plt.subplots()

# 		# self.fig.canvas.setFocusPolicy( QtCore.Qt.ClickFocus )
# 		# self.fig.canvas.setFocus()

# 		# self.fig.canvas.mpl_connect('scroll_event', self.onscroll)
# 		self.fig.canvas.mpl_connect('button_press_event', self.on_keypress)

# 		self.update()

# 	def on_keypress(self, evt):
# 		print "key_press_event"

# 	def onscroll(self, evt):
# 		print "asdfasdf"
# 		self.update()

# 	def update(self):
# 		print "update"
# 		self.ax.plot(np.random.rand(12), np.random.rand(12), 'go')
# 		xl = self.ax.set_xlabel('easy come, easy go')
# 		self.ax.set_title('Press a key')
# 		# plt.tight_layout()
# 		# self.ax1.cla()
# 		# self.ax1.axis('off')

# 		# self.ax1.set_title("Test")
# 		# self.ax1.set_ylabel("TESTEST")
# 		# self.ax1.imshow(np.zeros((100,100)),cmap="gray",interpolation='nearest', aspect='auto')

# 		# self.fig.canvas.draw()


class SimulatorExtendedState(object):
	""" Simulates using the given optimal policy """
	def __init__(self, Policy, PolicyReverse, Cost=None):
		self.Policy = Policy
		self.PolicyReverse = PolicyReverse
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
		if start_state.throttle > 0:
			print "FORWARDS"
			optimal_control = self.Policy.at(start_state, interpolate=interpolate)
			new_state = self.transition(start_state, optimal_control)
		elif start_state.throttle == 0:
			# print "STOP"
			optimal_control = self.Policy.at(start_state, interpolate=interpolate)
			new_state = start_state
		else:
			print "REVERSE"
			optimal_control = self.PolicyReverse.at(start_state, interpolate=interpolate)
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
		end, state_hist, action_hist = self.n_simulation_steps(50,start, make_history_array=True)
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


		
		# plt.plot(line_xs, line_ys)
		self.ax = plt.gca()

		if self.has_image:
			implot = self.ax.imshow(self.im)


		plt.show()

		# ax.quiver(X, Y, U, V, angles='xy', scale_units='xy', scale=1)
		# ax.set_xlim([0, self.Policy.d0])
		# ax.set_ylim([0, self.Policy.d1])
		# # ax.set_xticks(np.arange(-3, 3, 0.1))
		# # ax.set_yticks(np.arange(-3, 3., 0.1))
		# plt.grid()
		# plt.draw()
		# plt.show()


if __name__ == '__main__':
	# InteractivePlot()
	# plt.show()

	# exit()

	policy = NHPolicy("./build/serialized/policy.object")
	# policy = NHPolicy("./build/serialized/go_center_70x70x200.object")
	# sim = Simulator(policy)

	sim = SimulatorExtendedState(policy, NHPolicy("./build/serialized/policy_reverse.object"))
	# sim = SimulatorExtendedState(policy, True)

	sim.load_image("70x70_map.png")
	sim.interactive_plot()

	# for y in xrange(1,70, 4):
	# 	print
	# 	start = State(0,y,0)
	# 	end, hist = sim.n_simulation_steps(10,start, make_history_array=True)

	# 	# print "start: ", start
	# 	# for i in hist[1:]:
	# 	# 	print " - ", i
	# 	# print "end: ", end

	# 	sim.plot(hist)

	# for y in xrange(40):
	# 	print "---------------------------------"
	# 	x_i = np.random.random()*policy.d0
	# 	y_i = np.random.random()*policy.d1
	# 	t_i = np.random.random()*np.pi * 2.0
	# 	start = ExtendedState(x_i,y_i,t_i, 1)
	# 	end, state_hist, action_hist = sim.n_simulation_steps(20,start, make_history_array=True)

	# 	# print "start: ", start
	# 	# for i in hist[1:]:
	# 	# 	print " - ", i
	# 	# print "end: ", end

	# 	sim.plot(state_hist, action_hist)

	# for i in xrange(20):
	# 	# print policy.at(State(0,i,0))

	# 	start = State(0,i,0)
	# 	print sim.simulate(start)
	print "DONE"


# # angles = [-0.3,0.0,0.3]
# angles = [0.35, 0.17]
# arrows = []

# for a in xrange(len(angles)):
# 	for i in np.linspace(0.0,5.0,num=10):
# 		ss = State(0,0,0)
# 		arrows.append(make_arrow(motion_model(Action(angles[a], 1.0), ss, step_size=i)))
# 		arrows.append(make_arrow(motion_model(Action(angles[a], -1.0), ss, step_size=i)))

# soa = np.array(arrows)
# # soa = np.array([make_arrow(State(1,1,0)),make_arrow(State(1,1,np.pi/2.0)),make_arrow(State(1,1,np.pi)),make_arrow(State(1,1,3.0*np.pi/2.0)),make_arrow(State(1,1,0))])
# X, Y, U, V = zip(*soa)
# plt.figure()
# ax = plt.gca()
# ax.quiver(X, Y, U, V, angles='xy', scale_units='xy', scale=1)
# ax.set_xlim([-3, 3])
# ax.set_ylim([-3, 3])
# # ax.set_xticks(np.arange(-3, 3, 0.1))
# # ax.set_yticks(np.arange(-3, 3., 0.1))
# plt.grid()
# plt.draw()
# plt.show()