import numpy as np
import matplotlib.pyplot as plt
import ujson, time
import csv
from collections import namedtuple
from scipy.interpolate import RegularGridInterpolator

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
 			# print "arc radius", R
 			dt = sign * action.throttle * step / R
			dy = sign * (1.0-np.cos(dt))*R
			dx = sign * np.sin(dt)*R

			# print dt, dx, dy

 		c, s = np.cos(state.theta), np.sin(state.theta)
		dx_global = c*dx - s*dy
		dy_global = s*dx + c*dy

		return State(dx_global + state.x, dy_global + state.y, (dt + state.theta)%(2.0*np.pi))

class AckermannModel(MotionModel):
 	""" Standard kinematic bicycle motion model """
 	def __init__(self, world_scale=1.0):
 		self._EPSILON = 0.0000001
 		self._SCALE = world_scale

 		self.max_steer = 0.335

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

# class ObstacleSpace(object):
# 	""" A 3D configuration space for obstacles """
# 	def __init__(self, arg):
# 		super(ObstacleSpace, self).__init__()
# 		self.arg = arg
		

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

class ContinuousPolicy(object):
	""" Given a continuous state, determines the optimal policy by backing up that state """
	def __init__(self, cost, model, discretization=3.0, step=1.0):
		self.model = model
		self.cost = cost
		self.deflection_penalty = 0.0
		self.reversal_penalty = 1.5
		# self.reversal_penalty = 0.0
		self.disc = discretization
		self.step = step
		self.actions = self.enumerate()
		self.discount = 1.0
		self._EPSILON = 0.00000001

	def enumerate(self):
		actions = []
		angles = np.linspace(-self.model.max_steer, self.model.max_steer, num=self.disc*2+1)
		actions += map(lambda a: Action(a, 1.0), angles)
		actions += map(lambda a: Action(a, -1.0), angles)
		return actions

	def backup(self, state):
		# print "BACKING UP STATE", state
		# check if this is a goal state, if so, stop moving
		cost_at_state = self.cost.interp(state)
		optimal_action = None
		cost = self.cost.max_cost

		world_state = ExtendedState(state.x, state.y, rescale_t(state.theta, self.cost.d2), state.throttle)
		# print "current throttle: ", world_state.throttle

		if cost_at_state == 0.0:
			cost = 0.0
			optimal_action = Action(0,0)
		else:
			min_action = Action(0,0)
			min_cost = self.cost.max_cost*10.0

			for action in self.actions:

				alt_state = self.model.apply(action, world_state, self.step)
				alt_state_extended = ExtendedState(alt_state.x, alt_state.y, unscale_t(alt_state.theta, self.cost.d2), action.throttle)
				

				cost_at_alt_state = self.discount * self.cost.interp(alt_state_extended)# + self.step

				# print
				# print "action:", action
				# print "alt state:", alt_state
				# print "alt cost:", cost_at_alt_state
				# print "  - ", alt_state_extended.throttle, alt_state_extended.throttle != state.throttle, cost_at_alt_state
				if abs(state.throttle) > 0 and abs(action.throttle) > 0 and alt_state_extended.throttle != state.throttle:
					# print "     rp"
					cost_at_alt_state += self.reversal_penalty

				

				if cost_at_alt_state < min_cost-self._EPSILON:
					# print "       ***"
					min_cost = cost_at_alt_state
					min_action = action

			cost = min_cost
			optimal_action = min_action

		return cost, optimal_action

	def at(self, state):
		return self.backup(state)[1]
		
class Cost(object):
	"""docstring for NHPolicy"""
	def __init__(self, path):
		self.max_cost = 100000.0

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

		self.backwards_cost = np.zeros(self.d0 * self.d1 * self.d2, dtype=np.float32)
		self.forwards_cost  = np.zeros(self.d0 * self.d1 * self.d2, dtype=np.float32)
		
		print "..loading forwards cost"
		f_file = open(path + "/cost_forwards.object", 'r')
	
		i = 0
		for row in f_file:
			self.forwards_cost[i] = float(row)
			i += 1

		print "..loading backwards cost"
		b_file = open(path + "/cost_forwards.object", 'r')
		i = 0
		for row in b_file:
			self.backwards_cost[i] = float(row)
			i += 1

		f_file.close()
		b_file.close()
		summary_file.close()

		# convert cost function into 3d arrays for easier interpolation
		self.forwards_cost_3d = np.reshape(self.forwards_cost, (self.d0,self.d1,self.d2), order="F")
		self.backwards_cost_3d = np.reshape(self.backwards_cost, (self.d0,self.d1,self.d2), order="F")

		self.forwards_interp = RegularGridInterpolator((
			np.linspace(0, self.d0, self.d0, endpoint=False), 
			np.linspace(0, self.d1, self.d1, endpoint=False), 
			np.linspace(0, self.d2, self.d2, endpoint=False)), self.forwards_cost_3d)

		self.backwards_interp = RegularGridInterpolator((
			np.linspace(0, self.d0, self.d0, endpoint=False), 
			np.linspace(0, self.d1, self.d1, endpoint=False), 
			np.linspace(0, self.d2, self.d2, endpoint=False)), self.backwards_cost_3d)

	def at(self, state, interpolate=False):
		x = int(np.round(state.x))
		y = int(np.round(state.y))
		z = int(np.round(state.theta) % self.d2)
		if state.throttle > 0:
			cost = self.forwards_cost[x+self.d0*(y+self.d1*z)]
		else:
			cost = self.backwards_cost[x+self.d0*(y+self.d1*z)]
		return cost

	def interp(self, state, interpolate=False):
		x = state.x
		y = state.y

		if x < 0 or x > self.d0 - 1 or y < 0 or y > self.d1 - 1:
			return self.max_cost * 2.0

		z = state.theta % self.d2

		# need to wrap the theta axis manually to implement periodic boundary conditions
		if z > self.d2 - 1:
			pts = np.array([[x,y,self.d2 - 1], [x,y,0]])
			if state.throttle > 0:
				cost = self.forwards_interp(pts)
			else:
				cost = self.backwards_interp(pts)
			t = z % 1.0
			return t * cost[1] + (1.0 - t) * cost[0]

		pts = np.array([[x,y,z]])

		if state.throttle > 0:
			# print "Forwards:", state.throttle
			cost = self.forwards_interp(pts)[0]
		else:
			# print "Backwards:", state.throttle
			cost = self.backwards_interp(pts)[0]
		return cost

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
	def __init__(self, Policy, Model, Cost, step):
		self.Policy = Policy
		self.Cost = Cost
		self.Model = Model
		self.step = step
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

			if action.throttle == 0:
				break
		
		if make_history_array:
			return state, state_hist, action_hist
		else:
			return state

	def get_optimal_policy(self, state):
		# print "GETTING POLICY:", state
		return self.Policy.at(state)

	# run a single simulation step based on the transition function and optimal policy
	def simulate(self, start_state, interpolate=False, give_action=False):
		# print start_state
		# if start_state.throttle == 0:
		# 	# print "STOP"
		# 	optimal_control = self.get_optimal_policy(start_state)
		# 	new_state = start_state
		# else:
		# 	optimal_control = self.get_optimal_policy(start_state)
		# 	new_state = self.transition(start_state, optimal_control)

		optimal_control = self.get_optimal_policy(start_state)
		new_state = self.transition(start_state, optimal_control)

		if self.last_throttle == 0 and not optimal_control.throttle == self.last_throttle:
			print "Accelerate: 0 -->",  optimal_control.throttle
		elif not optimal_control.throttle == self.last_throttle:
			print "Reversal: ", self.last_throttle, "-->", optimal_control.throttle
		
		self.last_throttle = optimal_control.throttle

		if give_action:
			return new_state, optimal_control
		else:
			return new_state

	def transition(self, state, action):
		# print
		# print state.theta
		unscaled_state = State(state.x,state.y, rescale_t(state.theta, self.Cost.d2))
		# print unscaled_state.theta
		ns = self.Model.apply(action, unscaled_state, self.step)

		noise = [0.0,0.0,0.0]

		if True:
			noise = (np.random.rand(3) - np.array([0.5,0.5,0.5])) * np.array([0.7,0.5,0.1])

		# print ns.theta
		scaled_new_state = ExtendedState(ns.x+noise[0],ns.y+noise[1], unscale_t(ns.theta+noise[2], self.Cost.d2), action.throttle)
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
		print self.down_spot, "-->", up_spot, angle, unscale_t(angle, self.Cost.d2)
		self.do_the_thing(State(start_x, start_y, unscale_t(angle, self.Cost.d2)))
		
	def do_the_thing(self, state):
		self.last_throttle = 0
		start = ExtendedState(state.x,state.y,state.theta, 0)
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
			s = State(state.x, state.y, rescale_t(state.theta, self.Cost.d2))
			arrows.append(make_arrow(s))
			
			if i < len(action_hist):
				lx, ly = self.Model.make_curves(s, action_hist[i], self.step, 10)
				line_xs += lx
				line_ys += ly
		
		soa = np.array(arrows)
		# soa = np.array([make_arrow(State(1,1,0)),make_arrow(State(1,1,np.pi/2.0)),make_arrow(State(1,1,np.pi)),make_arrow(State(1,1,3.0*np.pi/2.0)),make_arrow(State(1,1,0))])
		X, Y, U, V = zip(*soa)
		plt.plot(line_xs, line_ys)
		self.ax.quiver(X, Y, U, V, angles='xy', scale_units='xy', scale=1)
		self.ax.set_xlim([0, self.Cost.d0])
		self.ax.set_ylim([0, self.Cost.d1])
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
			s = State(state.x, state.y, rescale_t(state.theta, self.Cost.d2))
			arrows.append(make_arrow(s))
			
			if i < len(action_hist):
				lx, ly = self.Model.make_curves(s, action_hist[i], self.step, 10)
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
		ax.set_xlim([0, self.Cost.d0])
		ax.set_ylim([0, self.Cost.d1])
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
		self.ax.set_xlim([0, self.Cost.d0])
		self.ax.set_ylim([0, self.Cost.d1])
		plt.show()

if __name__ == '__main__':
	step = 3.0
	world_scale = 41.67
	action_disc = 8

	t = time.time()
	# policy = NHPolicy2("./build/serialized/parallel_park")
	# policy = NHPolicy2("./build/serialized/parallel_park_basement_left")
	# precomputed_policy = NHPolicy2("./build/serialized/parallel_park_modified2")
	cost   = Cost("./build/serialized/parallel_park_modified2")
	model  = EmpericalModel(world_scale=world_scale)
	# model  = AckermannModel(world_scale=world_scale)
	policy = ContinuousPolicy(cost, model, discretization=action_disc, step=step)

	print "Loaded in:", time.time() - t

	# exit() 

	
	s = ExtendedState(0,0,0,0)
	# print policy.at(s)
	sim = SimulatorExtendedState(policy, model, cost, step)

	# sim.load_image("70x70_map.png")
	# sim.load_image("./maps/210x210_kitchen.png")
	# sim.load_image("./maps/small_kitchen.png")
	sim.load_image("./maps/kitchen_box3.png")
	# sim.load_image("./maps/parking_goal_map_left.png")
	# sim.load_image("160x160_parallel2.png")
	sim.interactive_plot()
	print "DONE"

# if __name__ == '__main__':
# 	world_scale = 1.0
# 	max_angle = 0.335
# 	disc = 19
# 	steering_angles = np.linspace(max_angle, -max_angle, num=disc)
# 	# steering_angles = np.linspace(max_angle, 0.0, num=disc)
# 	speeds = [-1, 1]
# 	actions = []

# 	emp = EmpericalModel(world_scale=world_scale)
# 	ack = AckermannModel(world_scale=world_scale)

# 	# sa = np.linspace(0.05, 0.15, num=18)
# 	# for a in sa:
# 	# 	print a, "-->", model.steering_arc_radius(a)

# 	for a in steering_angles:
# 		actions.append(Action(a, speeds[0]))
# 		actions.append(Action(a, speeds[1]))
# 		# actions.append(Action(a, -1.0))

# 	state1 = State(0,0,0)
# 	# state2 = State(1.25,0,np.pi/2.0)

# 	state2 = State(0,0,0)

# 	step = 1.5

# 	print list(steering_angles)

# 	for a in actions:
# 		print
# 		# print a

# 		lx, ly = ack.make_curves(state2, a, step, 20, scale_step=True)
# 		plt.plot(lx, ly, color="blue")


# 		# newthing = emp.apply(a, state2, 2.0)
# 		# print a.steering_angle, newthing.x, newthing.y, newthing.theta

# 		# lx, ly = emp.make_curves(state2, a, step, 20, scale_step=True)
# 		# plt.plot(lx, ly, color="black")

# 	mv = 2.0 * world_scale
# 	plt.xlim([-1.5*mv, 1.5*mv])
# 	plt.ylim([-mv, mv])
# 	plt.show()
