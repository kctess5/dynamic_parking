/*
* Copyright 2017 Corey H. Walsh (corey.walsh11@gmail.com)

* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at

*     http://www.apache.org/licenses/LICENSE-2.0

* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/

#include <gflags/gflags.h>

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

#include "./includes/parking.h"

// DEFINE_string(method, "RayMarching",
// 	"Which range method to use, one of:\n  BresenhamsLine (or bl)\n  RayMarching (or rm)\n  CDDTCast (or cddt)\n  PrunedCDDTCast (or pcddt)\n  GiantLUTCast (or glt)\n");
// DEFINE_string(map_path, "BASEMENT_MAP", "Path to map image, relative to current directory");
// DEFINE_string(cddt_save_path, "", "Path to serialize CDDT data structure to.");
// DEFINE_string(log_path, "", "Where to store high fidelity logs, does not save if not specified.");
// DEFINE_string(which_benchmark, "", "Which benchmark to run, one of:\n  random\n  grid\n");
// DEFINE_string(query, "", "Query point x,y,theta to ray cast from. example: --query=0,0,3.14");
// DEFINE_string(trace_path, "", "Path to output trace map of memory access pattern. Works for Bresenham's Line or Ray Marching.");
// DEFINE_string(lut_slice_path, "", "Path to output a slice of the LUT.");
DEFINE_string(test_flag, "hello", "Which LUT slice to output");

using namespace dp;

template <typename T>
void omni_action() {
	std::vector<int> dims = {200,200,4};
	std::vector<float> min_vals = {0.0,0.0, 0.0};
	std::vector<float> max_vals = {1.0,1.0, M_2PI};

	OmniValueIterator<T> *VI 
	    = new OmniValueIterator<T>(dims, min_vals, max_vals);

	OmniActionModel<T> *ActionModel = new OmniActionModel<T>(2.0, 4);
	VI->set_action_model(ActionModel);

	// set the desired states to zero cost
  //   for (int x = dims[0]/2-2; x < dims[0]/2+3; ++x)
		// for (int y = dims[1]/2-2; y < dims[1]/2+3; ++y)
		// 	for (int t = 0; t < dims[2]; ++t)
		// 		VI->getJ()->at(x,y,t) = 0.0;

	int steps = 105;

	auto start_time = std::chrono::high_resolution_clock::now();
	// policy iteration
	for (int i = 0; i < steps; ++i)
	{
		std::cout << "Step: " << i << std::endl;
		VI->step();

		VI->save_slice(0, 300.0, "./sequence/" + padded(i,3) + ".png");
		VI->save_policy(0,"./policy_sequence/" + padded(i,3) + ".png");
	}

	// auto end_time = std::chrono::high_resolution_clock::now();
	// std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
	// std::cout << "Done in: " << time_span.count() << std::endl;
	// std::cout << "  - per step: " << time_span.count() / steps << std::endl;
	
}

template <typename T>
void nonholonomic_action() {
	std::vector<int> dims = {70,70,180};
	std::vector<float> min_vals = {0.0,0.0, 0.0};
	std::vector<float> max_vals = {1.0,1.0, M_2PI};

	NonHolonomicIterator<T> *VI 
	    = new NonHolonomicIterator<T>(dims, min_vals, max_vals);

	NonHolonomicModel<T> *ActionModel = new NonHolonomicModel<T>(3.0, 4, 0.35, 1.0);

	// std::vector<NonHolonomicAction<T> > actions = ActionModel->enumerate();
	// for (int i = 0; i < actions.size(); ++i)
	// 	std::cout << actions[i].steering_angle << " " << actions[i].throttle << std::endl;

	VI->set_action_model(ActionModel);

	// set the desired states to zero cost
  //   for (int x = dims[0]/2-2; x < dims[0]/2+3; ++x)
		// for (int y = dims[1]/2-2; y < dims[1]/2+3; ++y)
		// 	for (int t = 0; t < dims[2]; ++t)
		// 		VI->getJ()->at(x,y,t) = 0.0;

	// float max_theta_err = 0.15;
	// // parallel park
	// for (int x = dims[0]/2-8; x < dims[0]/2+8; ++x)
	// 	for (int y = dims[1]/2-2; y < dims[1]/2+3; ++y)
	// 		for (int t = 0; t < dims[2]; ++t) {
	// 			float theta = M_2PI * (float)t / (float)dims[2];
	// 			if (theta <= max_theta_err || M_2PI - theta <= max_theta_err)
	// 				VI->getJ()->at(x,y,t) = 0.0;
				
	// 		}
				

	int steps = 26;
	int slice = 0;

	VI->save_slice(slice, 30.0, "./sequence/" + padded(0,3) + ".png");

	auto start_time = std::chrono::high_resolution_clock::now();
	// policy iteration
	for (int i = 0; i < steps; ++i)
	{
		std::cout << "Step: " << i << std::endl;
		VI->step();

		VI->save_slice(slice, 30.0, "./sequence/" + padded(i+1,3) + ".png");
		VI->save_policy(slice,"./policy_sequence/" + padded(i+1,3) + ".png");
	}

	// VI->getJ()->printSlice(5,3);

	for (int i = 0; i < dims[2]; ++i)
	{
		VI->save_slice(i, 50.0, "./slices/" + padded(i,3) + ".png");
		VI->save_policy(i,"./policy_slices/" + padded(i,3) + ".png");
	}

	// serialize policy and cost function to the file system for analysis elsewhere
	VI->serialize_cost("./serialized/cost.object");
	VI->serialize_policy("./serialized/policy.object");

	auto end_time = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
	std::cout << "Done in: " << time_span.count() << std::endl;
	std::cout << "  - per step: " << time_span.count() / steps << std::endl;
}

template <typename T>
void nonholonomic_action2() {
	// nonholonomic motion model with extended state space
	std::vector<int> dims = {210,210,60};
	std::vector<float> min_vals = {0.0,0.0, 0.0};
	std::vector<float> max_vals = {1.0,1.0, M_2PI};

	NonHolonomicIteratorLargeSpace<T> *VI 
	    = new NonHolonomicIteratorLargeSpace<T>(dims, min_vals, max_vals);

	NonHolonomicModel<T> *ActionModel = new NonHolonomicModel<T>(2.0, 8, 0.32, 10.4);
    // NonHolonomicModel<T> *ActionModel = new NonHolonomicModel<T>(2.0, 8, 0.32, 4.0);
	// this part is necessary for the fast apply statement
	ActionModel->memoize(dims[2]);
	// VI->load_map("../70x70_map.png");
	
	// goal info
	VI->load_map("../maps/210x210_kitchen2.png");
	
	// 3d occupancy obstacle info
	VI->load_3D_omap("../maps/dilated_kitchen/dilated_kitchen_summary.yaml");
	// return;

	// VI->load_map("../kitchen_2.png");
	// VI->load_map("../160x160_parallel2.png");

	// return;
// 
	// std::vector<NonHolonomicAction<T> > actions = ActionModel->enumerate();
	// for (int i = 0; i < actions.size(); ++i)
	// 	std::cout << actions[i].steering_angle << " " << actions[i].throttle << std::endl;

	VI->set_action_model(ActionModel);

	// std::cout << "blah: " << VI->is_goal(40, 35, 39) << std::endl;

	// std::cout << "blah1: " << VI->within_err(41, M_2PI / 2.0, 0.15) << std::endl;

	// std::cout << "blah1: " << VI->within_err(30, M_2PI / 2.0, 0.15) << std::endl;
	// std::cout << "blah2: " << VI->within_err(39, M_2PI / 2.0, 0.15) << std::endl;
	// std::cout << "blah3: " << VI->within_err(40, M_2PI / 2.0, 0.15) << std::endl;
	// std::cout << "blah4: " << VI->within_err(0, M_2PI, 0.15) << std::endl;
	// std::cout << "blah5: " << VI->within_err(1, M_2PI, 0.15) << std::endl;
	// std::cout << "blah6: " << VI->within_err(79, M_2PI, 0.15) << std::endl;
	// return;

	int steps = -1;
	// int steps = 40;
	int slice = 0;
	float max_cost = 200.0;

	VI->save_slice(slice, max_cost, "./sequence/" + padded(0,4) + ".png");

	auto start_time = std::chrono::high_resolution_clock::now();

	if (steps < 0) {
		int i = 0;
		bool useful_update = true;
		while (useful_update) {
			std::cout << "Step: " << i << std::endl;
			VI->step();
			VI->save_slice(slice, max_cost, "./sequence/" + padded(i+1,4) + ".png");
			VI->save_policy(slice,"./policy_sequence/" + padded(i+1,4) + ".png");
			i++;

			useful_update = VI->last_step_useful();
		}
	} else {
		// policy iteration
		for (int i = 0; i < steps; ++i)
		{
			std::cout << "Step: " << i << std::endl;
			VI->step();

			VI->save_slice(slice, max_cost, "./sequence/" + padded(i+1,4) + ".png");
			VI->save_policy(slice,"./policy_sequence/" + padded(i+1,4) + ".png");
		}
	}
	

	// VI->getJ()->printSlice(5,3);

	for (int i = 0; i < dims[2]; ++i)
	{
		VI->save_slice(i, max_cost, "./slices/" + padded(i,3) + ".png");
		VI->save_policy(i,"./policy_slices/" + padded(i,3) + ".png");
	}

	VI->save_full_policy("./serialized/parallel_park");

	// serialize policy and cost function to the file system for analysis elsewhere
	VI->serialize_cost("./serialized/cost.object");
	// VI->serialize_policy("./serialized/policy.object");
	// VI->serialize_policy2("./serialized/policy_reverse.object");

	auto end_time = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
	std::cout << "Done in: " << time_span.count() << std::endl;
	std::cout << "  - per step: " << time_span.count() / steps << std::endl;
}

int main(int argc, char *argv[])
{
	// set usage message
	std::string usage("This library provides fast 2D ray casting on occupancy grid maps.  Sample usage:\n\n");
	usage += "   ";
	usage += argv[0];
	google::SetUsageMessage(usage);
	gflags::ParseCommandLineFlags(&argc, &argv, true);

	std::cout << "start." << std::endl;

	// ValueIterator<float> vi = ValueIterator<float>();

	// omni_action<float>();
	// nonholonomic_action<float>();
	nonholonomic_action2<float>();

	std::cout << FLAGS_test_flag << std::endl;

	std::cout << "done." << std::endl;
	return 0;
}