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

void omni_action() {
	std::vector<int> dims = {120,120,20};
	std::vector<float> min_vals = {0.0,0.0, 0.0};
	std::vector<float> max_vals = {1.0,1.0, M_2PI};

	// initialize a value iterator
	// AbstractValueIterator<float, OmniAction<float>, OmniActionModel<float> > *VI 
	//     = new AbstractValueIterator<float, OmniAction<float>, OmniActionModel<float> >(dims, min_vals, max_vals);
	OmniValueIterator<float> *VI 
	    = new OmniValueIterator<float>(dims, min_vals, max_vals);

	OmniActionModel<float> *ActionModel = new OmniActionModel<float>(2.0, 3);
	VI->set_action_model(ActionModel);

	// set the desired states to zero cost
    for (int x = dims[0]/2-2; x < dims[0]/2+3; ++x)
		for (int y = dims[1]/2-2; y < dims[1]/2+3; ++y)
			for (int t = 0; t < dims[2]; ++t)
				VI->getJ()->at(x,y,t) = 0.0;

	int steps = 100;

	auto start_time = std::chrono::high_resolution_clock::now();
	// policy iteration
	for (int i = 0; i < steps; ++i)
	{
		std::cout << "Step: " << i << std::endl;
		VI->step();

		VI->save_slice(0, 120.0, "./sequence/" + padded(i,3) + ".png");
		VI->save_policy(0,"./policy_sequence/" + padded(i,3) + ".png");
	}

	// auto end_time = std::chrono::high_resolution_clock::now();
	// std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
	// std::cout << "Done in: " << time_span.count() << std::endl;
	// std::cout << "  - per step: " << time_span.count() / steps << std::endl;
	
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

	omni_action();

	std::cout << FLAGS_test_flag << std::endl;

	std::cout << "done." << std::endl;
	return 0;
}