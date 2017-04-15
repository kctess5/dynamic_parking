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

int main(int argc, char *argv[])
{
	// set usage message
	std::string usage("This library provides fast 2D ray casting on occupancy grid maps.  Sample usage:\n\n");
	usage += "   ";
	usage += argv[0];
	google::SetUsageMessage(usage);
	gflags::ParseCommandLineFlags(&argc, &argv, true);

	std::cout << "start." << std::endl;

	ValueIterator vi = ValueIterator();

	std::cout << FLAGS_test_flag << std::endl;

	std::cout << "done." << std::endl;
	return 0;
}