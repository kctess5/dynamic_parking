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

#ifndef PARKING_H
#define PARKING_H

#include "vendor/lodepng/lodepng.h"

#include <stdio.h>      /* printf */
#include <cstdlib>
#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <algorithm>    // std::min
#include <time.h>
#include <chrono>
#include <set>
#include <iomanip>      // std::setw
#include <unistd.h>
#include <stdexcept>
#include <sstream>
// #define NDEBUG
#include <cassert>
#include <tuple>
#include <limits>

#define _EPSILON 0.0000001
#define M_2PI 6.2831853071795864769252867665590

#define _INTERP_AVOID_MEMREADS 1

#define RESCALE(ind, disc, min_val, max_val) ((max_val - min_val) * (float)ind / (float)(disc) + min_val)
#define UNSCALE(val, disc, min_val, max_val) ((disc) * (val - min_val) /(max_val - min_val))

// No inline
#if _NO_INLINE == 1
#define ANIL __attribute__ ((noinline))
#else
#define ANIL 
#endif

namespace dp {
	template <typename T>
	struct Array3D
	{
		T*** ptr;
		T* pool;

		int d0, d1, d2;

		Array3D(int x, int y, int z) : d0(x), d1(y), d2(z) {
			pool = new T[d0 * d1 * d2];
			ptr = new T**[d2];
			for (int i = 0; i < d2; ++i) {
				ptr[i] = new T*[d1];
				for (int j = 0; j < d1; ++j)
					ptr[i][j] = &pool[d1 * (j + d2 * i)];
			}
		}
		~Array3D(){
			for (int i = 0; i < d2; ++i)
				delete [] ptr[i];
			free(pool);
		}

		/* these methods index using pointers */
		T at2(int x, int y, int z) const{
			assert(x >= 0 && x < d0); assert(y >= 0 && y < d1); assert(z >= 0 && z < d2);
			return ptr[z][y][x]; 
		}
		T & at2(int x, int y, int z) {
			assert(x >= 0 && x < d0); assert(y >= 0 && y < d1); assert(z >= 0 && z < d2);
			return ptr[z][y][x]; 
		}

		/* these methods directly index the memory pool */
		T at(int x, int y, int z) const{
			assert(x >= 0 && x < d0); assert(y >= 0 && y < d1); assert(z >= 0 && z < d2);
			return pool[x + d1 * (y + d2 * z)]; 
		}
		T & at(int x, int y, int z) {
			assert(x >= 0 && x < d0); assert(y >= 0 && y < d1); assert(z >= 0 && z < d2);
			return pool[x + d1 * (y + d2 * z)]; 
		}

		std::vector<int> size() { return {d0, d1, d2}; }
		int size(int dim) {
			switch(dim) {
				case 0: return d0;
				case 1: return d1;
				case 2: return d2;
				default: assert(dim >= 0 && dim <= 2);
					return -1;
					break;
			}
		}

		T interp(float x, float y, float z) const{
			assert(x >= 0 && x <= d0-1); assert(y >= 0 && y <= d1-1);
			float x_l = fmodf(x, 1.f);
			float y_l = fmodf(y, 1.f);
			float z_l = fmodf(z, 1.f);

			float one_minus_x_l = 1.0 - x_l;
			float one_minus_y_l = 1.0 - y_l;
			float one_minus_z_l = 1.0 - z_l;

			int x_b = x;
			int y_b = y;

			// std::cout << x_l << std::endl;
			// std::cout << x_b << std::endl;
			// std::cout << one_minus_x_l << std::endl;

			// wrap the third dimension, which corresponds to theta
			int z_b = int(z) % d2;
			int z_plus_1 = (z_b + 1) % d2;

			#if _INTERP_AVOID_MEMREADS == 1
				float t1 = one_minus_x_l * one_minus_y_l * one_minus_z_l;
				float t2 = x_l * one_minus_y_l * one_minus_z_l;
				float t3 = one_minus_x_l * y_l * one_minus_z_l;
				float t4 = one_minus_x_l * one_minus_y_l * z_l;
				float t5 = x_l * one_minus_y_l * z_l;
				float t6 = one_minus_x_l * y_l * z_l;
				float t7 = x_l * y_l * one_minus_z_l;
				float t8 = x_l * y_l * z_l;
				T Vxyz = 0.0;
				if (t1 > _EPSILON || t1 < -_EPSILON) Vxyz += t1 * at(x_b,y_b,z_b);
				if (t2 > _EPSILON || t2 < -_EPSILON) Vxyz += t2 * at(x_b+1,y_b,z_b);
				if (t3 > _EPSILON || t3 < -_EPSILON) Vxyz += t3 * at(x_b,y_b+1,z_b);
				if (t4 > _EPSILON || t4 < -_EPSILON) Vxyz += t4 * at(x_b,y_b,z_plus_1);
				if (t5 > _EPSILON || t5 < -_EPSILON) Vxyz += t5 * at(x_b+1,y_b,z_plus_1);
				if (t6 > _EPSILON || t6 < -_EPSILON) Vxyz += t6 * at(x_b,y_b+1,z_plus_1);
				if (t7 > _EPSILON || t7 < -_EPSILON) Vxyz += t7 * at(x_b+1,y_b+1,z_b);
				if (t8 > _EPSILON || t8 < -_EPSILON) Vxyz += t8 * at(x_b+1,y_b+1,z_plus_1);
			#else
				T v000 = at(x_b,y_b,z_b);
				T v100 = at(x_b+1,y_b,z_b);
				T v010 = at(x_b,y_b+1,z_b);
				T v001 = at(x_b,y_b,z_plus_1);
				T v101 = at(x_b+1,y_b,z_plus_1);
				T v011 = at(x_b,y_b+1,z_plus_1);
				T v110 = at(x_b+1,y_b+1,z_b);
				T v111 = at(x_b+1,y_b+1,z_plus_1);
				T Vxyz = v000 * one_minus_x_l * one_minus_y_l * one_minus_z_l +
					v100 * x_l * one_minus_y_l * one_minus_z_l +
					v010 * one_minus_x_l * y_l * one_minus_z_l +
					v001 * one_minus_x_l * one_minus_y_l * z_l +
					v101 * x_l * one_minus_y_l * z_l +
					v011 * one_minus_x_l * y_l * z_l +
					v110 * x_l * y_l * one_minus_z_l +
					v111 * x_l * y_l * z_l;
			#endif
			return Vxyz;
		}

		void fill(T val) {
			for (int x = 0; x < d0; ++x)
				for (int y = 0; y < d1; ++y)
					for (int z = 0; z < d2; ++z)
						at(x,y,z) = val;
		}

		void replace(Array3D * other) {
			memcpy(pool, other->pool, d0*d1*d2*sizeof(pool[0]));
		}
	};

	struct OmniAction
	{
		float direction;
		bool stop;
		OmniAction() : direction(-1.0), stop(true) {}
		OmniAction(float val) : stop(false) { direction = fmodf(val, M_2PI);}
	};

	class OmniActionModel
	{
	public:
		OmniActionModel(float ss, int disc) : step_size(ss), discretization(disc) {};
		~OmniActionModel() {};

		std::vector<OmniAction> enumerate() {
			std::vector<OmniAction> actions;
			float coeff = M_2PI / float(discretization);
			for (int i = 0; i < discretization; ++i)
				actions.push_back(OmniAction(coeff * i));
			return actions;
		};

		void apply(OmniAction action, float x, float y, float z, float &ax, float &ay, float &az) {
			// az = z + action.direction;
			// ax = x + step_size * cosf(az);
			// ay = y + step_size * sinf(az);

			// az = 0.0;
			az = fmodf(z + action.direction, M_2PI);
			// ax = x + step_size * cos(action.direction);
			// ay = y + step_size * sin(action.direction);
			ax = x + step_size * cos(az);
			ay = y + step_size * sin(az);
		}

		float cost(OmniAction action) { return step_size; }
	
	private:
		float step_size;
		int discretization;
	};

	class ValueIterator
	{
	public:
		ValueIterator() {
			dims = {200,200,1};
			min_vals = {0.0,0.0, 0.0};
			max_vals = {1.0,1.0, M_2PI};

			J = new Array3D<float>(dims[0], dims[1], dims[2]);
			J->fill(std::numeric_limits<float>::max() / 2.0);
			Policy = new Array3D<OmniAction>(dims[0], dims[1], dims[2]);

			action_model = new OmniActionModel(0.01 * dims[0], 4);
			
			OmniAction null_action;
			Policy->fill(null_action);

			// mark the center region as the goal
			// for (int x = 50; x < 55; ++x)
			// 	for (int y = 50; y < 55; ++y)
			// 		for (int t = 0; t < dims[2]; ++t)
			// 			J->at(x,y,t) = 0.0;

			for (int x = 47; x < 54; ++x)
				for (int y = 47; y < 54; ++y)
					for (int t = 0; t < dims[2]; ++t)
						J->at(x,y,t) = 0.0;

			// auto start_time = std::chrono::high_resolution_clock::now();

			// std::cout << RESCALE(25, 50, 0.0, 1.0) << std::endl;
			// std::cout << UNSCALE((RESCALE(25, 50, 0.0, 1.0)), 50, 0.0, 1.0) << std::endl;

			// for (int i = 0; i < 50; ++i)
			// {
			// 	std::cout << UNSCALE((RESCALE(i, 50, 0.0, 1.0)), 50, 0.0, 1.0) << std::endl;
			// }
			// std::cout << RESCALE(25, 50, 0.0, 1.0) << std::endl;
			// std::cout << UNSCALE((RESCALE(25, 50, 0.0, 1.0)), 50, 0.0, 1.0) << std::endl;

			// policy iteration
			for (int i = 0; i < 150; ++i)
			{
				std::cout << "Step: " << i << std::endl;
				step();

				// save_slice(0,1,1,"test.png");
			}
			save_slice(0,1,0,"test.png");
			// save_slice(0,1,1,"test2.png");
			save_policy(0,1,0,"policy.png");
			// save_policy(0,1,1,"policy2.png");

			// auto end_time = std::chrono::high_resolution_clock::now();
			// std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
			// std::cout << "Done in: " << time_span.count() << std::endl;


			// save_slice(0,1,0,"test.png");
		};
		~ValueIterator() {};

		float cost_at(float x, float y, float z) {
			if (x < 0 || x > dims[0] - 1 || y < 0 || y > dims[1] - 1) {
				// bad++;
				return std::numeric_limits<float>::max() / 2.0;
			} else {
				// good++;
				return J->interp(x,y,z);
			}
		}

		// run a single step of value iteration
		void step() {
			std::vector<OmniAction> actions = action_model->enumerate();

			// preallocate space for the various variables for speed
			float sx, sy, sz;
			float ax, ay, az;
			float ux, uy, uz;
			int i, min_i;
			float min_cost;
			float alt_cost;

			Array3D<float>* J_prime = new Array3D<float>(dims[0], dims[1], dims[2]);
			J_prime->fill(std::numeric_limits<float>::max() / 2.0);

			// int x,y,z;
			// for (x = 0; x < dims[0]; ++x) {
			// 	for (y = 0; y < dims[1]; ++y) {
			// 		for (z = 0; z < dims[2]; ++z) {
			// 			sx = RESCALE(x, dims[0], min_vals[0], max_vals[0]);
			// 			sy = RESCALE(y, dims[1], min_vals[1], max_vals[1]);
			// 			sz = RESCALE(z, dims[2], min_vals[2], max_vals[2]);

			// 			// the baseline is the "do nothing" cost
			// 			min_i = -1;
			// 			min_cost = cost_at(x, y, z);
			// 			for (i = 0; i < actions.size(); ++i) {
			// 				action_model->apply(actions[i], sx, sy, sz, ax, ay, az);

			// 				ux = UNSCALE(ax, dims[0], min_vals[0], max_vals[0]);
			// 				uy = UNSCALE(ay, dims[1], min_vals[1], max_vals[1]);
			// 				uz = UNSCALE(az, dims[2], min_vals[2], max_vals[2]);

			// 				alt_cost = cost_at(ux, uy, uz);
			// 				// alt_cost = cost_at(ux, uy, 0.0);
			// 				if (alt_cost < min_cost && alt_cost < 100000000) {
			// 					min_cost = alt_cost;
			// 					min_i = i;
			// 				}
			// 			}
			// 			// if the do nothing option is optimal, do not modify cost, only policy
			// 			if (min_i == -1) {
			// 				Policy->at(x,y,z) = OmniAction();
			// 				J_prime->at(x,y,z) = min_cost;
			// 			} else {
			// 				// std::cout << "TEST" << std::endl;
			// 				// otherwise, modify both the policy and the minimum cost to reflect the found results
			// 				Policy->at(x,y,z) = actions[min_i];
			// 				J_prime->at(x,y,z) = min_cost + action_model->cost(actions[min_i]);
			// 				// std::cout << "TEST: " << min_cost << "  " << action_model->cost(actions[min_i]) << std::endl;
			// 			}
			// 		}
			// 	}
			// }

			int x,y,z;
			for (x = 0; x < dims[0]; ++x) {
				for (y = 0; y < dims[1]; ++y) {
					for (z = 0; z < dims[2]; ++z) {
						// the baseline is the "do nothing" cost
						min_i = -1;
						min_cost = cost_at(x, y, z);
						for (i = 0; i < actions.size(); ++i) {
							action_model->apply(actions[i], x, y, z, ax, ay, az);
							alt_cost = cost_at(ax, ay, az);
							if (alt_cost < min_cost && alt_cost < 100000000) {
								min_cost = alt_cost;
								min_i = i;
							}
						}
						// if the do nothing option is optimal, do not modify cost, only policy
						if (min_i == -1) {
							Policy->at(x,y,z) = OmniAction();
							J_prime->at(x,y,z) = min_cost;
						} else {
							// std::cout << "TEST" << std::endl;
							// otherwise, modify both the policy and the minimum cost to reflect the found results
							Policy->at(x,y,z) = actions[min_i];
							J_prime->at(x,y,z) = min_cost + action_model->cost(actions[min_i]);
							// std::cout << "TEST: " << min_cost << "  " << action_model->cost(actions[min_i]) << std::endl;
						}
					}
				}
			}

			J->replace(J_prime);
			// std::cout << "good: " << good << std::endl;
			// std::cout << "bad: " << bad << std::endl;
		}

		bool save_slice(int dim1, int dim2, int ind, std::string filename) {
			std::vector<unsigned char> png;
			lodepng::State state; //optionally customize this one
			int width = dims[0];
			int height = dims[1];
			char image[width * height * 4];

			float max_val = 0.0;

			for (int y = 0; y < height; ++y) {
				for (int x = 0; x < width; ++x) {
					unsigned idx = 4 * y * width + 4 * x;
					float val = J->at(x,y,ind);
					if (val < 100000 and val > max_val)
						max_val = val;
				}
			}

			for (int y = 0; y < height; ++y) {
				for (int x = 0; x < width; ++x) {
					unsigned idx = 4 * y * width + 4 * x;
					float val = J->at(x,y,ind);

					if (val > 100000) {
						val = 255;
					} else {
						val = 150.0 * val / max_val;
					}

					image[idx + 2] = (char)val;
					image[idx + 1] = (char)val;
					image[idx + 0] = (char)val;
					image[idx + 3] = (char)255;
				}
			}
			unsigned error = lodepng::encode(png, reinterpret_cast<const unsigned char*> (image), width, height, state);
			if(!error) lodepng::save_file(png, filename);
			//if there's an error, display it
			if(error) std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
			return error;
		}

		bool save_policy(int dim1, int dim2, int ind, std::string filename) {
			std::vector<unsigned char> png;
			lodepng::State state; //optionally customize this one
			int width = dims[0];
			int height = dims[1];
			char image[width * height * 4];

			float max_val = 0.0;

			for (int y = 0; y < height; ++y) {
				for (int x = 0; x < width; ++x) {
					unsigned idx = 4 * y * width + 4 * x;
					float val = Policy->at(x,y,ind).direction;
					if (val < 100000 and val > max_val)
						max_val = val;
				}
			}

			for (int y = 0; y < height; ++y) {
				for (int x = 0; x < width; ++x) {
					unsigned idx = 4 * y * width + 4 * x;
					float val = Policy->at(x,y,ind).direction;

					if (val > 100000) {
						val = 255;
					} else if (val == -1.0) {
						val = 255;
					} else {
						val = 150.0 * val / max_val;
					}

					image[idx + 2] = (char)val;
					image[idx + 1] = (char)val;
					image[idx + 0] = (char)val;
					image[idx + 3] = (char)255;
				}
			}
			unsigned error = lodepng::encode(png, reinterpret_cast<const unsigned char*> (image), width, height, state);
			if(!error) lodepng::save_file(png, filename);
			//if there's an error, display it
			if(error) std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
			return error;
		}
	
	protected:

		int good = 0;
		int bad = 0;

		Array3D<float>* J;
		Array3D<OmniAction>* Policy;

		OmniActionModel * action_model;
		std::vector<int> dims;
		std::vector<float> min_vals;
		std::vector<float> max_vals;
	};
} 

#endif