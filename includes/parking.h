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

#define MAX_COST 1000000

#define DISCOUNT 0.999

#define CHANGE_PENALTY 0.05

#define _INTERP_AVOID_MEMREADS 1
#define _REPLACE_J 1

#define RESCALE(ind, disc, min_val, max_val) ((max_val - min_val) * (float)ind / (float)(disc) + min_val)
#define UNSCALE(val, disc, min_val, max_val) ((disc) * (val - min_val) /(max_val - min_val))

#define RESCALE_T(ind, disc) (M_2PI * (state_T)(ind) / (state_T)(disc))
#define UNSCALE_T(val, disc) ((disc) * (val) / M_2PI)


// No inline
#if _NO_INLINE == 1
#define ANIL __attribute__ ((noinline))
#else
#define ANIL 
#endif


std::string padded(int val, int num_zeros) {
  std::ostringstream ss;
  ss << std::setw(num_zeros) << std::setfill('0') << val;
  return ss.str();
}

namespace dp {
	template <typename T>
	struct Array3D
	{
		// T*** ptr;
		T* pool;

		int d0, d1, d2;

		Array3D(int x, int y, int z) : d0(x), d1(y), d2(z) {
			pool = new T[d0 * d1 * d2];
		}
		~Array3D(){
			free(pool);
		}

		/* these methods directly index the memory pool */
		T at(int x, int y, int z) const{
			assert(x >= 0 && x < d0); assert(y >= 0 && y < d1); assert(z >= 0 && z < d2);
			return pool[x+d0*(y+d1*z)]; 
		}
		T & at(int x, int y, int z) {
			assert(x >= 0 && x < d0); assert(y >= 0 && y < d1); assert(z >= 0 && z < d2);
			return pool[x+d0*(y+d1*z)]; 
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
			float x_l = fmod(x, 1.f);
			float y_l = fmod(y, 1.f);
			float z_l = fmod(z, 1.f);

			float one_minus_x_l = 1.0 - x_l;
			float one_minus_y_l = 1.0 - y_l;
			float one_minus_z_l = 1.0 - z_l;

			int x_b = x;
			int y_b = y;

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
				T v010 = 0;
				T v011 = 0;
				T v100 = 0;
				T v101 = 0;
				T v110 = 0;
				T v111 = 0;
				if (y_b+1 < d1) {
					v010 = at(x_b,y_b+1,z_b);
					v011 = at(x_b,y_b+1,z_plus_1);
				}
				if (x_b+1 < d0) {
					v100 = at(x_b+1,y_b,z_b);
					v101 = at(x_b+1,y_b,z_plus_1);
				}
				if (y_b+1 < d1 && x_b+1 < d0) {
					v110 = at(x_b+1,y_b+1,z_b);
					v111 = at(x_b+1,y_b+1,z_plus_1);
				}
				T v000 = at(x_b,y_b,z_b);
				T v001 = at(x_b,y_b,z_plus_1);

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
			for (int i = 0; i < d0*d1*d2; ++i)
				pool[i] = val;
		}

		void replace(Array3D * other) {
			memcpy(pool, other->pool, d0*d1*d2*sizeof(pool[0]));
		}

		void printSlice(int which, int places) {
			float factor = pow(10.0, places);
			for (int y = 0; y < d1; ++y){
				for (int x = 0; x < d0; ++x)
					std::cout << round(at(x,y,which) * factor) / factor << "  ";
				std::cout << std::endl;
			}
		}
	};

	template <typename T>
	struct OmniAction
	{
		T direction;
		bool stop;
		OmniAction() : direction(-1.0), stop(true) {}
		OmniAction(T val) : stop(false) { direction = fmod(val, M_2PI);}
	};

	template <typename T>
	class OmniActionModel
	{
	public:
		OmniActionModel(T ss, int disc) : step_size(ss), discretization(disc) {};
		~OmniActionModel() {};

		std::vector<OmniAction<T> > enumerate() {
			std::vector<OmniAction<T> > actions;
			T coeff = M_2PI / (T)discretization;
			for (int i = 0; i < discretization; ++i)
				actions.push_back(OmniAction<T>(coeff * i));
			// stop action
			actions.push_back(OmniAction<T>());
			return actions;
		};

		void apply(OmniAction<T> action, T x, T y, T z, T &ax, T &ay, T &az) {

			if (action.direction < 0) {
				ax = x; ay = y; az = z;
			} else {
				// az = action.direction;
				az = fmod(z + action.direction, M_2PI);
				ax = x + step_size * cos(az);
				ay = y + step_size * sin(az);
			}
			// az = z + action.direction;
			// ax = x + step_size * cosf(az);
			// ay = y + step_size * sinf(az);

			// az = 0.0;
			// ax = x + step_size * cos(action.direction);
			// ay = y + step_size * sin(action.direction);
			
			// az = fmodf(z + action.direction, M_2PI);
			// az = action.direction;
			// az = action.direction + z;
		}

		T cost(OmniAction<T> action) { 
			if (action.direction < 0) return 0.0;
			if (action.direction == 0) return step_size;
			return step_size + CHANGE_PENALTY; 
		}
	
	private:
		T step_size;
		int discretization;
	};

	template <typename state_T>
	class ValueIterator
	{
	public:
		ValueIterator() {
			dims = {200,200,10};
			min_vals = {0.0,0.0, 0.0};
			max_vals = {1.0,1.0, M_2PI};

			J = new Array3D<state_T>(dims[0], dims[1], dims[2]);
			J->fill(MAX_COST);
			Policy = new Array3D<OmniAction<state_T> >(dims[0], dims[1], dims[2]);

			#if _REPLACE_J == 1
			J_prime = new Array3D<state_T>(dims[0], dims[1], dims[2]);
			#endif

			action_model = new OmniActionModel<state_T>(2.0, 8);
			
			OmniAction<state_T> null_action;
			Policy->fill(null_action);

			// for (int x = 47; x < 54; ++x)
			// 	for (int y = 47; y < 54; ++y)
			// 		for (int t = 0; t < dims[2]; ++t)
			// 			J->at(x,y,t) = 0.0;

			for (int x = dims[0]/2-2; x < dims[0]/2+3; ++x)
				for (int y = dims[1]/2-2; y < dims[1]/2+3; ++y)
					for (int t = 0; t < dims[2]; ++t)
						J->at(x,y,t) = 0.0;

			// for (int x = 97; x < 104; ++x)
			// 	for (int y = 97; y < 104; ++y)
			// 		for (int t = 0; t < dims[2]; ++t)
			// 			J->at(x,y,t) = 0.0;

			// for (int x = 6; x < 10; ++x)
			// 	for (int y = 6; y < 10; ++y)
			// 		for (int t = 0; t < dims[2]; ++t)
			// 			J->at(x,y,t) = 0.0;

			// for (int t = 0; t < dims[2]; ++t){
			// 	J->at(4,4,t) = 0.0;
			// 	J->at(4,5,t) = 0.0;
			// 	J->at(5,4,t) = 0.0;
			// 	J->at(5,5,t) = 0.0;
			// }

			int steps = 80;

			auto start_time = std::chrono::high_resolution_clock::now();
			// policy iteration
			for (int i = 0; i < steps; ++i)
			{
				std::cout << "Step: " << i << std::endl;
				step();

				save_slice(0,1,0,"./sequence/" + padded(i,3) + ".png");
				save_policy(0,1,0,"./policy_sequence/" + padded(i,3) + ".png");
			}

			auto end_time = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
			std::cout << "Done in: " << time_span.count() << std::endl;
			std::cout << "  - per step: " << time_span.count() / steps << std::endl;

			// std::cout << "slice 1" << std::endl;
			// J->printSlice(0, 3);

			// if (dims[2] > 1) {
			// 	std::cout << "slice 2" << std::endl;
			// 	J->printSlice(1, 3);
			// }

			// std::cout << "slice 1" << std::endl;
			// J->printSlice(0, 3);

			// std::cout << "slice 2" << std::endl;
			// J->printSlice(1, 3);

			// int places = 3;
			// float factor = pow(10.0, places);

			// for (int y = 0; y < dims[1]; ++y)
			// {
			// 	for (int x = 0; x < dims[0]; ++x)
			// 	{
			// 		std::cout << round(J->at(x,y,0) * factor) / factor << "  ";
			// 	}
			// 	std::cout << std::endl;
			// }
			save_slice(0,1,0,"test.png");
			save_policy(0,1,0,"policy.png");
			if (dims[2] > 1) {
				save_slice(0,1,1,"test2.png");
				save_policy(0,1,1,"policy2.png");
			}


			


			// save_slice(0,1,0,"test.png");
		};
		~ValueIterator() {};

		state_T cost_at(state_T x, state_T y, state_T z) {
			if (x < 0 || x > dims[0] - 1 || y < 0 || y > dims[1] - 1) {
				// out of bounds cost
				return MAX_COST * 2.0;
			} else {
				return J->interp(x,y,z);
			}
		}

		// run a single step of value iteration
		void step() {
			std::vector<OmniAction<state_T> > actions = action_model->enumerate();
			OmniAction<state_T> stop = OmniAction<state_T>();

			// preallocate space for the various variables for speed
			state_T sx, sy, sz;
			state_T ax, ay, az;
			state_T ux, uy, uz;
			int i, min_i;
			state_T min_cost;
			state_T alt_cost;
			J_prime->fill(MAX_COST);

			// bool debug = 1;
			int x,y,z;
			for (x = 0; x < dims[0]; ++x) {
				for (y = 0; y < dims[1]; ++y) {
					for (z = 0; z < dims[2]; ++z) {
						min_i = -1;
						min_cost = MAX_COST*10.0;
						// if (debug && z == 0) std::cout << "state: ("<<x<<", "<<y<<", "<<z<<"), cost: " << min_cost << std::endl;

						for (i = 0; i < actions.size(); ++i) {
							action_model->apply(actions[i], x, y, RESCALE_T(z, dims[2]), ax, ay, az);
							az = UNSCALE_T(az, dims[2]);
							
							alt_cost = cost_at(ax, ay, az) + action_model->cost(actions[i]);
							// if (debug && z == 0) std::cout << "     + "<<actions[i].direction<<" --> ("<<ax<< ", "<<ay<<", "<<az<<") @ "<<alt_cost;
							if (alt_cost < min_cost-_EPSILON) {
								min_cost = alt_cost;
								min_i = i;
								// if (debug && z == 0) std::cout << "***";
							}
							// if (debug && z == 0) std::cout <<std::endl;
						}
						
						// This should never happen - the algorithm must always make a choice
						assert(min_i >= 0);

						// modify the policy and the cost to reflect the iteration
						Policy->at(x,y,z) = actions[min_i];
						J_prime->at(x,y,z) = min_cost;
					}
				}
			}
			J->replace(J_prime);
		}

		bool save_slice(int dim1, int dim2, int ind, std::string filename) {
			std::vector<unsigned char> png;
			lodepng::State state; //optionally customize this one
			int width = dims[0];
			int height = dims[1];
			char image[width * height * 4];

			float max_val = 0.0;

			// for (int y = 0; y < height; ++y) {
			// 	for (int x = 0; x < width; ++x) {
			// 		unsigned idx = 4 * y * width + 4 * x;
			// 		float val = J->at(x,y,ind);
			// 		if (val < 100 and val > max_val)
			// 			max_val = val;
			// 	}
			// }

			max_val = 120.0;

			for (int y = 0; y < height; ++y) {
				for (int x = 0; x < width; ++x) {
					unsigned idx = 4 * y * width + 4 * x;
					float val = J->at(x,y,ind);

					if (val > 1000) {
						image[idx + 2] = (char)255;
						image[idx + 1] = (char)255;
						image[idx + 0] = (char)255;
						image[idx + 3] = (char)255;
					} else {
						val = 150.0 * val / max_val;
						image[idx + 2] = (char)val;
						image[idx + 1] = (char)val;
						image[idx + 0] = (char)val;
						image[idx + 3] = (char)255;
					}

					
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
					if (val < 100 and val > max_val)
						max_val = val;
				}
			}

			for (int y = 0; y < height; ++y) {
				for (int x = 0; x < width; ++x) {
					unsigned idx = 4 * y * width + 4 * x;
					float val = Policy->at(x,y,ind).direction;

					if (val == -1.0) {
						image[idx + 2] = (char)255;
						image[idx + 1] = (char)0;
						image[idx + 0] = (char)0;
						image[idx + 3] = (char)255;
					} else {
						if (val > 100) {
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

		Array3D<state_T>* J;
		Array3D<OmniAction<state_T> >* Policy;

		#if _REPLACE_J == 1
		Array3D<state_T>* J_prime;
		#endif

		OmniActionModel<state_T> * action_model;
		std::vector<int> dims;
		std::vector<float> min_vals;
		std::vector<float> max_vals;
	};
} 

#endif