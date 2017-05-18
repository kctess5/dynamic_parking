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
#include "vendor/distance_transform.h"


// #include <stdlib.h>
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
#include <bitset>
#include <limits>

#define _EPSILON 0.0001
#define M_2PI 6.2831853071795864769252867665590

#define MAX_COST 100000
#define DISCOUNT 1.0
#define CHANGE_PENALTY 0.05
#define REVERSAL_PENALTY 10.0
#define DEFLECTION_PENALTY 0.0

#define _USE_INTERPOLATION 1
#define _INTERP_AVOID_MEMREADS 1
#define _USE_FAST_APPLY 1
#define _USE_HEURISTIC_PRUNING 1

#define RESCALE(ind, disc, min_val, max_val) ((max_val - min_val) * (float)ind / (float)(disc) + min_val)
#define UNSCALE(val, disc, min_val, max_val) ((disc) * (val - min_val) /(max_val - min_val))

#define RESCALE_T(ind, disc) (M_2PI * (state_T)(ind) / (state_T)(disc))
#define UNSCALE_T(val, disc) ((disc) * (val) / M_2PI)

// these defines are for yaml/JSON serialization
#define T1 "  "
#define T2 T1 T1
#define T3 T1 T1 T1
#define T4 T2 T2

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

// templated sign function
template <typename T> T sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

namespace dp {
	// You could also take an existing vector as a parameter.
	std::vector<std::string> split(std::string str, char delimiter) {
	  std::vector<std::string> internal;
	  std::stringstream ss(str); // Turn the string into a stream.
	  std::string tok;
	  
	  while(getline(ss, tok, delimiter)) {
	    internal.push_back(tok);
	  }
	  
	  return internal;
	}

	float rgb2gray(float r, float g, float b) {
		return 0.229 * r + 0.587 * g + 0.114 * b;
	}
	// occupancy grid map, useful for computing the heuristic pruning method
	struct OMap
	{
		bool has_error;
		unsigned width;  // x axis
		unsigned height; // y axis
		std::vector<std::vector<bool> > grid;
		std::string fn; // filename

		OMap(int w, int h) : width(w), height(h), fn(""), has_error(false) {
			for (int i = 0; i < w; ++i) {
				std::vector<bool> y_axis;
				for (int q = 0; q < h; ++q) y_axis.push_back(false);
				grid.push_back(y_axis);
			}
		}

		OMap(std::string filename) : OMap(filename, 128) {}
		OMap(std::string filename, float threshold) : fn(filename), has_error(false) {
			unsigned error;
			unsigned char* image;

			error = lodepng_decode32_file(&image, &width, &height, filename.c_str());
			if(error) {
				printf("ERROR %u: %s\n", error, lodepng_error_text(error));
				has_error = true;
				return;
			}

			for (int i = 0; i < width; ++i) {
				std::vector<bool> y_axis;
				for (int q = 0; q < height; ++q) y_axis.push_back(false);
				grid.push_back(y_axis);
			}

			for (int y = 0; y < height; ++y) {
				for (int x = 0; x < width; ++x) {
					unsigned idx = 4 * y * width + 4 * x;
					int r = image[idx + 2];
					int g = image[idx + 1];
					int b = image[idx + 0];
					int gray = (int) rgb2gray(r,g,b);
					if (gray < threshold) grid[x][y] = true;
				}
			}
		}

		bool get(int x, int y) { return grid[x][y]; }
		bool isOccupied(int x, int y) { return grid[x][y]; }

		void fill(bool val) {
			for (int x = 0; x < width; ++x)
				for (int y = 0; y < height; ++y)
					grid[x][y] = val;
		}

		bool save(std::string filename) {
			std::vector<unsigned char> png;
			lodepng::State state; //optionally customize this one
			// char image = new char[width * height * 4] = 0;
			char image[width * height * 4];

			for (int y = 0; y < height; ++y) {
				for (int x = 0; x < width; ++x) {
					unsigned idx = 4 * y * width + 4 * x;
					
					image[idx + 2] = (char)255;
					image[idx + 1] = (char)255;
					image[idx + 0] = (char)255;
					image[idx + 3] = (char)255;

					if (grid[x][y]) {
						image[idx + 0] = 0;
						image[idx + 1] = 0;
						image[idx + 2] = 0;
					}
				}
			}
			unsigned error = lodepng::encode(png, reinterpret_cast<const unsigned char*> (image), width, height, state);
			if(!error) lodepng::save_file(png, filename);
			//if there's an error, display it
			if(error) std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
			return error;
		}
	};

	struct DistanceTransform
	{
		unsigned width;
		unsigned height;
		std::vector<std::vector<float> > grid;

		float get(int x, int y) { return grid[x][y]; }

		DistanceTransform() : width(0), height(0) {}

		DistanceTransform(int w, int h) : width(w), height(h) {
			// allocate space in the vectors
			for (int i = 0; i < width; ++i) {
				std::vector<float> y_axis;
				for (int q = 0; q < height; ++q) y_axis.push_back(1.0);
				grid.push_back(y_axis);
			}
		}

		// computes the distance transform of a given OMap
		DistanceTransform(OMap *map) {
			width = map->width;
			height = map->height;

			std::vector<std::size_t> grid_size({width, height});
		    dt::MMArray<float, 2> f(grid_size.data());
		    dt::MMArray<std::size_t, 2> indices(grid_size.data());

		    for (std::size_t i = 0; i < width; ++i)
		        for (std::size_t j = 0; j < height; ++j)
		        	if (map->isOccupied(i,j)) f[i][j] = 0.0f;
		        	else f[i][j] = std::numeric_limits<float>::max();
		    
			dt::DistanceTransform::distanceTransformL2(f, f, indices, false);

			// allocate space in the vectors
			for (int i = 0; i < width; ++i) {
				std::vector<float> y_axis;
				for (int q = 0; q < height; ++q) y_axis.push_back(1.0);
				grid.push_back(y_axis);
			}

			// store to array
			for (int y = 0; y < height; y++) {
				for (int x = 0; x < width; x++) {
					grid[x][y] = f[x][y];
				}
			}
		}

		void fill(float val) {
			for (int x = 0; x < width; ++x)
				for (int y = 0; y < height; ++y)
					grid[x][y] = val;
		}

		void init(OMap *map) {
			assert(width == map->width);
			assert(height == map->height);

			std::vector<std::size_t> grid_size({width, height});
		   dt::MMArray<float, 2> f(grid_size.data());
		   dt::MMArray<std::size_t, 2> indices(grid_size.data());

		   for (std::size_t i = 0; i < width; ++i)
		   	for (std::size_t j = 0; j < height; ++j)
		  			if (map->isOccupied(i,j)) f[i][j] = 0.0f;
		        	else f[i][j] = std::numeric_limits<float>::max();
		    
			dt::DistanceTransform::distanceTransformL2(f, f, indices, false);

			// store to array
			for (int y = 0; y < height; y++) {
				for (int x = 0; x < width; x++) {
					grid[x][y] = f[x][y];
				}
			}
		}

		bool save(std::string filename) {
			std::vector<unsigned char> png;
			lodepng::State state; 
			char image[width * height * 4];

			float scale = 0;
			for (int y = 0; y < height; ++y) {
				for (int x = 0; x < width; ++x) {
					scale = std::max(grid[x][y], scale);
				}
			}
			scale *= 1.0 / 255.0;
			for (int y = 0; y < height; ++y) {
				for (int x = 0; x < width; ++x) {
					unsigned idx = 4 * y * width + 4 * x;
					// std::cout << (int)(grid[x][y] / scale) << " " << grid[x][y] / scale << std::endl;
					// image[idx + 2] = std::min(255, (int)grid[x][y]);
					// image[idx + 1] = std::min(255, (int)grid[x][y]);
					// image[idx + 0] = std::min(255, (int)grid[x][y]);
					image[idx + 2] = (int)(grid[x][y] / scale);
					image[idx + 1] = (int)(grid[x][y] / scale);
					image[idx + 0] = (int)(grid[x][y] / scale);
					image[idx + 3] = (char)255;
				}
			}
			unsigned error = lodepng::encode(png, reinterpret_cast<const unsigned char*> (image), width, height, state);
			if(!error) lodepng::save_file(png, filename);
			//if there's an error, display it
			if(error) std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
			return error;
		}

		int memory() {
			return width*height*sizeof(float);
		}
	};
	inline int positive_modulo(int i, int n) {
	   return (i % n + n) % n;
	}
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
			// if (x == 42 && y == 36) std::cout << "thing: " << x+d0*(y+d1*z) << std::endl;
			// if (x+d0*(y+d1*z) == 2922) std::cout << "(" << x << ","<< y << ","<< z << ")" << std::endl;
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

		T interp_xy(float x, float y, float z)const {
			// round theta, interpolate x and y values
			int pz = positive_modulo(z, d2);
			float x_l = fmod(x, 1.f);
			float y_l = fmod(y, 1.f);

			int px = (int)x; // floor of x
			int py = (int)y; // floor of y

			float one_minus_x_l = 1.0 - x_l;
			float one_minus_y_l = 1.0 - y_l;

			return at(px, py, pz) * one_minus_x_l * one_minus_y_l
			     + at(px+1, py, pz) * x_l * one_minus_y_l
			     + at(px, py+1, pz) * one_minus_x_l * y_l
			     + at(px+1, py+1, pz) * x_l * y_l;
		}

		T interp(float x, float y, float z) const{

			#if _USE_INTERPOLATION == 3
				// int pz = positive_modulo(z, d2);
				// prevent the fmod from returning a negative number
				while(z < 0) z += (float)d2;
				
				float z_l = fmod(z, 1.f);
				int pz = (int)z; // floor of y
				float one_minus_z_l = 1.0 - z_l;

				return interp_xy(x, y, pz) * z_l + interp_xy(x, y, pz+1) * one_minus_z_l;


			#elif _USE_INTERPOLATION == 2
				// round theta, interpolate x and y values

				int pz = positive_modulo(z, d2);
				float x_l = fmod(x, 1.f);
				float y_l = fmod(y, 1.f);

				int px = (int)x; // floor of x
				int py = (int)y; // floor of y

				float one_minus_x_l = 1.0 - x_l;
				float one_minus_y_l = 1.0 - y_l;

				return at(px, py, pz) * one_minus_x_l * one_minus_y_l
				     + at(px+1, py, pz) * x_l * one_minus_y_l
				     + at(px, py+1, pz) * one_minus_x_l * y_l
				     + at(px+1, py+1, pz) * x_l * y_l;

			#elif _USE_INTERPOLATION == 1
				while(z < 0) z += (float)d2;
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
			#else
				int xd = round(x);
				int yd = round(y);
				// int zd = ((int) (round(z))) % d2;
				int zd = positive_modulo(z, d2);

				if (xd < 0 || xd >= d0 || yd < 0 || yd >= d1) return MAX_COST * 2.0;
				else return at(xd,yd,zd);
			#endif
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
		std::string to_string() {return "{dir:"<<direction<<",stop:"<<stop<<"}"; } 
		std::string to_string2() {return direction<<" "<<stop; } 
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
			// actions.push_back(OmniAction<T>());
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

	template <typename T>
	struct NonHolonomicAction
	{
		T steering_angle;
		int throttle; // backwards: -1, stop: 0, forwards: 0
		NonHolonomicAction() : steering_angle(0), throttle(0) {}
		NonHolonomicAction(T steer, int t) : steering_angle(steer), throttle(t) {}
		std::string to_string() {
			std::stringstream ss;
			// ss << "1";
			ss<<"{\"a\":"<<steering_angle<<",\"t\":"<<throttle<<"}";
			return ss.str();
		} 
		std::string to_string2() {
			std::stringstream ss;
			// ss << "1";
			ss<<steering_angle<<" "<<throttle;
			return ss.str();
		} 
	};

	/*
	This class uses ackerman dynamics as the state space
	*/
	template <typename T>
	class NonHolonomicModel
	{
	public:
		NonHolonomicModel(T ss, int disc, T ms, T L) : step_size(ss), discretization(disc), max_steer(ms), wheelbase(L) {};
		~NonHolonomicModel() {};

		void memoize(int theta_disc) {
			memo = new T[3 * discretization];
			trigTable = new T[3 * theta_disc];

			T step = max_steer / discretization;

			// forward direction angular sweep
			T angle = max_steer;
			for (int i = 0; i < discretization; ++i) {
				int index = round(angle * (float)discretization / max_steer) - 1;
				T tanAngle = tan(angle);
				memo[3*index] = step_size * tanAngle / wheelbase;
				memo[3*index+1] = (1.0-cos(memo[3*index]))*wheelbase / tanAngle;
				memo[3*index+2] = wheelbase * sin(memo[3*index]) / tanAngle;
				angle -= step;
			}

			for (int i = 0; i < theta_disc; ++i)
			{
				T angle = (M_2PI * (T)(i) / (T)(theta_disc));
				trigTable[i*2] = cos(angle);
				trigTable[i*2+1] = sin(angle);
			}
		}

		std::vector<NonHolonomicAction<T> > enumerate() {
			std::vector<NonHolonomicAction<T> > actions;

			T step = max_steer / discretization;

			// forward direction angular sweep
			T angle = max_steer;
			for (int i = 0; i < discretization; ++i) {
				actions.push_back(NonHolonomicAction<T>(angle, 1));
				angle -= step;
			}
			// forward direction forwards
			actions.push_back(NonHolonomicAction<T>(0.0, 1));
			angle = -step;
			for (int i = 0; i < discretization; ++i) {
				actions.push_back(NonHolonomicAction<T>(angle, 1));
				angle -= step;
			}

			// backward direction angular sweep
			angle = max_steer;
			for (int i = 0; i < discretization; ++i) {
				actions.push_back(NonHolonomicAction<T>(angle, -1));
				angle -= step;
			}

			// backward direction straight
			actions.push_back(NonHolonomicAction<T>(0.0, -1));
			angle = -step;
			for (int i = 0; i < discretization; ++i) {
				actions.push_back(NonHolonomicAction<T>(angle, -1));
				angle -= step;
			}

			// null action - 0 throttle
			// actions.push_back(NonHolonomicAction<T>(0,0));

			return actions;
		};

		// in this case, z is not scaled, it is in indicies
		void applyFast(NonHolonomicAction<T> action, T x, T y, T z, int zi, T &ax, T &ay, T &az) {
		// void applyFast(NonHolonomicAction<T> action, T x, T y, T z, T &ax, T &ay, T &az) {
			if (action.throttle == 0) {
				ax = x; ay = y; az = z;
			} else {
				T dt = 0.0;
				T dy = 0.0;
				T dx = 0.0;
				if (action.steering_angle == 0) {
					dx = action.throttle * step_size;
				} else {
					int index = 3*(round(fabs(action.steering_angle) * (float)discretization / max_steer) - 1);

					// std::cout << "a: " << action.steering_angle << " i: " << index << " " << fabs(action.steering_angle) << std::endl;
					dt = memo[index];
					dy = memo[index+1];
					dx = memo[index+2];

					if (action.throttle < 0) {
						dt = -dt;
						dx = -dx;
					}

					if (action.steering_angle < 0) {
						dy = -dy;
						dt = -dt;
					}
				}

				T c = trigTable[zi*2];
				T s = trigTable[zi*2+1];

				// std::cout << c << " " << s << " " << dx << " " << dy << " " << dt << " " << std::endl;
				ax = c*dx - s*dy + x;
				ay = s*dx + c*dy + y;
				az = dt + z;
			}
		}

		// TODO this could reuse some memory
		void apply(NonHolonomicAction<T> action, T x, T y, T z, T &ax, T &ay, T &az) {
			if (action.throttle == 0) {
				ax = x; ay = y; az = z;
			} else {
				// compute local coordinate frame deltas
				T tanAngle = tan(action.steering_angle);

				T dt = 0.0;
				T dy = 0.0;
				T dx = 0.0;
				if (tanAngle > _EPSILON || tanAngle < -_EPSILON) {
					dt = action.throttle * step_size * tanAngle / wheelbase;
					dy = (1.0-cos(dt))*wheelbase / tanAngle;
					dx = wheelbase * sin(dt) / tanAngle;
				} else if (action.throttle > 0) {
					dx = step_size;
				} else if (action.throttle < 0) {
					dx = -step_size;
				}
				
				// convert the above to global coordinate frame
				T c = cos(z);
				T s = sin(z);

				// std::cout << c << " " << s << " " << dx << " " << dy << " " << dt << " " << std::endl;
				ax = c*dx - s*dy + x;
				ay = s*dx + c*dy + y;
				az = dt + z;
			}
		}

		T cost(NonHolonomicAction<T> action) { 
			// if (action.throttle == 0) return 0.1;
			if (action.throttle == 0) return 0.0;
			if (action.steering_angle == 0) return step_size;

			return step_size + DEFLECTION_PENALTY * step_size * fabs(action.steering_angle) / max_steer;

			// return step_size;// + CHANGE_PENALTY; 
		}

		T max_affect_distance() { return step_size + 0.5; }
	
	private:
		T step_size;
		T max_steer;
		T wheelbase;
		int discretization;

		T* memo;
		T* trigTable;
	};

	template <typename T>
	class EmpericalNonHolonomicModel
	{
	public:
		EmpericalNonHolonomicModel(T ss, int disc, T scale) : step_size(ss), discretization(disc), world_scale(scale) {};
		~EmpericalNonHolonomicModel() {};

		std::vector<NonHolonomicAction<T> > enumerate() {
			std::vector<NonHolonomicAction<T> > actions;

			T step = max_steer / discretization;

			// forward direction angular sweep
			T angle = max_steer;
			for (int i = 0; i < discretization; ++i) {
				actions.push_back(NonHolonomicAction<T>(angle, 1));
				angle -= step;
			}
			// forward direction forwards
			actions.push_back(NonHolonomicAction<T>(0.0, 1));
			angle = -step;
			for (int i = 0; i < discretization; ++i) {
				actions.push_back(NonHolonomicAction<T>(angle, 1));
				angle -= step;
			}

			// backward direction angular sweep
			angle = max_steer;
			for (int i = 0; i < discretization; ++i) {
				actions.push_back(NonHolonomicAction<T>(angle, -1));
				angle -= step;
			}

			// backward direction straight
			actions.push_back(NonHolonomicAction<T>(0.0, -1));
			angle = -step;
			for (int i = 0; i < discretization; ++i) {
				actions.push_back(NonHolonomicAction<T>(angle, -1));
				angle -= step;
			}

			return actions;
		};

		void memoize(int theta_disc) {
			memo = new T[3 * discretization];
			trigTable = new T[3 * theta_disc];

			T step = max_steer / discretization;

			// forward direction angular sweep
			T angle = max_steer;
			for (int i = 0; i < discretization; ++i) {
				int index = round(angle * (float)discretization / max_steer) - 1;
				
				// dt = action.throttle * step_size * tanAngle / wheelbase;
				// dy = (1.0-cos(dt))*wheelbase / tanAngle;
				// dx = wheelbase * sin(dt) / tanAngle;

				// // T tanAngle = tan(angle);
				// memo[3*index] = step_size * tanAngle / wheelbase;
				// memo[3*index+1] = (1.0-cos(memo[3*index]))*wheelbase / tanAngle;
				// memo[3*index+2] = wheelbase * sin(memo[3*index]) / tanAngle;
				T R = steering_arc_radius(angle);

				memo[3*index]   = step_size / R;
				memo[3*index+1] = (1.0-cos(memo[3*index]))*R;
				memo[3*index+2] = sin(memo[3*index])*R;
				
				angle -= step;
			}

			for (int i = 0; i < theta_disc; ++i)
			{
				T angle = (M_2PI * (T)(i) / (T)(theta_disc));
				trigTable[i*2] = cos(angle);
				trigTable[i*2+1] = sin(angle);
			}
		}

		// given a positive steering angle, returns the expected steering arc radius
		//    angle (0, ackerman_transition_begin) uses ackermann geometry
		//    angle (ackerman_transition_begin, ackerman_transition_end) interpolates ackermann geometry and emperical model
		//    angle (ackerman_transition_end, max_steer) uses emperical model
		//    angle (max_steer, inf) uses emperical model with max_steer as the angle
		T steering_arc_radius(T angle) {
			T unscaled_radius = 0.0;
			if (angle <= ackerman_transition_begin) {
				unscaled_radius = wheelbase / tan(angle);
			} else if (angle >= ackerman_transition_end) {
				angle = std::min(angle, max_steer);
				unscaled_radius = 1.0 / (angle * angle * poly_a + angle * poly_b + poly_c);
			} else {
				// interpolation value
				T t = (angle - ackerman_transition_begin) / (ackerman_transition_end - ackerman_transition_begin);
				unscaled_radius = (1.0 - t) * (wheelbase / tan(angle)) + t * (1.0 / (angle * angle * poly_a + angle * poly_b + poly_c));
			}

			// this case should never happen
			assert(unscaled_radius != 0.0);
			return unscaled_radius * world_scale;
		}

		// in this case, z is not scaled, it is in indicies
		void applyFast(NonHolonomicAction<T> action, T x, T y, T z, int zi, T &ax, T &ay, T &az) {
			if (action.throttle == 0) {
				ax = x; ay = y; az = z;
			} else {
				T dt = 0.0;
				T dy = 0.0;
				T dx = 0.0;
				if (action.steering_angle == 0) {
					dx = action.throttle * step_size;
				} else {
					int index = 3*(round(fabs(action.steering_angle) * (float)discretization / max_steer) - 1);

					// std::cout << "a: " << action.steering_angle << " i: " << index << " " << fabs(action.steering_angle) << std::endl;
					dt = memo[index];
					dy = memo[index+1];
					dx = memo[index+2];

					if (action.throttle < 0) {
						dt = -dt;
						dx = -dx;
					}

					if (action.steering_angle < 0) {
						dy = -dy;
						dt = -dt;
					}
				}

				T c = trigTable[zi*2];
				T s = trigTable[zi*2+1];

				// std::cout << c << " " << s << " " << dx << " " << dy << " " << dt << " " << std::endl;
				ax = c*dx - s*dy + x;
				ay = s*dx + c*dy + y;
				az = dt + z;
			}
		}

		void apply(NonHolonomicAction<T> action, T x, T y, T z, T &ax, T &ay, T &az) {
			if (action.throttle == 0) {
				ax = x; ay = y; az = z;
			} else {
				// compute local coordinate frame deltas
				T dt = 0.0;
				T dy = 0.0;
				T dx = 0.0;
				if (action.steering_angle > _EPSILON || action.steering_angle < -_EPSILON) {
					T R = steering_arc_radius(fabs(action.steering_angle));
					T sign = sgn(action.steering_angle);
					dt = sign * action.throttle * step_size / R;
					dy = sign * (1.0-cos(dt))*R;
					dx = sign * sin(dt)*R;
				} else {
					dx = action.throttle * step_size;
				}

				// convert the above to global coordinate frame
				T c = cos(z);
				T s = sin(z);

				// std::cout << c << " " << s << " " << dx << " " << dy << " " << dt << " " << std::endl;
				ax = c*dx - s*dy + x;
				ay = s*dx + c*dy + y;
				az = dt + z;
			}
		}

		T cost(NonHolonomicAction<T> action) { 
			if (action.throttle == 0) return 0.0;
			if (action.steering_angle == 0) return step_size;
			return step_size + DEFLECTION_PENALTY * step_size * fabs(action.steering_angle) / max_steer;
		}

		T max_affect_distance() { return step_size + 0.5; }
	
	private:
		// variables
		T world_scale;
		T step_size;
		int discretization;

		// model constants
		T max_steer = 0.335;
		T wheelbase = 0.23;
		T poly_a = -4.035;
		T poly_b =  5.153;
		T poly_c = -0.018;
		T ackerman_transition_begin = 0.095;
		T ackerman_transition_end = 0.105;
		
		T* memo;
		T* trigTable;
	};

	template <typename state_T, typename action_T, typename action_model_T>
	class AbstractValueIterator
	{
	public:
		AbstractValueIterator(std::vector<int> &d, std::vector<float> &min_v, std::vector<float> &max_v) : dims(d), min_vals(min_v), max_vals(max_v) {
			J = new Array3D<state_T>(dims[0], dims[1], dims[2]);
			J->fill(MAX_COST);
			J_prime = new Array3D<state_T>(dims[0], dims[1], dims[2]);
			
			Policy = new Array3D<action_T >(dims[0], dims[1], dims[2]);
			Policy->fill(action_T()); // null action
		};
		~AbstractValueIterator() {};

		state_T cost_at(state_T x, state_T y, state_T z) {
			if (x < 0 || x > dims[0] - 1 || y < 0 || y > dims[1] - 1) {
				// out of bounds cost
				return MAX_COST * 2.0;
			} else {
				return J->interp(x,y,z);
			}
		}

		void set_action_model(action_model_T * am) { action_model = am; }

		// run a single step of value iteration
		virtual void step() {
			std::vector<action_T > actions = action_model->enumerate();

			// preallocate space for the various variables for speed
			state_T sx, sy, sz;
			state_T ax, ay, az;
			state_T ux, uy, uz;
			int i, min_i;
			state_T min_cost;
			state_T alt_cost;
			J_prime->fill(MAX_COST);

			bool debug = 0;
			int x,y,z;

			// if (debug && y > 8 && y < 12) std::cout << "HELLO" << std::endl;
			for (x = 0; x < dims[0]; ++x) {
				for (y = 0; y < dims[1]; ++y) {
					for (z = 0; z < dims[2]; ++z) {
						if (is_goal(x,y,z)) {
							// std::cout << "goal state: " << x << " " << y << " " << z << std::endl;
							// short circuit the action model if this is a goal state
							Policy->at(x,y,z) = action_T(); // null action
							J_prime->at(x,y,z) = terminal_cost(x,y,z);
						} else {
							min_i = -1;
							min_cost = MAX_COST*10.0;
							// if (debug && y > 13 && y < 15) std::cout << "state: ("<<x<<", "<<y<<", "<<z<<"), cost: " << min_cost << std::endl;
							// perform one step of value iteration
							for (i = 0; i < actions.size(); ++i) {
								action_model->apply(actions[i], x, y, RESCALE_T(z, dims[2]), ax, ay, az);
								az = UNSCALE_T(az, dims[2]);
								
								alt_cost = DISCOUNT * cost_at(ax, ay, az) + action_model->cost(actions[i]);
								// if (debug && y > 13 && y < 15) std::cout << "     + "<<debug_action(actions[i])<<" --> ("<<ax<< ", "<<ay<<", "<<az<<") @ "<<alt_cost;
								if (alt_cost < min_cost-_EPSILON) {
									min_cost = alt_cost;
									min_i = i;
									// if (debug && y > 13 && y < 15) std::cout << "***";
								}
								// if (debug && y > 13 && y < 15) std::cout <<std::endl;
							}
							
							// This should never happen - the algorithm must always make a choice
							assert(min_i >= 0);

							// modify the policy and the cost to reflect the iteration
							Policy->at(x,y,z) = actions[min_i];
							J_prime->at(x,y,z) = min_cost;
						}
					}
				}
			}
			J->replace(J_prime);
		}

		// return true if the given state is a goal state
		virtual bool is_goal(state_T x, state_T y, state_T z) {
			return false;
		}

		virtual state_T terminal_cost(state_T x, state_T y, state_T z) {
			return 0.0;
		}

		virtual std::string debug_action(action_T action) { return " "; }

		bool save_policy(int ind, std::string filename) {
			std::vector<unsigned char> png;
			lodepng::State state; //optionally customize this one
			int width = this->dims[0];
			int height = this->dims[1];
			char image[width * height * 4];

			for (int y = 0; y < height; ++y) {
				for (int x = 0; x < width; ++x) {
					unsigned idx = 4 * y * width + 4 * x;

					std::tuple<int, int, int> color = get_policy_color(x,y,ind);
					image[idx + 2] = (char)std::get<2>(color);
					image[idx + 1] = (char)std::get<1>(color);
					image[idx + 0] = (char)std::get<0>(color);
					image[idx + 3] = (char)255;
				}
			}
			unsigned error = lodepng::encode(png, reinterpret_cast<const unsigned char*> (image), width, height, state);
			if(!error) lodepng::save_file(png, filename);
			//if there's an error, display it
			if(error) std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
			return error;
		}

		bool save_slice(int ind, std::string filename) {
			int width = dims[0];
			int height = dims[1];
			float max_val = 0.0;
			for (int y = 0; y < height; ++y) {
				for (int x = 0; x < width; ++x) {
					float val = J->at(x,y,ind);
					if (val < 1000 and val > max_val)
						max_val = val;
				}
			}
			save_slice(ind, max_val, filename);
		}
		bool save_slice(int ind, float max_val, std::string filename) {
			std::vector<unsigned char> png;
			lodepng::State state; //optionally customize this one
			int width = dims[0];
			int height = dims[1];
			char image[width * height * 4];
			for (int y = 0; y < height; ++y) {
				float val;
				for (int x = 0; x < width; ++x) {
					unsigned idx = 4 * y * width + 4 * x;
					val = J->at(x,y,ind);

					if (val > 1000) {
						image[idx + 2] = (char)255;
						image[idx + 1] = (char)255;
						image[idx + 0] = (char)255;
						image[idx + 3] = (char)255;
					} else {
						
						val = 175.0 * val / max_val;
						image[idx + 2] = (char)val;
						image[idx + 1] = (char)val;
						image[idx + 0] = (char)val;
						image[idx + 3] = (char)255;
					}
				}
				// std::cout << val << "  " << J->at(200,y,ind) << std::endl;
				
			}
			unsigned error = lodepng::encode(png, reinterpret_cast<const unsigned char*> (image), width, height, state);
			if(!error) lodepng::save_file(png, filename);
			//if there's an error, display it
			if(error) std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
			return error;
		}

		bool serialize_cost(std::string filename) {
			std::cout << "serializing cost to file system" << std::endl;
			std::stringstream serialized;

			serialized << "{\"cost\": {" << std::endl;
			serialized << T1 << "\"data\": [";
			for (int i = 0; i < J->d0 * J->d1 * J->d2; ++i)
			{
				if (i > 0) serialized <<  ", ";
				serialized << J->pool[i]; 
			}
			serialized << "]," << std::endl;
			serialized << T1 << "\"d0\": " << J->d0 << "," << std::endl;
			serialized << T1 << "\"d1\": " << J->d1 << "," << std::endl;
			serialized << T1 << "\"d2\": " << J->d2 << std::endl;
			serialized << "}}" << std::endl;


			std::ofstream file;  
			file.open(filename);
			file << serialized.str();
			file.close();
		}

		bool load_map(std::string filename) {
			unsigned error;
			unsigned char* image;

			unsigned int width, height;

			error = lodepng_decode32_file(&image, &width, &height, filename.c_str());
			if(error) {
				printf("ERROR %u: %s\n", error, lodepng_error_text(error));
				return false;
			}

			ObstacleMap = new Array3D<bool>(dims[0], dims[1], 1);
			GoalMap = new Array3D<bool>(dims[0], dims[1], 1);

			std::cout << height << width << std::endl;
			std::cout << dims[0] << dims[1] << std::endl;

			for (int y = 0; y < height; ++y) {
				for (int x = 0; x < width; ++x) {
					unsigned idx = 4 * y * width + 4 * x;
					int r = image[idx + 2];
					int g = image[idx + 1];
					int b = image[idx + 0];

					if (r == 0 && g == 255 && b == 0)
						GoalMap->at(x, y, 0) = 1;
					else
						GoalMap->at(x, y, 0) = 0;

					if (r == 0 && g == 0 && b == 0)
						ObstacleMap->at(x, y, 0) = 1;
					else
						ObstacleMap->at(x, y, 0) = 0;
				}
			}
		}

		bool serialize_policy(std::string filename) {
			std::cout << "serializing to file system" << std::endl;
			std::stringstream serialized;

			serialized << "{\"policy\": {" << std::endl;
			serialized << T1 << "\"data\": [";
			for (int i = 0; i < J->d0 * J->d1 * J->d2; ++i)
			{
				if (i > 0) serialized <<  ", ";
				serialized << Policy->pool[i].to_string(); 
			}
			serialized << "]," << std::endl;
			serialized << T1 << "\"d0\": " << J->d0 << "," << std::endl;
			serialized << T1 << "\"d1\": " << J->d1 << "," << std::endl;
			serialized << T1 << "\"d2\": " << J->d2 << std::endl;
			serialized << "}}" << std::endl;


			std::ofstream file;  
			file.open(filename);
			file << serialized.str();
			file.close();
		}

		// override this function to make save_policy give something interesting
		virtual std::tuple<int, int, int> get_policy_color(int x, int y, int t) { return std::make_tuple(0, 0, 0); }

		Array3D<state_T>* getJ() { return J; }
		Array3D<state_T>* getPolicy() { return Policy; }
	protected:
		Array3D<state_T>* J_prime;
		Array3D<state_T>* J;
		Array3D<action_T >* Policy;

		Array3D<bool>* ObstacleMap;
		Array3D<bool>* GoalMap;
		
		action_model_T * action_model;
		std::vector<int> dims;
		std::vector<float> min_vals;
		std::vector<float> max_vals;
	};

	template <typename state_T>
	class OmniValueIterator : public AbstractValueIterator<state_T, OmniAction<state_T>, OmniActionModel<state_T> >
	{
	public:
		OmniValueIterator(std::vector<int> &d, std::vector<float> &min_v, std::vector<float> &max_v) 
			: AbstractValueIterator<state_T, OmniAction<state_T>, OmniActionModel<state_T> >(d, min_v, max_v) {}
		
		std::tuple<int, int, int> get_policy_color(int x, int y, int t) {
			float direction = this->Policy->at(x,y,t).direction;
			if (direction < 0)
				return std::make_tuple(255, 0, 0);
			int val = (255.0 * direction / M_2PI);
			return std::make_tuple(val, val, val);
		}

		bool is_goal(state_T x, state_T y, state_T z) {
			return this->dims[0]/2-2 <= x && x < this->dims[0]/2+2 
			    && this->dims[1]/2-2 <= y && y < this->dims[1]/2+2;
		}

		state_T terminal_cost(state_T x, state_T y, state_T z) {
			return abs(this->dims[0]/2 - x) + abs(this->dims[1]/2 - y);
		}
	};

	template <typename state_T, typename actionModel_T>
	class NonHolonomicIterator : public AbstractValueIterator<state_T, NonHolonomicAction<state_T>, actionModel_T >
	{
	public:
		NonHolonomicIterator(std::vector<int> &d, std::vector<float> &min_v, std::vector<float> &max_v) 
			: AbstractValueIterator<state_T, NonHolonomicAction<state_T>, actionModel_T >(d, min_v, max_v) {}

		virtual std::string debug_action(NonHolonomicAction<state_T> action) { return " |" + std::to_string(action.steering_angle) + "  " + std::to_string(action.throttle) + "| "; }
			
		bool is_goal(state_T x, state_T y, state_T z) {
			float max_theta_err = 0.18;
			float theta = M_2PI * (float)z / (float)this->dims[2];
			return this->dims[0]/2-8 <= x && x < this->dims[0]/2+8 
			    && this->dims[1]/2-3 <= y && y < this->dims[1]/2+3
			    && (theta <= max_theta_err || M_2PI - theta <= max_theta_err);
		}

		std::tuple<int, int, int> get_policy_color(int x, int y, int t) {
			float throttle = this->Policy->at(x,y,t).throttle;
			if (throttle == 0)
				return std::make_tuple(255, 0, 0);
			if (throttle > 0) {
				float angle = this->Policy->at(x,y,t).steering_angle;
				int val = (75.0 * angle / 0.35);
				// return std::make_tuple(0, val, 0);
				return std::make_tuple(0, 175 + val, 0);
			} else {
				float angle = this->Policy->at(x,y,t).steering_angle;
				int val = (75.0 * angle / 0.35);
				// return std::make_tuple(0, 0, val);
				return std::make_tuple(0, 0, 175 + val);
			}
		}
	};

	template <typename state_T, typename actionModel_T>
	class NonHolonomicIteratorLargeSpace : public NonHolonomicIterator<state_T, actionModel_T>
	{
	public:
		NonHolonomicIteratorLargeSpace(std::vector<int> &d, std::vector<float> &min_v, std::vector<float> &max_v) 
			: NonHolonomicIterator<state_T, actionModel_T>(d, min_v, max_v) {
				// these are for backwards motion
				J2 = new Array3D<state_T>(this->dims[0], this->dims[1], this->dims[2]);
				J2->fill(MAX_COST);
				J2_prime = new Array3D<state_T>(this->dims[0], this->dims[1], this->dims[2]);

				Policy2 = new Array3D<NonHolonomicAction<state_T> >(this->dims[0], this->dims[1], this->dims[2]);
				Policy2->fill(NonHolonomicAction<state_T>()); // null action

				#if _USE_HEURISTIC_PRUNING
				usefulMap = new OMap(this->dims[0], this->dims[1]);
				usefulDist = new DistanceTransform(this->dims[0], this->dims[1]);

				usefulDist->fill(0.0);
				#endif
			}

		bool load_3D_omap(std::string summary) {

			std::ifstream summary_file(summary);
			std::string line;

			int discretization = 0;
			std::vector<std::string> files;
			if (summary_file.is_open()) {
				int num_left = 0;
				bool collecting_file_list = false;
				while ( getline (summary_file,line) ) {

					// parse out the discretization
					if (line.find("discretization:") != std::string::npos) discretization = std::stoi(split(line, ':')[1]);
					
					if (line.find("files") != std::string::npos) {
						collecting_file_list = discretization > 0;
						num_left = discretization;
						continue;
					} 

					if (collecting_file_list) {
						files.push_back(split(line, '-')[1]);
						num_left -= 1;
						if (num_left == 0) collecting_file_list = false;
					}
				}
				summary_file.close();
			} else {
				std::cout << "Unable to open file" << std::endl;
				return false;
			} 

			if (discretization <= 0) {
				std::cout << "Relevant data not found in file." << std::endl;
				return false;
			}

			// this->ObstacleMap = new Array3D<bool>();
			std::string filename;
			std::cout << "Loading " << discretization << " slices of occupancy map" << std::endl;
			for (int i = 0; i < files.size(); ++i)
			{
				filename = files[i];
				std::cout << " - Loading: " << files[i] << std::endl;

				unsigned error;
				unsigned char* image;

				unsigned int width, height;

				error = lodepng_decode32_file(&image, &width, &height, filename.c_str());
				if(error) {
					printf("ERROR %u: %s\n", error, lodepng_error_text(error));
					return false;
				}

				if (i == 0) {
					std::cout << " - Allocating space in the ObstacleMap: " << width << " * " << height << " * " << discretization << std::endl;
					this->ObstacleMap = new Array3D<bool>(width, height, discretization);
				}

				for (int y = 0; y < height; ++y) {
					for (int x = 0; x < width; ++x) {
						unsigned idx = 4 * y * width + 4 * x;
						int r = image[idx + 2];
						int g = image[idx + 1];
						int b = image[idx + 0];

						if (r == 0 && g == 0 && b == 0)
							this->ObstacleMap->at(x, y, i) = 1;
						else
							this->ObstacleMap->at(x, y, i) = 0;
					}
				}
			}
		}

		// save all information necessary to execute the policy
		//   - raw files containing forwards and backwards policies
		//   - metadata file
		void save_full_policy(std::string path) {
			std::cout << "Serializing full policy to: " << path << std::endl;
			std::stringstream policy1;
			std::stringstream policy2;
			std::stringstream cost1;
			std::stringstream cost2;
			std::stringstream summary;

			char *full_path = realpath(path.c_str(), NULL);

			summary << "{\"policy\": {" << std::endl;
			summary << T1 << "\"d0\": " << J2->d0 << "," << std::endl;
			summary << T1 << "\"d1\": " << J2->d1 << "," << std::endl;
			summary << T1 << "\"d2\": " << J2->d2 << "," << std::endl;
			summary << T1 << "\"forwards_policy\": \"" << full_path << "/policy_forwards.object\"," << std::endl;
			summary << T1 << "\"backwards_policy\": \"" << full_path << "/policy_backwards.object\"," << std::endl;
			summary << T1 << "\"forwards_cost\": \"" << full_path << "/cost_forwards.object\"," << std::endl;
			summary << T1 << "\"backwards_cost\": \"" << full_path << "/cost_backwards.object\"" << std::endl;
			summary << "}}" << std::endl;


			for (int i = 0; i < J2->d0 * J2->d1 * J2->d2; ++i)
			{
				policy1 << this->Policy->pool[i].to_string2() << "\n"; 
				policy2 << Policy2->pool[i].to_string2() << "\n"; 

				cost1 << this->J->pool[i] << "\n"; 
				cost2 << J2->pool[i] << "\n"; 
			}

			std::ofstream file;  
			file.open(path + "/summary.json");
			file << summary.str();
			file.close();

			std::ofstream file2;  
			file2.open(path + "/policy_forwards.object");
			file2 << policy1.str();
			file2.close();

			std::ofstream file3;  
			file3.open(path + "/policy_backwards.object");
			file3 << policy2.str();
			file3.close();

			std::ofstream file4;  
			file4.open(path + "/cost_forwards.object");
			file4 << cost1.str();
			file4.close();

			std::ofstream file5;  
			file5.open(path + "/cost_backwards.object");
			file5 << cost2.str();
			file5.close();

			free(full_path);
		}

		bool within_err(int z, float target, float err) {
			float diff = M_2PI * (float) z / (float)this->dims[2] - target;
			while (diff < 0) diff += M_2PI;
			return diff <= err || M_2PI - diff <= err;
		}

		bool is_goal(int x, int y, int z) {
			if (x >= this->GoalMap->d0 || y >= this->GoalMap->d1) return false;
			float max_theta_err = 0.25;
			return this->GoalMap->at(x, y, 0) == 1 && (within_err(z, 0.0, max_theta_err) || within_err(z, 3.1415, max_theta_err));

			// return this->GoalMap->at(x, y, 0) == 1 && within_err(z, 0.0, max_theta_err);
			// return this->GoalMap->at(x, y, 0) == 1 && within_err(z, 3.1415, max_theta_err);
		}

		bool is_obstacle(state_T x, state_T y, state_T z) {
			int px = round(x);
			int py = round(y);
			if (px >= this->ObstacleMap->d0 || py >= this->ObstacleMap->d1) return false;

			// std::cout << x << "  " << y << "  " << round(x) << "  " <<  round(y) << std::endl;
			if (this->ObstacleMap->d2 > 1) {
				state_T theta = RESCALE_T(z, this->dims[2]);
				int omap_slice = round(UNSCALE_T(theta, this->ObstacleMap->d2));
				omap_slice = positive_modulo(omap_slice, this->ObstacleMap->d2);
				return this->ObstacleMap->at(round(x), round(y), omap_slice) == 1;
			} else {
				return this->ObstacleMap->at(px, py, 0) == 1;
			}
		}
		// bool is_obstacle(int x, int y) {
		// 	int px = round(x);
		// 	int py = round(y);
		// 	if (x >= this->ObstacleMap->d0 || y >= this->ObstacleMap->d1) return false;
		// 	return this->ObstacleMap->at(px, py, 0) == 1;
		// }


		state_T cost_at2(state_T x, state_T y, state_T z) {
			if (x < 0 || x > this->dims[0] - 1 || y < 0 || y > this->dims[1] - 1 || is_obstacle(x,y,z)) {
			// if (x < 0 || x > this->dims[0] - 1 || y < 0 || y > this->dims[1] - 1 || is_obstacle(x,y)) {
				// out of bounds cost
				return MAX_COST * 2.0;
			} else {
				return J2->interp(x,y,z);
			}
		}

		state_T cost_at(state_T x, state_T y, state_T z) {
			if (x < 0 || x > this->dims[0] - 1 || y < 0 || y > this->dims[1] - 1 || is_obstacle(x,y,z)) {
			// if (x < 0 || x > this->dims[0] - 1 || y < 0 || y > this->dims[1] - 1 || is_obstacle(x,y)) {
				// out of bounds cost
				return MAX_COST * 2.0;
			} else {
				return this->J->interp(x,y,z);
			}
		}

		bool serialize_policy2(std::string filename) {
			std::cout << "serializing to file system" << std::endl;
			std::stringstream serialized;

			serialized << "{\"policy\": {" << std::endl;
			serialized << T1 << "\"data\": [";
			for (int i = 0; i < J2->d0 * J2->d1 * J2->d2; ++i)
			{
				if (i > 0) serialized <<  ", ";
				serialized << Policy2->pool[i].to_string(); 
			}
			serialized << "]," << std::endl;
			serialized << T1 << "\"d0\": " << J2->d0 << "," << std::endl;
			serialized << T1 << "\"d1\": " << J2->d1 << "," << std::endl;
			serialized << T1 << "\"d2\": " << J2->d2 << std::endl;
			serialized << "}}" << std::endl;


			std::ofstream file;  
			file.open(filename);
			file << serialized.str();
			file.close();
		}

		bool useful_update(state_T old_cost, state_T new_cost) {
			return (old_cost - new_cost > 0.5 || old_cost - new_cost < -0.5) && new_cost < 1000;
		}


		// run a single step of value iteration
		void step() {
			std::vector<NonHolonomicAction<state_T> > actions = this->action_model->enumerate();

			// preallocate space for the various variables for speed
			state_T sx, sy, sz;
			state_T ax, ay, az;
			state_T ux, uy, uz;
			int i, min_i;
			state_T min_cost;
			state_T alt_cost;
			// this->J_prime->fill(MAX_COST);

			this->J_prime->replace(this->J);
			this->J2_prime->replace(this->J2);

			bool debug = 0;
			int x,y,z;

			#if _USE_HEURISTIC_PRUNING == 1
			// clear binary usefulness buffer
			usefulMap->fill(false);
			#endif

			

			// good = 0;
			// bad = 0;

			was_useful = false;

			// if (debug && y > 8 && y < 12) std::cout << "HELLO" << std::endl;
			for (x = 0; x < this->dims[0]; ++x) {
				for (y = 0; y < this->dims[1]; ++y) {
					#if _USE_HEURISTIC_PRUNING == 1
					if (usefulDist->grid[x][y] > this->action_model->max_affect_distance()) continue;
					#endif

					for (z = 0; z < this->dims[2]; ++z) {
						if (this->is_goal(x,y,z)) {
							if (iter > 0) continue;
							// short circuit the action model if this is a goal state
							this->Policy->at(x,y,z) = NonHolonomicAction<state_T>(); // null action
							this->Policy2->at(x,y,z) = NonHolonomicAction<state_T>(); // null action
							this->J_prime->at(x,y,z) = this->terminal_cost(x,y,z);
							this->J2_prime->at(x,y,z) = this->terminal_cost(x,y,z);

							if (iter == 0) {
								#if _USE_HEURISTIC_PRUNING == 1
								usefulMap->grid[x][y] = true;
								#endif
								good ++;
								was_useful = true;
							} else {
								bad ++;
							}
						} else {
							// value iteration for the forwards direction
							min_i = -1;
							min_cost = MAX_COST*10.0;
							// if (debug && y > 13 && y < 15) std::cout << "state: ("<<x<<", "<<y<<", "<<z<<"), cost: " << min_cost << std::endl;
							// perform one step of value iteration
							for (i = 0; i < actions.size(); ++i) {
								#if _USE_FAST_APPLY == 1
								this->action_model->applyFast(actions[i], x, y, RESCALE_T(z, this->dims[2]), z, ax, ay, az);
								az = UNSCALE_T(az, this->dims[2]);
								#else
								this->action_model->apply(actions[i], x, y, RESCALE_T(z, this->dims[2]), ax, ay, az);
								az = UNSCALE_T(az, this->dims[2]);
								#endif

								if (actions[i].throttle > 0) {
									// forwards motion, unpenalized
									alt_cost = DISCOUNT * this->cost_at(ax, ay, az) + this->action_model->cost(actions[i]);
								} else {
									// backwards motion, penalized
									alt_cost = DISCOUNT * this->cost_at2(ax, ay, az) + this->action_model->cost(actions[i]) + REVERSAL_PENALTY;
								}

								// if (debug && y > 13 && y < 15) std::cout << "     + "<<debug_action(actions[i])<<" --> ("<<ax<< ", "<<ay<<", "<<az<<") @ "<<alt_cost;
								if (alt_cost < min_cost-_EPSILON) {
									min_cost = alt_cost;
									min_i = i;
									// if (debug && y > 13 && y < 15) std::cout << "***";
								}
								// if (debug && y > 13 && y < 15) std::cout <<std::endl;
							}
							
							// This should never happen - the algorithm must always make a choice
							assert(min_i >= 0);

							// modify the policy and the cost to reflect the iteration
							this->Policy->at(x,y,z) = actions[min_i];
							this->J_prime->at(x,y,z) = min_cost;

							// if (x == 39 && y == 34) std::cout << "SHOULD NEVER HAPPEN" << std::endl;

							if (useful_update(this->J->at(x,y,z), min_cost)) {
								#if _USE_HEURISTIC_PRUNING == 1
								usefulMap->grid[x][y] = true;
								#endif
								was_useful = true;
								good ++;
							} else {
								bad ++;
							}

							// value iteration for the backwards motion direction
							min_i = -1;
							min_cost = MAX_COST*10.0;
							// if (debug && y > 13 && y < 15) std::cout << "state: ("<<x<<", "<<y<<", "<<z<<"), cost: " << min_cost << std::endl;
							// perform one step of value iteration
							for (i = 0; i < actions.size(); ++i) {
								#if _USE_FAST_APPLY == 1
								this->action_model->applyFast(actions[i], x, y, RESCALE_T(z, this->dims[2]), z, ax, ay, az);
								az = UNSCALE_T(az, this->dims[2]);
								#else
								this->action_model->apply(actions[i], x, y, RESCALE_T(z, this->dims[2]), ax, ay, az);
								az = UNSCALE_T(az, this->dims[2]);
								#endif
								
								if (actions[i].throttle > 0) {
									// forwards motion, penalized
									alt_cost = DISCOUNT * this->cost_at(ax, ay, az) + this->action_model->cost(actions[i]) + REVERSAL_PENALTY;
								} else {
									// backwards motion, unpenalized
									alt_cost = DISCOUNT * this->cost_at2(ax, ay, az) + this->action_model->cost(actions[i]);
								}

								// if (debug && y > 13 && y < 15) std::cout << "     + "<<debug_action(actions[i])<<" --> ("<<ax<< ", "<<ay<<", "<<az<<") @ "<<alt_cost;
								if (alt_cost < min_cost-_EPSILON) {
									min_cost = alt_cost;
									min_i = i;
									// if (debug && y > 13 && y < 15) std::cout << "***";
								}
								// if (debug && y > 13 && y < 15) std::cout <<std::endl;
							}
							
							// This should never happen - the algorithm must always make a choice
							assert(min_i >= 0);

							// modify the policy and the cost to reflect the iteration
							this->Policy2->at(x,y,z) = actions[min_i];
							this->J2_prime->at(x,y,z) = min_cost;

							if (useful_update(this->J2->at(x,y,z), min_cost)) {
								#if _USE_HEURISTIC_PRUNING == 1
								usefulMap->grid[x][y] = true;
								#endif
								was_useful = true;
								good ++;
							} else {
								bad ++;
							}
						}
					}
				}
			}

			iter ++;

			#if _USE_HEURISTIC_PRUNING == 1
			usefulDist->init(usefulMap);
			usefulMap->save("./useful_sequence/" + padded(iter,3) + ".png");
			#endif

			this->J->replace(this->J_prime);
			this->J2->replace(this->J2_prime);

			std::cout << "    g: " << good  << "  b: " << bad << std::endl;
		}

		bool last_step_useful() { return was_useful; }

	protected:
		Array3D<state_T>* J2_prime;
		Array3D<state_T>* J2;
		Array3D<NonHolonomicAction<state_T> >* Policy2;

		#if _USE_HEURISTIC_PRUNING
		OMap * usefulMap;
		DistanceTransform * usefulDist;
		#endif

		int good = 0;
		int bad = 0;
		int iter = 0;

		// whether or not the last step was useful
		bool was_useful = true;
	};
}



#endif