#pragma once

#include "utilities.hpp"

class LBM;
struct Info { // contains redundant information for console printing
	LBM* lbm = nullptr;
	bool allow_rendering = false; // allows interactive redering if true
	double dt=1.0, dt_smooth=1.0, runtime=0.0, runtime_last=0.0; // for printing simulation info
	ulong steps=max_ulong, steps_last=0ull; // runtime_last and steps_last are there if multiple run() commands are executed consecutively
	uint host_allocation=17u, device_transfer=0u; // CPU memory allocation and GPU transfer per LBM step in Byte/cell
	uint cpu_mem_required=0u, gpu_mem_required=0u; // all in MB
	string collision = "";
	void initialize(LBM* lbm);
	void append(const ulong steps, const ulong t);
	void update(const double dt);
	double time() const; // returns either elapsed time or remaining time
	void print_logo() const;
	void print_initialize(); // enables interactive rendering
	void print_update() const;
	void print_finalize(); // disables interactive rendering
};
extern Info info; // declared in info.cpp