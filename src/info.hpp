#pragma once

#include "utilities.hpp"
#include <mutex>

class LBM;
struct Info { // contains redundant information for console printing
	LBM* lbm = nullptr;
	double runtime_lbm=0.0, runtime_total=0.0f, runtime_total_last=0.0; // lbm (compute) and total (compute + rendering + data evaluation) runtime
	double runtime_lbm_timestep_last=1.0, runtime_lbm_timestep_smooth=1.0; // for printing simulation info
	Clock clock; // for measuring total runtime
	ulong steps=max_ulong, steps_last=0ull; // runtime_total_last and steps_last are there if multiple run() commands are executed consecutively
	uint cpu_mem_required=0u, gpu_mem_required=0u; // all in MB
	string collision = "";
	std::mutex allow_printing; // to prevent threading conflicts when continuously printing updates to console
	void append(const ulong steps, const ulong total_steps, const ulong t);
	void update(const double dt);
	double time() const; // returns either elapsed time or remaining time
	void print_logo() const;
	void print_initialize(LBM* lbm); // enables interactive rendering
	void print_update() const;
	void print_finalize(); // disables interactive rendering
};
extern Info info; // declared in info.cpp