#pragma once

#include "utilities.hpp"

class LBM;
struct Info { // contains redundant information for console printing
	LBM* lbm = nullptr;
	bool allow_rendering = false; // allows interactive redering if true
	double runtime_lbm=0.0, runtime_total=0.0f; // lbm (compute) and total (compute + rendering + data evaluation) runtime
	double runtime_lbm_timestep_last=1.0, runtime_lbm_timestep_smooth=1.0, runtime_lbm_last=0.0; // for printing simulation info
	Clock clock; // for measuring total runtime
	ulong steps=max_ulong, steps_last=0ull; // runtime_lbm_last and steps_last are there if multiple run() commands are executed consecutively
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