#pragma once

#include "utilities.hpp"

class Units { // contains the 4 base units [m], [kg], [s], [K] for unit conversions and vtk output
private:
	float unit_m=1.0f, unit_kg=1.0f, unit_s=1.0f, unit_K=1.0f; // 1 lattice unit times [m]/[kg]/[s]/[K] is [meter]/[kilogram]/[second]/[Kelvin]
public:
	void set_m_kg_s(const float x, const float u, const float rho/*=1.0f*/, const float si_x, const float si_u, const float si_rho) { // length x, velocity u, density rho in both simulation and SI units
		unit_m = si_x/x; // length si_x = x*[m]
		unit_kg = si_rho/rho*cb(unit_m); // density si_rho = rho*[kg/m^3]
		unit_s = u/si_u*unit_m; // velocity si_u = u*[m/s]
		print_info("Unit Conversion: 1 cell = "+to_string(1000.0f*this->si_x(1.0f), 3u)+" mm, 1 s = "+to_string(this->t(1.0f))+" time steps");
	}
	void set_m_kg_s(const float m, const float kg, const float s) { // do unit conversion manually
		this->unit_m = m;
		this->unit_kg = kg;
		this->unit_s = s;
		print_info("Unit Conversion: 1 cell = "+to_string(1000.0f*this->si_x(1.0f), 3u)+" mm, 1 s = "+to_string(this->t(1.0f))+" time steps");
	}
	void set_m_kg_s_K(const float x, const float u, const float rho/*=1.0f*/, const float T/*=1.0f*/, const float si_x, const float si_u, const float si_rho, const float si_T) { // length x, velocity u, density rho, temperature T in both simulation and SI units
		unit_m = si_x/x; // length si_x = x*[m]
		unit_kg = si_rho/rho*cb(unit_m); // density si_rho = rho*[kg/m^3]
		unit_s = u/si_u*unit_m; // velocity si_u = u*[m/s]
		unit_K = si_T/T; // length si_T = T*[K]
		print_info("Unit Conversion: 1 cell = "+to_string(1000.0f*this->si_x(1.0f), 3u)+" mm, 1 s = "+to_string(this->t(1.0f))+" time steps");
	}
	void set_m_kg_s_K(const float m, const float kg, const float s, const float K) { // do unit conversion manually
		this->unit_m = m;
		this->unit_kg = kg;
		this->unit_s = s;
		this->unit_K = K;
		print_info("Unit Conversion: 1 cell = "+to_string(1000.0f*this->si_x(1.0f), 3u)+" mm, 1 s = "+to_string(this->t(1.0f))+" time steps");
	}

	// the following methods convert SI units into simulation units (have to be called after set_m_kg_s(...);)
	float x(const float si_x) const { return si_x/unit_m; } // length si_x = x*[m]
	float m(const float si_m) const { return si_m/unit_kg; } // mass si_m = m*[kg]
	ulong t(const float si_t) const { return to_ulong(si_t/unit_s); } // time si_t = t*[s]
	float frequency(const float si_frequency) const { return si_frequency*unit_s; } // frequency si_frequency = frequency*[1/s]
	float omega(const float si_omega) const { return si_omega*unit_s; } // frequency si_omega = omega/[s]
	float u(const float si_u) const { return si_u*unit_s/unit_m; } // velocity si_u = u*[m/s]
	float rho(const float si_rho) const { return si_rho*cb(unit_m)/unit_kg; } // density si_rho = rho*[kg/m^3]
	float Q(const float si_Q) const { return si_Q*unit_s/cb(unit_m); } // flow rate si_Q = Q*[m^3/s]
	float nu(const float si_nu) const { return si_nu*unit_s/sq(unit_m); } // kinematic shear viscosity si_nu = nu*[m^2/s]
	float mu(const float si_mu) const { return si_mu*unit_s*unit_m/unit_kg; } // dynamic shear viscosity si_mu = mu*[kg/(m*s)]
	float g(const float si_g) const { return si_g/unit_m*sq(unit_s); } // gravitational acceleration si_g = g*[m/s^2]
	float f(const float si_f) const { return si_f*sq(unit_m*unit_s)/unit_kg; } // force per volume si_f = f*[kg/(m*s)^2]
	float f(const float si_rho, const float si_g) const { return si_rho*si_g*sq(unit_m*unit_s)/unit_kg; } // force per volume f = rho*g = si_rho/[kg/m^3]*si_g/[m/s^2] = si_rho*si_g*[(m*s)^2/kg]
	float F(const float si_F) const { return si_F*sq(unit_s)/(unit_kg*unit_m); } // force si_F = F*[kg*m/s^2]
	float M(const float si_M) const { return si_M*sq(unit_s)/(unit_kg*sq(unit_m)); } // torque si_M = M*[kg*m^2/s^2]
	float sigma(const float si_sigma) const { return si_sigma*sq(unit_s)/unit_kg; } // surface tension si_sigma = sigma*[kg/s^2]
	float T(const float si_T) { return si_T/unit_K; } // temperature si_T = T*[K]
	float alpha(const float si_alpha) { return si_alpha*unit_s/sq(unit_m); } // thermal diffusion coefficient si_alpha = alpha*[m^2/s]
	float beta(const float si_beta) { return si_beta*unit_K; } // (volumetric) thermal expansion coefficient si_beta = beta*[1/K]

	// the following methods convert simulation units into SI units (have to be called after set_m_kg_s(...);)
	float si_x(const uint x) const { return (float)x*unit_m; } // length si_x = x*[m]
	float si_x(const float x) const { return x*unit_m; } // length si_x = x*[m]
	float si_m(const float m) const { return m*unit_kg; } // mass si_m = m*[kg]
	float si_t(const ulong t) const { return (float)t*unit_s; } // time si_t = t*[s]
	float si_frequency(const float frequency) const { return frequency/unit_s; } // frequency si_frequency = frequency*[1/s]
	float si_V(const float V) const { return V*cb(unit_m); } // volume si_V = V*[m^3]
	float si_u(const float u) const { return u*unit_m/unit_s; } // velocity si_u = u*[m/s]
	float si_rho(const float rho) const { return rho*unit_kg/cb(unit_m); } // density si_rho = rho*[kg/m^3]
	float si_p(const float p) const { return p*unit_kg/(unit_m*sq(unit_s)); } // pressure si_p = p*[kg/(m*s^2)]
	float si_Q(const float Q) const { return Q*cb(unit_m)/unit_s; } // flow rate si_Q = Q*[m^3/s]
	float si_nu(const float nu) const { return nu*sq(unit_m)/unit_s; } // kinematic shear viscosity si_nu = nu*[m^2/s]
	float si_g(const float g) const { return g*unit_m/sq(unit_s); } // gravitational acceleration si_g = g*[m/s^2]
	float si_f(const float f) const { return f*unit_kg/sq(unit_m*unit_s); } // force per volume si_f = f*[kg/(m*s)^2]
	float si_F(const float F) const { return F*unit_kg*unit_m/sq(unit_s); } // force si_F = F*[kg*m/s^2]
	float si_M(const float M) const { return M*unit_kg*sq(unit_m)/sq(unit_s); } // torque si_M = M*[kg*m^2/s^2]
	float si_sigma(const float sigma) const { return sigma*unit_kg/sq(unit_s); } // surface tension si_sigma = sigma*[kg/s^2]
	float si_T(const float T) { return T*unit_K; } // temperature si_T = T*[K]
	float si_alpha(const float alpha) { return alpha*sq(unit_m)/unit_s; } // thermal diffusion coefficient si_alpha = alpha*[m^2/s]
	float si_beta(const float beta) { return beta/unit_K; } // (volumetric) thermal expansion coefficient si_beta = beta*[1/K]

	// other conversions in simulation units (can be called before set_m_kg_s(...);)
	float Re(const float si_Re) const { return si_Re; } // Reynolds number Re = x*u/nu = [1] no unit
	float Re(const float x, const float u, const float nu) const { return x*u/nu; } // Reynolds number Re = x*u/nu = [1] no unit
	float Re(const float x, const float u, const float mu, const float rho) const { return x*u*rho/mu; } // Reynolds number Re = x*u*rho/mu = [1] no unit
	float We(const float x, const float u, const float rho, const float sigma) const { return x*sq(u)*rho/sigma; } // Weber number We = x*u^2*rho/sigma = [1] no unit
	float Fr(const float x, const float u, const float g) const { return u/sqrt(x*g); } // Froude number Fr = u/sqrt(x*g) = [1] no unit
	float Ca(const float u, const float mu, const float sigma) const { return u*mu/sigma; } // Capillary number Ca = u*mu/sigma = [1] no unit
	float Ca(const float u, const float rho, const float nu, const float sigma) const { return u*rho*nu/sigma; } // Capillary number Ca = u*rho*nu/sigma = [1] no unit
	float Bo(const float x, const float rho, const float g, const float sigma) const { return sq(x)*rho*g/sigma; } // Bond number Bo = x^2*rho*g/sigma = [1] no unit
	float Mo(const float rho, const float delta_rho, const float g, const float sigma, const float mu) const { return g*delta_rho*sq(sq(mu))/(cb(sigma)*sq(rho)); } // Morton number g*delta_rho*mu^4/(sigma^3*rho^2)
	float Ga(const float x, const float nu, const float g) { return cb(x)*g/sq(nu); } // Galilei number Ga = x^3*g/nu^2 = [1] no unit
	float Ga(const float x, const float rho, const float nu, const float f) { return cb(x)*f/(sq(nu)*rho); } // Galilei number Ga = x^3*g/nu^2 = [1] no unit
	float Ma(const float u) const { return u/0.57735027f; } // Mach number Ma = u/c = [1] no unit, c = 1/sqrt(3) is lattice speed of sound
	float rho_from_p(const float p) const { return 3.0f*p; } // density rho = p/c^2 = 3*p = [kg/(m*s^2)], p is pressure, c = 1/sqrt(3) is lattice speed of sound
	float rho_laplace(const float sigma, const float R) { return 6.0f*sigma/R; } // sphere laplace pressure p = 2*sigma/R, density rho = 3*p
	float rho_hydrostatic(const float f, const float z, const float h) { return 3.0f*f*(h-z)+1.0f; } // hydrostatic pressure p = rho*g*h+p0 = f*h+1/3, density rho = 3*p, force per volume f = rho*g, rho0 = 1
	float nu_from_mu(const float mu, const float rho) const { return mu/rho; } // kinematic shear viscosity nu = mu/rho = [m^2/s]
	float nu_from_tau(const float tau) const { return (tau-0.5f)/3.0f; } // kinematic shear viscosity nu = (tau-1/2)/3 = [m^2/s], LBM relaxation time tau
	float nu_from_Re(const float Re, const float x, const float u) const { return x*u/Re; } // kinematic shear viscosity nu = x*u/Re = [m^2/s]
	float f_from_F(const float F, const float V) const { return F/V; } // force per volume f = F/V = [kg/(m*s)^2]
	float f_from_g(const float g, const float rho) const { return rho*g; } // force per volume f = rho*g = [kg/(m*s)^2]
	float g_from_f(const float f, const float rho) const { return f/rho; } // force per volume f = rho*g = [kg/(m*s)^2]
	float u_from_Re(const float Re, const float x, const float nu) const { return Re*nu/x; } // velocity u = Re*nu/x, Reynolds number Re = x*u/nu = [1] no unit
	float u_from_Re(const float Re, const float x, const float mu, const float rho) const { return Re*mu/(x*rho); } // velocity u = Re*nu/x, Reynolds number Re = x*u*rho/mu = [1] no unit
	float u_from_Ma(const float Ma) const { return 0.57735027f*Ma; } // velocity u = c*Ma, Mach number Ma = u/c = [1] no unit, c = 1/sqrt(3) is lattice speed of sound
	float u_from_We(const float We, const float x, const float sigma, const float rho) const { return sqrt(We*sigma/(x*rho)); } // velocity u = sqrt(We*sigma/(x*rho)) = [m/s], Weber number We = x*u^2*rho/sigma = [1] no unit
	float u_from_Fr(const float Fr, const float x, const float g) const { return Fr*sqrt(x*g); } // velocity u = Fr*sqrt(x*g) = [m/s], Froude number Fr = u/sqrt(x*g) = [1] no unit
	float u_from_Ca(const float Ca, const float sigma, const float nu, const float rho) const { return Ca*sigma/(rho*nu); } // velocity u = Ca*sigma/(rho*nu) = [m/s], Capillary number Ca = u*rho*nu/sigma = [1] no unit
	float u_from_Ca(const float Ca, const float sigma, const float mu) const { return Ca*sigma/mu; } // velocity u = Ca*sigma/(rho*nu) = [m/s], Capillary number Ca = u*mu/sigma = [1] no unit
	float u_from_f_Poiseuille_2D(const float f, const float rho, const float nu, const float R) const { return f*sq(R)/(2.0f*rho*nu); } // center velocity in 2D Poiseuille flow in channel u = f/(2*rho*nu)*R^2, force per volume f = 2*u*rho*nu/R^2
	float u_from_f_Poiseuille_3D(const float f, const float rho, const float nu, const float R) const { return f*sq(R)/(4.0f*rho*nu); } // center velocity in 3D Poiseuille flow in channel u = f/(4*rho*nu)*R^2, force per volume f = 4*u*rho*nu/R^2
	float u_from_f_Poiseuille_2D(const float Q, const float R) const { return 0.75f*Q/R; } // center velocity in 2D Poiseuille flow in channel u = 3/4*Q/R, force per volume f = 2*u*rho*nu/R^2, flow rate Q = 2/3*f*R^3/(rho*nu) = 4/3*u*R
	float u_from_f_Poiseuille_3D(const float Q, const float R) const { return 2.0f/pif*Q/sq(R); } // center velocity in 3D Poiseuille flow in channel u = 2/pi*Q/R^2, force per volume f = 4*u*rho*nu/R^2, flow rate Q = pi/8*f*R^4/(rho*nu) = pi/2*R^2*u
	float f_from_u_Poiseuille_2D(const float u, const float rho, const float nu, const float R) const { return 2.0f*u*rho*nu/sq(R); } // force per volume for 2D Poiseuille flow f = 2*u*rho*nu/R^2, u is center velocity, R is channel radius
	float f_from_u_Poiseuille_3D(const float u, const float rho, const float nu, const float R) const { return 4.0f*u*rho*nu/sq(R); } // force per volume for 3D Poiseuille flow in cylinder f = 4*u*rho*nu/R^2, u is center velocity, R is cylinder radius
	float f_from_u_rectangular_duct(const float w, const float h, const float rho, const float nu, const float u) { // force per volume f from center velocity u in rectangular channel with (-w/2<=y<=w/2, -h/2<=z<=h/2)
		const int N = 23u; // ux will become NaN for N>23 due to cosh(x) blowing up
		const double A = (cb(pi)*(double)rho*(double)nu*(double)u)/(4.0*sq((double)h));
		const double piw2h=pi*0.5*(double)w/(double)h, pi2=pi*0.5;
		double sum = 0.0;
		for(int i=2*N-1; i>=1; i-=2) { // calculate terms in reverse (smallest first) to mitigate numerical loss of significance
			const double n = (double)i;
			sum += (1.0-1.0/cosh(n*piw2h))/cb(n)*sin(n*pi2);
		}
		return (float)(A/sum); // u = A*sum
	}
	float3 u_Stokes(const float3& x, const float3& u0, const float R) const { // returns velocity at point x for laminar flow with velocity u0 around stationary sphere with radius R in origin
		const float r=1.0f/length(x), r3=cb(r), r5=sq(r)*r3, R2=sq(R);
		return u0-(1.0f/r>R ? (u0*(r+R2*r3/3.0f)+x*dot(u0, x)*(r3-R2*r5))*0.75f*R : u0);
	}
	float rho_Stokes(const float3& x, const float3& u0, const float R, const float rho, const float nu) const { // returns density at point x for laminar flow with velocity u0 around stationary sphere with radius R in origin
		const float r = length(x);
		return r>R ? rho-4.5f*rho*nu*R*dot(u0, x)/cb(r) : rho;
	}
	float f_Stokes(const float rho, const float u, const float nu, const float R, const float V) const { return 6.0f*pif*rho*nu*R*u/V; } // returns Stokes drag force density for laminar flow past stationary sphere with radius R in fluid volume V, F = 6*pi*mu*R*u, f = F/V
	float F_Stokes(const float rho, const float u, const float nu, const float R) const { return 6.0f*pif*rho*nu*R*u; } // returns Stokes drag force for laminar flow past stationary sphere with radius R, F = 6*pi*mu*R*u

	// other conversions in SI units (can be called before set_m_kg_s(...);)
	float si_Re(const float Re) const { return Re; } // Reynolds number Re = x*u/nu = [1] no unit
	float si_Re(const float si_x, const float si_u, const float si_nu) const { return si_x*si_u/si_nu; } // Reynolds number Re = x*u/nu = [1] no unit
	float si_Re(const float si_x, const float si_u, const float si_mu, const float si_rho) const { return si_x*si_u*si_rho/si_mu; } // Reynolds number Re = x*u*rho/mu = [1] no unit
	float si_We(const float si_x, const float si_u, const float si_rho, const float si_sigma) const { return si_x*sq(si_u)*si_rho/si_sigma; } // Weber number We = x*u^2*rho/sigma = [1] no unit
	float si_Fr(const float si_x, const float si_u, const float si_g) const { return si_u/sqrt(si_x*si_g); } // Froude number Fr = u/sqrt(x*g) = [1] no unit
	float si_Ca(const float si_u, const float si_mu, const float si_sigma) const { return si_u*si_mu/si_sigma; } // Capillary number Ca = u*mu/sigma = [1] no unit
	float si_Ca(const float si_u, const float si_rho, const float si_nu, const float si_sigma) const { return si_u*si_rho*si_nu/si_sigma; } // Capillary number Ca = u*rho*nu/sigma = [1] no unit
	float si_Bo(const float si_x, const float si_rho, const float si_g, const float si_sigma) const { return sq(si_x)*si_rho*si_g/si_sigma; } // Bond number Bo = x^2*rho*g/sigma = [1] no unit
	float si_Mo(const float si_rho, const float si_delta_rho, const float si_g, const float si_sigma, const float si_mu) const { return Mo(si_rho, si_delta_rho, si_g, si_sigma, si_mu); } // Morton number g*delta_rho*mu^4/(sigma^3*rho^2)
	float si_Ga(const float si_x, const float si_nu, const float si_g) { return cb(si_x)*si_g/sq(si_nu); } // Galilei number Ga = x^3*g/nu^2 = [1] no unit
	float si_Ga(const float si_x, const float si_rho, const float si_nu, const float si_f) { return cb(si_x)*si_f/(sq(si_nu)*si_rho); } // Galilei number Ga = x^3*g/nu^2 = [1] no unit
	float si_nu_from_si_mu(const float si_mu, const float si_rho) const { return si_mu/si_rho; } // kinematic shear viscosity nu = mu/rho = [m^2/s]
	float si_nu_from_si_Re(const float si_Re, const float si_x, const float si_u) const { return si_x*si_u/si_Re; } // kinematic shear viscosity nu = x*u/Re = [m^2/s]
	float si_mu_from_si_nu(const float si_nu, const float si_rho) const { return si_nu*si_rho; } // dynamic viscosity mu = nu*rho = [kg/(m*s)]
	float si_f_from_si_g(const float si_g, const float si_rho) const { return si_rho*si_g; } // force per volume f = rho*g = [kg/(m*s)^2]
	float si_g_from_si_f(const float si_f, const float si_rho) const { return si_f/si_rho; } // force per volume f = rho*g = [kg/(m*s)^2]
	float si_u_from_si_Re(const float si_Re, const float si_x, const float si_nu) const { return si_Re*si_nu/si_x; } // velocity u = Re*nu/x, Reynolds number Re = x*u/nu = [1] no unit
	float si_u_from_si_Re(const float si_Re, const float si_x, const float si_mu, const float si_rho) const { return si_Re*si_mu/(si_x*si_rho); } // velocity u = Re*nu/x, Reynolds number Re = x*u*rho/mu = [1] no unit
	float si_u_from_si_We(const float si_We, const float si_x, const float si_sigma, const float si_rho) const { return sqrt(si_We*si_sigma/(si_x*si_rho)); } // velocity u = sqrt(We*sigma/(x*rho)) = [m/s], Weber number We = x*u^2*rho/sigma = [1] no unit
	float si_u_from_si_Fr(const float si_Fr, const float si_x, const float si_g) const { return si_Fr*sqrt(si_x*si_g); } // velocity u = Fr*sqrt(x*g) = [m/s], Froude number Fr = u/sqrt(x*g) = [1] no unit
	float si_u_from_si_h(const float si_h, const float si_g) const { return sqrt(2.0f*si_h*si_g); } // free fall drop height h = 1/2*g*t^2 -> t = sqrt(2*h/g), u = a*t = g*sqrt(2*h/g) = sqrt(2*h*g) = [m/s]
	float si_u_Poiseuille_2D(const float si_Q, const float si_R) const { return 0.75f*si_Q/si_R; } // center velocity in 2D Poiseuille flow in channel u = 3/4*Q/R, force per volume f = 2*u*rho*nu/R^2, flow rate Q = 2/3*f*R^3/(rho*nu) = 4/3*u*R
	float si_u_Poiseuille_3D(const float si_Q, const float si_R) const { return 2.0f/pif*si_Q/sq(si_R); } // center velocity in 3D Poiseuille flow in channel u = 2/pi*Q/R^2, force per volume f = 4*u*rho*nu/R^2, flow rate Q = pi/8*f*R^4/(rho*nu) = pi/2*R^2*u
};
extern Units units; // declared in lbm.cpp