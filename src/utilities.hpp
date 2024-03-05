#pragma once

#define UTILITIES_REGEX
#define UTILITIES_FILE
#define UTILITIES_PNG
#define UTILITIES_CONSOLE_COLOR
#define UTILITIES_CONSOLE_INPUT
#define UTILITIES_CONSOLE_DITHER_LOOKUP
#define CONSOLE_WIDTH 79
//#define UTILITIES_NO_CPP17

#pragma warning(disable:26451)
#pragma warning(disable:6386)
#include <cmath>
#include <vector>
#ifdef UTILITIES_REGEX
#include <regex> // contains <string>, <vector>, <algorithm> and others
#else // UTILITIES_REGEX
#include <string>
#endif // UTILITIES_REGEX
#include <iostream>
#include <thread> // contains <chrono>
#include <functional> // for parallel_for(...)
#undef min
#undef max
using std::string;
using std::vector;
using std::thread;
typedef unsigned char uchar;
typedef unsigned short ushort;
typedef unsigned int uint;
typedef int64_t slong;
typedef uint64_t ulong;
#define pif 3.1415927f
#define pi 3.141592653589793
#define min_char ((char)-128)
#define max_char ((char)127)
#define max_uchar ((uchar)255)
#define min_short ((short)-32768)
#define max_short ((short)32767)
#define max_ushort ((ushort)65535)
#define min_int -2147483648
#define max_int 2147483647
#define max_uint 4294967295u
#define min_slong -9223372036854775808ll
#define max_slong 9223372036854775807ll
#define max_ulong 18446744073709551615ull
#define min_float 1.401298464E-45f
#define max_float 3.402823466E38f
#define epsilon_float 1.192092896E-7f
#define inf_float as_float(0x7F800000)
#define nan_float as_float(0xFFFFFFFF)
#define min_double 4.9406564584124654E-324
#define max_double 1.7976931348623158E308
#define epsilon_double 2.2204460492503131E-16
#define inf_double as_double(0x7FF0000000000000)
#define nan_double as_double(0xFFFFFFFFFFFFFFFF)

inline void parallel_for(const uint N, const uint threads, std::function<void(uint, uint)> lambda) { // usage: parallel_for(N, threads, [&](uint n, uint t) { ... });
	vector<thread> thread_array(threads);
	for(uint t=0u; t<threads; t++) thread_array[t] = thread([=]() {
		for(ulong n=(ulong)N*(ulong)t/(ulong)threads; n<(ulong)N*(ulong)(t+1u)/(ulong)threads; n++) lambda((uint)n, t);
	});
	for(uint t=0u; t<threads; t++) thread_array[t].join();
}
inline void parallel_for(const uint N, const uint threads, std::function<void(uint)> lambda) { // usage: parallel_for(N, threads, [&](uint n) { ... });
	vector<thread> thread_array(threads);
	for(uint t=0u; t<threads; t++) thread_array[t] = thread([=]() {
		for(ulong n=(ulong)N*(ulong)t/(ulong)threads; n<(ulong)N*(ulong)(t+1u)/(ulong)threads; n++) lambda((uint)n);
	});
	for(uint t=0u; t<threads; t++) thread_array[t].join();
}
inline void parallel_for(const uint N, std::function<void(uint)> lambda) { // usage: parallel_for(N, [&](uint n) { ... });
	parallel_for(N, (uint)thread::hardware_concurrency(), lambda);
}
inline void parallel_for(const ulong N, const uint threads, std::function<void(ulong, uint)> lambda) { // usage: parallel_for(N, threads, [&](ulong n, uint t) { ... });
	vector<thread> thread_array(threads);
	for(uint t=0u; t<threads; t++) thread_array[t] = thread([=]() {
		for(ulong n=N*(ulong)t/(ulong)threads; n<N*(ulong)(t+1u)/(ulong)threads; n++) lambda(n, t);
	});
	for(uint t=0u; t<threads; t++) thread_array[t].join();
}
inline void parallel_for(const ulong N, const uint threads, std::function<void(ulong)> lambda) { // usage: parallel_for(N, threads, [&](ulong n) { ... });
	vector<thread> thread_array(threads);
	for(uint t=0u; t<threads; t++) thread_array[t] = thread([=]() {
		for(ulong n=N*(ulong)t/(ulong)threads; n<N*(ulong)(t+1u)/(ulong)threads; n++) lambda(n);
	});
	for(uint t=0u; t<threads; t++) thread_array[t].join();
}
inline void parallel_for(const ulong N, std::function<void(ulong)> lambda) { // usage: parallel_for(N, [&](ulong n) { ... });
	parallel_for(N, (uint)thread::hardware_concurrency(), lambda);
}

class Clock {
private:
	typedef std::chrono::high_resolution_clock clock;
	std::chrono::time_point<clock> t;
public:
	inline Clock() { start(); }
	inline void start() { t = clock::now(); }
	inline double stop() const { return std::chrono::duration_cast<std::chrono::duration<double>>(clock::now()-t).count(); }
};
inline void sleep(const double t) {
	if(t>0.0) std::this_thread::sleep_for(std::chrono::milliseconds((int)(1E3*t+0.5)));
}

inline float as_float(const uint x) {
	return *(float*)&x;
}
inline uint as_uint(const float x) {
	return *(uint*)&x;
}
inline double as_double(const ulong x) {
	return *(double*)&x;
}
inline ulong as_ulong(const double x) {
	return *(ulong*)&x;
}

inline float half_to_float(const ushort x) { // IEEE-754 16-bit floating-point format (without infinity): 1-5-10, exp-15, +-131008.0, +-6.1035156E-5, +-5.9604645E-8, 3.311 digits
	const uint e = (x&0x7C00)>>10; // exponent
	const uint m = (x&0x03FF)<<13; // mantissa
	const uint v = as_uint((float)m)>>23; // evil log2 bit hack to count leading zeros in denormalized format
	return as_float((x&0x8000)<<16 | (e!=0)*((e+112)<<23|m) | ((e==0)&(m!=0))*((v-37)<<23|((m<<(150-v))&0x007FE000))); // sign : normalized : denormalized
}
inline ushort float_to_half(const float x) { // IEEE-754 16-bit floating-point format (without infinity): 1-5-10, exp-15, +-131008.0, +-6.1035156E-5, +-5.9604645E-8, 3.311 digits
	const uint b = as_uint(x)+0x00001000; // round-to-nearest-even: add last bit after truncated mantissa
	const uint e = (b&0x7F800000)>>23; // exponent
	const uint m = b&0x007FFFFF; // mantissa; in line below: 0x007FF000 = 0x00800000-0x00001000 = decimal indicator flag - initial rounding
	return (b&0x80000000)>>16 | (e>112)*((((e-112)<<10)&0x7C00)|m>>13) | ((e<113)&(e>101))*((((0x007FF000+m)>>(125-e))+1)>>1) | (e>143)*0x7FFF; // sign : normalized : denormalized : saturate
}
inline float half_to_float_custom(const ushort x) { // custom 16-bit floating-point format, 1-4-11, exp-15, +-1.99951168, +-6.10351562E-5, +-2.98023224E-8, 3.612 digits
	const uint e = (x&0x7800)>>11; // exponent
	const uint m = (x&0x07FF)<<12; // mantissa
	const uint v = as_uint((float)m)>>23; // evil log2 bit hack to count leading zeros in denormalized format
	return as_float((x&0x8000)<<16 | (e!=0)*((e+112)<<23|m) | ((e==0)&(m!=0))*((v-37)<<23|((m<<(150-v))&0x007FF000))); // sign : normalized : denormalized
}
inline ushort float_to_half_custom(const float x) { // custom 16-bit floating-point format, 1-4-11, exp-15, +-1.99951168, +-6.10351562E-5, +-2.98023224E-8, 3.612 digits
	const uint b = as_uint(x)+0x00000800; // round-to-nearest-even: add last bit after truncated mantissa
	const uint e = (b&0x7F800000)>>23; // exponent
	const uint m = b&0x007FFFFF; // mantissa; in line below: 0x007FF800 = 0x00800000-0x00000800 = decimal indicator flag - initial rounding
	return (b&0x80000000)>>16 | (e>112)*((((e-112)<<11)&0x7800)|m>>12) | ((e<113)&(e>100))*((((0x007FF800+m)>>(124-e))+1)>>1) | (e>127)*0x7FFF; // sign : normalized : denormalized : saturate
}

inline float sq(const float x) {
	return x*x;
}
inline float cb(const float x) {
	return x*x*x;
}
inline float pow(const float x, const uint n) {
	float r = 1.0f;
	for(uint i=0u; i<n; i++) {
		r *= x;
	}
	return r;
}
inline float sign(const float x) {
	return x>=0.0f ? 1.0f : -1.0f;
}
inline float clamp(const float x, const float a, const float b) {
	return fmin(fmax(x, a), b);
}
inline float rsqrt(const float x) {
	return 1.0f/sqrt(x);
}
inline float ln(const float x) {
	return log(x); // natural logarithm
}
inline int log2_fast(const float x) {
	return (as_uint(x)>>23)-127;
}
inline float degrees(const float radians) {
	return (180.0f/pif)*radians;
}
inline float radians(const float degrees) {
	return (pif/180.0f)*degrees;
}
inline float find_zero(float f(float), float min=-1.0f, float max=1.0f, float offset=0.0f) { // find zero of function f(x) in [min, max] by nested intervals
	const uint n = log2_fast((max-min)*1.67772162E7f); // evil log2 hack: log2(x)=(as_uint(x)>>23)-127
	float m = 0.5f*(min+max);
	const float s = 2.0f*(float)(f(min)>offset)-1.0f;
	offset *= s;
	for(uint i=0u; i<=n; i++) {
		if(s*f(m)>offset) min = m;
		else /**********/ max = m;
		m = 0.5f*(min+max);
	}
	return m;
}
inline float integrate(float f(float), const float min, const float max, const uint N=1000000u) { // do numeric integration of function f(x)dx in [min, max] by trapezoidal rule
	const float dx = (max-min)/(float)N;
	double s = 0.0;
	for(uint i=1u; i+1<N; i++) {
		s += (double)f(min+(float)i*dx);
	}
	return (float)(0.5*((double)f(min)+(double)f(max))+s)*dx;
}
inline float derivative(float f(float), const float x) { // compute central first derivative of f(x) at x
	const float dx = fmax(fabs(x), 1.0f)*2.0f*epsilon_float;
	return (float)((double)f(x+dx)-(double)f(x-dx))*0.5f/dx;
}
inline float second_derivative(float f(float), const float x) { // compute central second derivative of f(x) at x
	const float dx = fmax(fabs(x), 1.0f)*2.0f*epsilon_float;
	return (float)((double)f(x+dx)-2.0*(double)f(x)+(double)f(x-dx))/sq(dx);
}
inline bool converged(const float x_penultimate, const float x_last, const float x, const float threshold=1E-5f) {
	const float slope = fabs(x-x_last);
	const float curvature = fabs(x-2.0f*x_last+x_penultimate);
	const float t = (x_penultimate+x_last+x)/3.0f*threshold; // small threshold relative to x
	return slope<t&&curvature<t; // convergence is reached if both slope and curvature are close to 0
}
inline float fmin(const uint n, const float* const x) {
	float r = x[0];
	for(uint i=1u; i<n; i++) {
		r = fmin(r, x[i]);
	}
	return r;
}
inline float fmax(const uint n, const float* const x) {
	float r = x[0];
	for(uint i=1u; i<n; i++) {
		r = fmax(r, x[i]);
	}
	return r;
}
inline float average(const uint n, const float* const x) {
	double a = 0.0;
	for(uint i=0u; i<n; i++) {
		a += (double)x[i];
	}
	return (float)(a/(double)n);
}
inline float standard_deviation(const uint n, const float* const x) {
	const double a = (double)average(n, x);
	double s = 0.0;
	for(uint i=0u; i<n; i++) {
		const double b = (double)x[i]-a;
		s += b*b;
	}
	return (float)sqrt(s/(double)n);
}
inline float random(uint& seed, const float x=1.0f) {
	seed = (1103515245u*seed+12345u)%2147483648u; // standard C99 LCG
	return x*(float)seed*4.65661287E-10f;
}
inline float random_symmetric(uint& seed, const float x=1.0f) {
	return 2.0f*x*(random(seed, 1.0f)-0.5f);
}
inline void lu_solve(float* M, float* x, float* b, const int N) { // solves system of N linear equations M*x=b
	for(int i=0; i<N; i++) { // decompose M in M=L*U
		for(int j=i+1; j<N; j++) {
			M[N*j+i] /= M[N*i+i];
			for(int k=i+1; k<N; k++) M[N*j+k] -= M[N*j+i]*M[N*i+k];
		}
	}
	for(int i=0; i<N; i++) { // find solution of L*y=b
		x[i] = b[i];
		for(int k=0; k<i; k++) x[i] -= M[N*i+k]*x[k];
	}
	for(int i=N-1; i>=0; i--) { // find solution of U*x=y
		for(int k=i+1; k<N; k++) x[i] -= M[N*i+k]*x[k];
		x[i] /= M[N*i+i];
	}
}

inline double sq(const double x) {
	return x*x;
}
inline double cb(const double x) {
	return x*x*x;
}
inline double pow(const double x, const uint n) {
	double r = 1.0;
	for(uint i=0u; i<n; i++) {
		r *= x;
	}
	return r;
}
inline double sign(const double x) {
	return x>=0.0 ? 1.0 : -1.0;
}
inline double clamp(const double x, const double a, const double b) {
	return fmin(fmax(x, a), b);
}
inline double rsqrt(const double x) {
	return 1.0/sqrt(x);
}
inline double ln(const double x) {
	return log(x); // natural logarithm
}
inline int log2_fast(const double x) {
	return (as_ulong(x)>>52)-1023;
}
inline double degrees(const double radians) {
	return (180.0/pi)*radians;
}
inline double radians(const double degrees) {
	return (pi/180.0)*degrees;
}
inline double find_zero(double f(double), double min=-1.0, double max=1.0, double offset=0.0) { // find zero of function f(x) in [min, max] by nested intervals
	const uint n = log2_fast((max-min)*9.007199254740992E15); // evil log2 hack: log2(x)=(as_ulong(x)>>52)-1023
	double m = 0.5*(min+max);
	const double s = 2.0*(double)(f(min)>offset)-1.0;
	offset *= s;
	for(uint i=0u; i<=n; i++) {
		if(s*f(m)>offset) min = m;
		else /**********/ max = m;
		m = 0.5*(min+max);
	}
	return m;
}
inline double integrate(double f(double), const double min, const double max, const uint N=1000000u) { // do numeric integration of function f(x)dx in [min, max] by trapezoidal rule
	const double dx = (max-min)/(double)N;
	double s = 0.0;
	for(uint i=1u; i+1<N; i++) {
		s += f(min+(double)i*dx);
	}
	return (0.5*(f(min)+f(max))+s)*dx;
}
inline double derivative(double f(double), const double x) { // compute central first derivative of f(x) at x
	const double dx = fmax(fabs(x), 1.0)*2.0*epsilon_double;
	return (f(x+dx)-f(x-dx))*0.5/dx;
}
inline double second_derivative(double f(double), const double x) { // compute central second derivative of f(x) at x
	const double dx = fmax(fabs(x), 1.0)*2.0*epsilon_double;
	return (f(x+dx)-2.0*f(x)+f(x-dx))/sq(dx);
}
inline bool converged(const double x_penultimate, const double x_last, const double x, const double threshold=1E-5) {
	const double slope = fabs(x-x_last);
	const double curvature = fabs(x-2.0*x_last+x_penultimate);
	const double t = (x_penultimate+x_last+x)/3.0*threshold; // small threshold relative to x
	return slope<t&&curvature<t; // convergence is reached if both slope and curvature are close to 0
}
inline double fmin(const uint n, const double* const x) {
	double r = x[0];
	for(uint i=1u; i<n; i++) {
		r = fmin(r, x[i]);
	}
	return r;
}
inline double fmax(const uint n, const double* const x) {
	double r = x[0];
	for(uint i=1u; i<n; i++) {
		r = fmax(r, x[i]);
	}
	return r;
}
inline double average(const uint n, const double* const x) {
	double a = 0.0;
	for(uint i=0u; i<n; i++) {
		a += x[i];
	}
	return (float)(a/(double)n);
}
inline double standard_deviation(const uint n, const double* const x) {
	const double a = average(n, x);
	double s = 0.0;
	for(uint i=0u; i<n; i++) {
		s += sq(x[i]-a);
	}
	return sqrt(s/(double)n);
}

inline int sq(const int x) {
	return x*x;
}
inline int cb(const int x) {
	return x*x*x;
}
inline int pow(const int x, const uint n) {
	int r = 1;
	for(uint i=0u; i<n; i++) {
		r *= x;
	}
	return r;
}
inline int sign(const int x) {
	return 1-2*(x>>31&1);
}
inline int min(const int x, const int y) {
	return x<y?x:y;
}
inline int max(const int x, const int y) {
	return x>y?x:y;
}
inline int clamp(const int x, const int a, const int b) {
	return min(max(x, a), b);
}

inline uint sq(const uint x) {
	return x*x;
}
inline uint cb(const uint x) {
	return x*x*x;
}
inline uint pow(const uint x, const uint n) {
	uint r = 1u;
	for(uint i=0u; i<n; i++) {
		r *= x;
	}
	return r;
}
inline uint min(const uint x, const uint y) {
	return x<y?x:y;
}
inline uint max(const uint x, const uint y) {
	return x>y?x:y;
}
inline uint clamp(const uint x, const uint a, const uint b) {
	return min(max(x, a), b);
}
inline uint gcd(uint x, uint y) { // greatest common divisor
	if(x*y==0u) return 0u;
	uint t;
	while(y!=0u) {
		t = y;
		y = x%y;
		x = t;
	}
	return x;
}
inline uint lcm(const uint x, const uint y) { // least common multiple
	return x*y==0u ? 0u : x*y/gcd(x, y);
}
inline uint log2_fast(const uint x) {
	return (uint)((as_uint((float)x)>>23)-127);
}

inline slong sq(const slong x) {
	return x*x;
}
inline slong cb(const slong x) {
	return x*x*x;
}
inline slong pow(const slong x, const uint n) {
	slong r = 1ll;
	for(uint i=0u; i<n; i++) {
		r *= x;
	}
	return r;
}
inline slong sign(const slong x) {
	return 1ll-2ll*(x>>63&1ll);
}
inline slong min(const slong x, const slong y) {
	return x<y?x:y;
}
inline slong max(const slong x, const slong y) {
	return x>y?x:y;
}
inline slong clamp(const slong x, const slong a, const slong b) {
	return min(max(x, a), b);
}

inline ulong sq(const ulong x) {
	return x*x;
}
inline ulong cb(const ulong x) {
	return x*x*x;
}
inline ulong pow(const ulong x, const uint n) {
	ulong r = 1ull;
	for(uint i=0u; i<n; i++) {
		r *= x;
	}
	return r;
}
inline ulong min(const ulong x, const ulong y) {
	return x<y?x:y;
}
inline ulong max(const ulong x, const ulong y) {
	return x>y?x:y;
}
inline ulong clamp(const ulong x, const ulong a, const ulong b) {
	return min(max(x, a), b);
}
inline ulong gcd(ulong x, ulong y) { // greatest common divisor
	if(x*y==0ull) return 0ull;
	ulong t;
	while(y!=0ull) {
		t = y;
		y = x%y;
		x = t;
	}
	return x;
}
inline ulong lcm(const ulong x, const ulong y) { // least common multiple
	return x*y==0ull ? 0ull : x*y/gcd(x, y);
}
inline uint log2_fast(const ulong x) {
	return (uint)((as_ulong((double)x)>>52)-1023);
}

inline int to_int(const float x) {
	return (int)(x+0.5f-(float)(x<0.0f));
}
inline int to_int(const double x) {
	return (int)(x+0.5-(double)(x<0.0));
}
inline uint to_uint(const float x) {
	return (uint)fmax(x+0.5f, 0.5f);
}
inline uint to_uint(const double x) {
	return (uint)fmax(x+0.5, 0.5);
}
inline slong to_slong(const float x) {
	return (slong)(x+0.5f);
}
inline slong to_slong(const double x) {
	return (slong)(x+0.5);
}
inline ulong to_ulong(const float x) {
	return (ulong)fmax(x+0.5f, 0.5f);
}
inline ulong to_ulong(const double x) {
	return (ulong)fmax(x+0.5, 0.5);
}

inline float reverse_bytes(const float x) { // reverse the order of bytes
	float r;
	char* cx = (char*)&x;
	char* cr = (char*)&r;
	cr[0] = cx[3];
	cr[1] = cx[2];
	cr[2] = cx[1];
	cr[3] = cx[0];
	return r;
}
inline double reverse_bytes(const double x) { // reverse the order of bytes
	double r;
	char* cx = (char*)&x;
	char* cr = (char*)&r;
	cr[0] = cx[7];
	cr[1] = cx[6];
	cr[2] = cx[5];
	cr[3] = cx[4];
	cr[4] = cx[3];
	cr[5] = cx[2];
	cr[6] = cx[1];
	cr[7] = cx[0];
	return r;
}
inline char reverse_bytes(const char x) { // reverse the order of bytes
	return x;
}
inline uchar reverse_bytes(const uchar x) { // reverse the order of bytes
	return x;
}
inline short reverse_bytes(const short x) { // reverse the order of bytes
	short r;
	char* cx = (char*)&x;
	char* cr = (char*)&r;
	cr[0] = cx[1];
	cr[1] = cx[0];
	return r;
}
inline ushort reverse_bytes(const ushort x) { // reverse the order of bytes
	ushort r;
	char* cx = (char*)&x;
	char* cr = (char*)&r;
	cr[0] = cx[1];
	cr[1] = cx[0];
	return r;
}
inline int reverse_bytes(const int x) { // reverse the order of bytes
	int r;
	char* cx = (char*)&x;
	char* cr = (char*)&r;
	cr[0] = cx[3];
	cr[1] = cx[2];
	cr[2] = cx[1];
	cr[3] = cx[0];
	return r;
}
inline uint reverse_bytes(const uint x) { // reverse the order of bytes
	uint r;
	char* cx = (char*)&x;
	char* cr = (char*)&r;
	cr[0] = cx[3];
	cr[1] = cx[2];
	cr[2] = cx[1];
	cr[3] = cx[0];
	return r;
}
inline slong reverse_bytes(const slong x) { // reverse the order of bytes
	slong r;
	char* cx = (char*)&x;
	char* cr = (char*)&r;
	cr[0] = cx[7];
	cr[1] = cx[6];
	cr[2] = cx[5];
	cr[3] = cx[4];
	cr[4] = cx[3];
	cr[5] = cx[2];
	cr[6] = cx[1];
	cr[7] = cx[0];
	return r;
}
inline ulong reverse_bytes(const ulong x) { // reverse the order of bytes
	ulong r;
	char* cx = (char*)&x;
	char* cr = (char*)&r;
	cr[0] = cx[7];
	cr[1] = cx[6];
	cr[2] = cx[5];
	cr[3] = cx[4];
	cr[4] = cx[3];
	cr[5] = cx[2];
	cr[6] = cx[1];
	cr[7] = cx[0];
	return r;
}

struct int3 { // 3D vector type
	int x, y, z;
	int3(const int x, const int y, const int z) {
		this->x = x;
		this->y = y;
		this->z = z;
	}
	int3(const uint x, const uint y, const uint z) {
		this->x = (int)x;
		this->y = (int)y;
		this->z = (int)z;
	}
	int3(const float x, const float y, const float z) {
		this->x = to_int(x);
		this->y = to_int(y);
		this->z = to_int(z);
	}
	int3(const double x, const double y, const double z) {
		this->x = to_int(x);
		this->y = to_int(y);
		this->z = to_int(z);
	}
	int3(const int x) {
		this->x = this->y = this->z = x;
	}
	int3(const uint x) {
		this->x = this->y = this->z = (int)x;
	}
	int3(const float x) {
		this->x = this->y = this->z = to_int(x);
	}
	int3(const double x) {
		this->x = this->y = this->z = to_int(x);
	}
	int3() = default;
	int3(const int3& v) = default; // copy constructor (elementwise copy)
	int3(int3&& v) = default; // move constructor
	~int3() = default;
	inline int3& operator=(const int3& v) = default; // copy assignment (elementwise copy)
	inline int3& operator=(int3&& v) = default; // move assignment
	inline int3& operator+=(const int3& v) { // elementwise addition
		this->x += v.x;
		this->y += v.y;
		this->z += v.z;
		return *this;
	}
	inline int3& operator-=(const int3& v) { // elementwise subtraction
		this->x -= v.x;
		this->y -= v.y;
		this->z -= v.z;
		return *this;
	}
	inline int3& operator+=(const int x) { // elementwise addition
		this->x += x;
		this->y += x;
		this->z += x;
		return *this;
	}
	inline int3& operator-=(const int x) { // elementwise subtraction
		this->x -= x;
		this->y -= x;
		this->z -= x;
		return *this;
	}
	inline int3& operator*=(const int x) { // elementwise multiplication
		this->x *= x;
		this->y *= x;
		this->z *= x;
		return *this;
	}
	inline int3& operator/=(const int x) { // elementwise division
		this->x /= x;
		this->y /= x;
		this->z /= x;
		return *this;
	}
	inline int3& operator+() {
		return *this;
	}
	inline const int3& operator+() const {
		return *this;
	}
	inline int3 operator-() const { // elementwise negation
		return int3(-this->x, -this->y, -this->z);
	}
	inline int3 operator+(const int3& v) const { // elementwise addition
		return int3(this->x+v.x, this->y+v.y, this->z+v.z);
	}
	inline int3 operator-(const int3& v) const { // elementwise subtraction
		return int3(this->x-v.x, this->y-v.y, this->z-v.z);
	}
	inline int operator*(const int3& v) const { // dot product
		return this->x*v.x+this->y*v.y+this->z*v.z;
	}
	inline bool operator==(const int3& v) const { // check for equality
		return this->x==v.x&&this->y==v.y&&this->z==v.z;
	}
	inline bool operator!=(const int3& v) const { // check for inequality
		return this->x!=v.x||this->y!=v.y||this->z!=v.z;
	}
	inline bool operator>(const int3& v) const { // compare length
		return sq(this->x)+sq(this->y)+sq(this->z)>sq(v.x)+sq(v.y)+sq(v.z);
	}
	inline bool operator<(const int3& v) const { // compare length
		return sq(this->x)+sq(this->y)+sq(this->z)<sq(v.x)+sq(v.y)+sq(v.z);
	}
	inline bool operator>=(const int3& v) const { // compare length
		return sq(this->x)+sq(this->y)+sq(this->z)>=sq(v.x)+sq(v.y)+sq(v.z);
	}
	inline bool operator<=(const int3& v) const { // compare length
		return sq(this->x)+sq(this->y)+sq(this->z)<=sq(v.x)+sq(v.y)+sq(v.z);
	}
};
inline int3 operator+(const int3& v, const int x) { // elementwise addition
	return int3(v.x+x, v.y+x, v.z+x);
}
inline int3 operator+(const int x, const int3& v) { // elementwise addition
	return v+x;
}
inline int3 operator-(const int3& v, const int x) { // elementwise subtraction
	return v+(-x);
}
inline int3 operator-(const int x, const int3& v) { // elementwise subtraction
	return int3(x-v.x, x-v.y, x-v.z);
}
inline int3 operator*(const int3& v, const int x) { // elementwise multiplication
	return int3(v.x*x, v.y*x, v.z*x);
}
inline int3 operator*(const int x, const int3& v) { // elementwise multiplication
	return v*x;
}
inline int3 operator/(const int3& v, const int x) { // elementwise division
	return int3(v.x/x, v.y/x, v.z/x);
}

struct uint3 { // 3D vector type
	uint x, y, z;
	uint3(const uint x, const uint y, const uint z) {
		this->x = x;
		this->y = y;
		this->z = z;
	}
	uint3(const int x, const int y, const int z) {
		this->x = (uint)x;
		this->y = (uint)y;
		this->z = (uint)z;
	}
	uint3(const float x, const float y, const float z) {
		this->x = to_uint(x);
		this->y = to_uint(y);
		this->z = to_uint(z);
	}
	uint3(const double x, const double y, const double z) {
		this->x = to_uint(x);
		this->y = to_uint(y);
		this->z = to_uint(z);
	}
	uint3(const uint x) {
		this->x = this->y = this->z = x;
	}
	uint3(const int x) {
		this->x = this->y = this->z = (uint)x;
	}
	uint3(const float x) {
		this->x = this->y = this->z = to_uint(x);
	}
	uint3(const double x) {
		this->x = this->y = this->z = to_uint(x);
	}
	uint3() = default;
	uint3(const uint3& v) = default; // copy constructor (elementwise copy)
	uint3(uint3&& v) = default; // move constructor
	~uint3() = default;
	inline uint3& operator=(const uint3& v) = default; // copy assignment (elementwise copy)
	inline uint3& operator=(uint3&& v) = default; // move assignment
	inline uint3& operator+=(const uint3& v) { // elementwise addition
		this->x += v.x;
		this->y += v.y;
		this->z += v.z;
		return *this;
	}
	inline uint3& operator-=(const uint3& v) { // elementwise subtraction
		this->x -= v.x;
		this->y -= v.y;
		this->z -= v.z;
		return *this;
	}
	inline uint3& operator+=(const uint x) { // elementwise addition
		this->x += x;
		this->y += x;
		this->z += x;
		return *this;
	}
	inline uint3& operator-=(const uint x) { // elementwise subtraction
		this->x -= x;
		this->y -= x;
		this->z -= x;
		return *this;
	}
	inline uint3& operator*=(const uint x) { // elementwise multiplication
		this->x *= x;
		this->y *= x;
		this->z *= x;
		return *this;
	}
	inline uint3& operator/=(const uint x) { // elementwise division
		this->x /= x;
		this->y /= x;
		this->z /= x;
		return *this;
	}
	inline uint3& operator+() {
		return *this;
	}
	inline const uint3& operator+() const {
		return *this;
	}
	inline uint3 operator+(const uint3& v) const { // elementwise addition
		return uint3(this->x+v.x, this->y+v.y, this->z+v.z);
	}
	inline uint3 operator-(const uint3& v) const { // elementwise subtraction
		return uint3(this->x-v.x, this->y-v.y, this->z-v.z);
	}
	inline uint operator*(const uint3& v) const { // dot product
		return this->x*v.x+this->y*v.y+this->z*v.z;
	}
	inline bool operator==(const uint3& v) const { // check for equality
		return this->x==v.x&&this->y==v.y&&this->z==v.z;
	}
	inline bool operator!=(const uint3& v) const { // check for inequality
		return this->x!=v.x||this->y!=v.y||this->z!=v.z;
	}
	inline bool operator>(const uint3& v) const { // compare length
		return sq(this->x)+sq(this->y)+sq(this->z)>sq(v.x)+sq(v.y)+sq(v.z);
	}
	inline bool operator<(const uint3& v) const { // compare length
		return sq(this->x)+sq(this->y)+sq(this->z)<sq(v.x)+sq(v.y)+sq(v.z);
	}
	inline bool operator>=(const uint3& v) const { // compare length
		return sq(this->x)+sq(this->y)+sq(this->z)>=sq(v.x)+sq(v.y)+sq(v.z);
	}
	inline bool operator<=(const uint3& v) const { // compare length
		return sq(this->x)+sq(this->y)+sq(this->z)<=sq(v.x)+sq(v.y)+sq(v.z);
	}
};
inline uint3 operator+(const uint3& v, const uint x) { // elementwise addition
	return uint3(v.x+x, v.y+x, v.z+x);
}
inline uint3 operator+(const uint x, const uint3& v) { // elementwise addition
	return v+x;
}
inline uint3 operator-(const uint3& v, const uint x) { // elementwise subtraction
	return uint3(v.x-x, v.y-x, v.z-x);
}
inline uint3 operator-(const uint x, const uint3& v) { // elementwise subtraction
	return uint3(x-v.x, x-v.y, x-v.z);
}
inline uint3 operator*(const uint3& v, const uint x) { // elementwise multiplication
	return uint3(v.x*x, v.y*x, v.z*x);
}
inline uint3 operator*(const uint x, const uint3& v) { // elementwise multiplication
	return v*x;
}
inline uint3 operator/(const uint3& v, const uint x) { // elementwise division
	return uint3(v.x/x, v.y/x, v.z/x);
}

struct float3x3; // forward-declare float3x3
struct float3 { // 3D vector type
	float x, y, z;
	float3(const float x, const float y, const float z) {
		this->x = x;
		this->y = y;
		this->z = z;
	}
	float3(const double x, const double y, const double z) {
		this->x = (float)x;
		this->y = (float)y;
		this->z = (float)z;
	}
	float3(const int x, const int y, const int z) {
		this->x = (float)x;
		this->y = (float)y;
		this->z = (float)z;
	}
	float3(const uint x, const uint y, const uint z) {
		this->x = (float)x;
		this->y = (float)y;
		this->z = (float)z;
	}
	float3(const float x) {
		this->x = this->y = this->z = x;
	}
	float3(const double x) {
		this->x = this->y = this->z = (float)x;
	}
	float3(const int x) {
		this->x = this->y = this->z = (float)x;
	}
	float3(const uint x) {
		this->x = this->y = this->z = (float)x;
	}
	float3(const float3x3& m); // forward-declare float3x3 constructor
	float3() = default;
	float3(const float3& v) = default; // copy constructor (elementwise copy)
	float3(float3&& v) = default; // move constructor
	~float3() = default;
	inline float3& operator=(const float3& v) = default; // copy assignment (elementwise copy)
	inline float3& operator=(float3&& v) = default; // move assignment
	inline float3& operator=(const float3x3& m); // forward-declare float3x3 copy
	inline float3& operator+=(const float3& v) { // elementwise addition
		this->x += v.x;
		this->y += v.y;
		this->z += v.z;
		return *this;
	}
	inline float3& operator-=(const float3& v) { // elementwise subtraction
		this->x -= v.x;
		this->y -= v.y;
		this->z -= v.z;
		return *this;
	}
	inline float3& operator+=(const float x) { // elementwise addition
		this->x += x;
		this->y += x;
		this->z += x;
		return *this;
	}
	inline float3& operator-=(const float x) { // elementwise subtraction
		this->x -= x;
		this->y -= x;
		this->z -= x;
		return *this;
	}
	inline float3& operator*=(const float x) { // elementwise multiplication
		this->x *= x;
		this->y *= x;
		this->z *= x;
		return *this;
	}
	inline float3& operator/=(const float x) { // elementwise division
		this->x /= x;
		this->y /= x;
		this->z /= x;
		return *this;
	}
	inline float3& operator+() {
		return *this;
	}
	inline const float3& operator+() const {
		return *this;
	}
	inline float3 operator-() const { // elementwise negation
		return float3(-this->x, -this->y, -this->z);
	}
	inline float3 operator+(const float3& v) const { // elementwise addition
		return float3(this->x+v.x, this->y+v.y, this->z+v.z);
	}
	inline float3 operator-(const float3& v) const { // elementwise subtraction
		return float3(this->x-v.x, this->y-v.y, this->z-v.z);
	}
	inline float operator*(const float3& v) const { // dot product
		return this->x*v.x+this->y*v.y+this->z*v.z;
	}
	inline bool operator==(const float3& v) const { // check for equality
		return this->x==v.x&&this->y==v.y&&this->z==v.z;
	}
	inline bool operator!=(const float3& v) const { // check for inequality
		return this->x!=v.x||this->y!=v.y||this->z!=v.z;
	}
	inline bool operator>(const float3& v) const { // compare length
		return sq(this->x)+sq(this->y)+sq(this->z)>sq(v.x)+sq(v.y)+sq(v.z);
	}
	inline bool operator<(const float3& v) const { // compare length
		return sq(this->x)+sq(this->y)+sq(this->z)<sq(v.x)+sq(v.y)+sq(v.z);
	}
	inline bool operator>=(const float3& v) const { // compare length
		return sq(this->x)+sq(this->y)+sq(this->z)>=sq(v.x)+sq(v.y)+sq(v.z);
	}
	inline bool operator<=(const float3& v) const { // compare length
		return sq(this->x)+sq(this->y)+sq(this->z)<=sq(v.x)+sq(v.y)+sq(v.z);
	}
};
inline float3 operator+(const float3& v, const float x) { // elementwise addition
	return float3(v.x+x, v.y+x, v.z+x);
}
inline float3 operator+(const float x, const float3& v) { // elementwise addition
	return v+x;
}
inline float3 operator-(const float3& v, const float x) { // elementwise subtraction
	return v+(-x);
}
inline float3 operator-(const float x, const float3& v) { // elementwise subtraction
	return float3(x-v.x, x-v.y, x-v.z);
}
inline float3 operator*(const float3& v, const float x) { // elementwise multiplication
	return float3(v.x*x, v.y*x, v.z*x);
}
inline float3 operator*(const float x, const float3& v) { // elementwise multiplication
	return v*x;
}
inline float3 operator/(const float3& v, const float x) { // elementwise division
	return float3(v.x/x, v.y/x, v.z/x);
}
inline float3 cross(const float3& v1, const float3& v2) {
	return float3(v1.y*v2.z-v1.z*v2.y, v1.z*v2.x-v1.x*v2.z, v1.x*v2.y-v1.y*v2.x);
}
inline float dot(const float3& v1, const float3& v2) {
	return v1.x*v2.x+v1.y*v2.y+v1.z*v2.z;
}
inline float length(const float3& v) {
	return sqrt(sq(v.x)+sq(v.y)+sq(v.z));
}
inline float3 normalize(const float3& v, const float s=1.0f) {
	const float n = s/length(v);
	return float3(v.x*n, v.y*n, v.z*n);
}
inline float angle(const float3& v1, const float3& v2) {
	return acos(dot(v1, v2)/(length(v1)*length(v2)));
}
inline float3 xyz(const float3& rtp) { // spherical coordinates (r, theta, phi) -> cartesian coordiantes (x, y, z)
	const float rsint=rtp.x*sin(rtp.y);
	return float3(rsint*cos(rtp.z), rsint*sin(rtp.z), rtp.x*cos(rtp.y));
}
inline float3 rtp(const float3& xyz) { // cartesian coordiantes (x, y, z) -> spherical coordinates (r, theta, phi)
	const float r = length(xyz);
	return float3(r, (float)acos(xyz.z/r), (float)atan2(xyz.y, xyz.x));
}
inline float plane_distance(const float3& p, const float3& plane_point, const float3& plane_normal) { // returns distance of point p from plane
	return dot(plane_normal, p-plane_point);
}

struct float3x3 { // 3D matrix type
	float xx, xy, xz, yx, yy, yz, zx, zy, zz;
	float3x3(float x) { // create diagonal matrix
		this->xx = this->yy = this->zz = x;
		this->xy = this->xz = this->yx = this->yz = this->zx = this->zy = 0.0f;
	}
	float3x3(const float xx, const float yy, const float zz) { // create diagonal matrix
		this->xx = xx;
		this->yy = yy;
		this->zz = zz;
		this->xy = this->xz = this->yx = this->yz = this->zx = this->zy = 0.0f;
	}
	float3x3(const float xx, const float xy, const float xz, const float yx, const float yy, const float yz, const float zx, const float zy, const float zz) { // create matrix
		this->xx = xx; this->xy = xy; this->xz = xz;
		this->yx = yx; this->yy = yy; this->yz = yz;
		this->zx = zx; this->zy = zy; this->zz = zz;
	}
	float3x3(const float3& v) { // create diagonal matrix from vector
		this->xx = v.x;
		this->yy = v.y;
		this->zz = v.z;
		this->xy = this->xz = this->yx = this->yz = this->zx = this->zy = 0.0f;
	}
	float3x3(const float3& v1, const float3& v2) { // create matrix from dyadic product of two vectors
		this->xx = v1.x*v2.x; this->xy = v1.x*v2.y; this->xz = v1.x*v2.z;
		this->yx = v1.y*v2.x; this->yy = v1.y*v2.y; this->yz = v1.y*v2.z;
		this->zx = v1.z*v2.x; this->zy = v1.z*v2.y; this->zz = v1.z*v2.z;
	}
	float3x3(const float3& v1, const float3& v2, const float3& v3) { // create matrix from three column vectors
		this->xx = v1.x; this->xy = v2.x; this->xz = v3.x;
		this->yx = v1.y; this->yy = v2.y; this->yz = v3.y;
		this->zx = v1.z; this->zy = v2.z; this->zz = v3.z;
	}
	float3x3(const float3& v, const float r) { // create rotation matrix for a rotation around the (normalized) axis v with angle r in radians
		const float sr=sin(r), cr=cos(r);
		this->xx = sq(v.x)+(1.0f-sq(v.x))*cr; this->xy = v.x*v.y*(1.0f-cr)-v.z*sr; this->xz = v.x*v.z*(1.0f-cr)+v.y*sr;
		this->yx = v.x*v.y*(1.0f-cr)+v.z*sr; this->yy = sq(v.y)+(1.0f-sq(v.y))*cr; this->yz = v.y*v.z*(1.0f-cr)-v.x*sr;
		this->zx = v.x*v.z*(1.0f-cr)-v.y*sr; this->zy = v.y*v.z*(1.0f-cr)+v.x*sr; this->zz = sq(v.z)+(1.0f-sq(v.z))*cr;
	}
	float3x3() = default;
	float3x3(const float3x3& m) = default; // copy constructor (elementwise copy)
	float3x3(float3x3&& x) = default; // move constructor
	~float3x3() = default;
	inline float3x3& operator=(const float3x3& m) = default; // copy assignment (elementwise copy)
	inline float3x3& operator=(float3x3&& m) = default; // move assignment
	inline float3x3& operator+=(const float3x3& m) { // elementwise addition
		this->xx += m.xx; this->xy += m.xy; this->xz += m.xz;
		this->yx += m.yx; this->yy += m.yy; this->yz += m.yz;
		this->zx += m.zx; this->zy += m.zy; this->zz += m.zz;
		return *this;
	}
	inline float3x3& operator-=(const float3x3& m) { // elementwise subtraction
		this->xx -= m.xx; this->xy -= m.xy; this->xz -= m.xz;
		this->yx -= m.yx; this->yy -= m.yy; this->yz -= m.yz;
		this->zx -= m.zx; this->zy -= m.zy; this->zz -= m.zz;
		return *this;
	}
	inline float3x3& operator*=(const float3x3& m) { // matrix product
		this->xx = this->xx*m.xx+this->xy*m.yx+this->xz*m.zx; this->xy = this->xx*m.xy+this->xy*m.yy+this->xz*m.zy; this->xz = this->xx*m.xz+this->xy*m.yz+this->xz*m.zz;
		this->yx = this->yx*m.xx+this->yy*m.yx+this->yz*m.zx; this->yy = this->yx*m.xy+this->yy*m.yy+this->yz*m.zy; this->yz = this->yx*m.xz+this->yy*m.yz+this->yz*m.zz;
		this->zx = this->zx*m.xx+this->zy*m.yx+this->zz*m.zx; this->zy = this->zx*m.xy+this->zy*m.yy+this->zz*m.zy; this->zz = this->zx*m.xz+this->zy*m.yz+this->zz*m.zz;
		return *this;
	}
	inline float3x3& operator+=(const float x) { // elementwise addition
		this->xx += x; this->xy += x; this->xz += x;
		this->yx += x; this->yy += x; this->yz += x;
		this->zx += x; this->zy += x; this->zz += x;
		return *this;
	}
	inline float3x3& operator-=(const float x) { // elementwise subtraction
		this->xx -= x; this->xy -= x; this->xz -= x;
		this->yx -= x; this->yy -= x; this->yz -= x;
		this->zx -= x; this->zy -= x; this->zz -= x;
		return *this;
	}
	inline float3x3& operator*=(const float x) { // elementwise multiplication
		this->xx *= x; this->xy *= x; this->xz *= x;
		this->yx *= x; this->yy *= x; this->yz *= x;
		this->zx *= x; this->zy *= x; this->zz *= x;
		return *this;
	}
	inline float3x3& operator/=(const float x) { // elementwise division
		this->xx /= x; this->xy /= x; this->xz /= x;
		this->yx /= x; this->yy /= x; this->yz /= x;
		this->zx /= x; this->zy /= x; this->zz /= x;
		return *this;
	}
	inline float3x3& operator+() {
		return *this;
	}
	inline const float3x3& operator+() const {
		return *this;
	}
	inline float3x3 operator-() const { // elementwise negation
		return float3x3(-this->xx, -this->xy, -this->xz, -this->yx, -this->yy, -this->yz, -this->zx, -this->zy, -this->zz);
	}
	inline float3x3 operator+(const float3x3& m) const { // elementwise addition
		return float3x3(
			this->xx+m.xx, this->xy+m.xy, this->xz+m.xz,
			this->yx+m.yx, this->yy+m.yy, this->yz+m.yz,
			this->zx+m.zx, this->zy+m.zy, this->zz+m.zz
		);
	}
	inline float3x3 operator-(const float3x3& m) const { // elementwise subtraction
		return float3x3(
			this->xx-m.xx, this->xy-m.xy, this->xz-m.xz,
			this->yx-m.yx, this->yy-m.yy, this->yz-m.yz,
			this->zx-m.zx, this->zy-m.zy, this->zz-m.zz
		);
	}
	inline float3x3 operator*(const float3x3& m) const { // matrix product
		return float3x3(
			this->xx*m.xx+this->xy*m.yx+this->xz*m.zx, this->xx*m.xy+this->xy*m.yy+this->xz*m.zy, this->xx*m.xz+this->xy*m.yz+this->xz*m.zz,
			this->yx*m.xx+this->yy*m.yx+this->yz*m.zx, this->yx*m.xy+this->yy*m.yy+this->yz*m.zy, this->yx*m.xz+this->yy*m.yz+this->yz*m.zz,
			this->zx*m.xx+this->zy*m.yx+this->zz*m.zx, this->zx*m.xy+this->zy*m.yy+this->zz*m.zy, this->zx*m.xz+this->zy*m.yz+this->zz*m.zz
		);
	}
	inline float3x3 operator^(const uint n) const { // matrix power
		float3x3 r = float3x3(float3(1.0f)); // create unit matrix
		for(uint i=0u; i<n; i++) r = r*(*this);
		return r;
	}
};
inline float3x3 operator+(const float3x3& m, const float x) { // elementwise addition
	return float3x3(
		m.xx+x, m.xy+x, m.xz+x,
		m.yx+x, m.yy+x, m.yz+x,
		m.zx+x, m.zy+x, m.zz+x
	);
}
inline float3x3 operator+(const float x, const float3x3& m) { // elementwise addition
	return m+x;
}
inline float3x3 operator-(const float3x3& m, const float x) { // elementwise subtraction
	return m+(-x);
}
inline float3x3 operator-(const float x, const float3x3& m) { // elementwise subtraction
	return float3x3(
		x-m.xx, x-m.xy, x-m.xz,
		x-m.yx, x-m.yy, x-m.yz,
		x-m.zx, x-m.zy, x-m.zz
	);
}
inline float3x3 operator*(const float3x3& m, const float x) { // elementwise multiplication
	return float3x3(
		m.xx*x, m.xy*x, m.xz*x,
		m.yx*x, m.yy*x, m.yz*x,
		m.zx*x, m.zy*x, m.zz*x
	);
}
inline float3x3 operator*(const float x, const float3x3& m) { // elementwise multiplication
	return m*x;
}
inline float3x3 operator/(const float3x3& m, const float x) { // elementwise division
	return float3x3(
		m.xx/x, m.xy/x, m.xz/x,
		m.yx/x, m.yy/x, m.yz/x,
		m.zx/x, m.zy/x, m.zz/x
	);
}
inline float3 operator*(const float3& v, const float3x3& m) { // multiply vector with matrix
	return float3(v.x*m.xx+v.y*m.yx+v.z*m.zx, v.x*m.xy+v.y*m.yy+v.z*m.zy, v.x*m.xz+v.y*m.yz+v.z*m.zz);
}
inline float3 operator*(const float3x3& m, const float3& v) { // multiply matrix with vector
		return float3(m.xx*v.x+m.xy*v.y+m.xz*v.z, m.yx*v.x+m.yy*v.y+m.yz*v.z, m.zx*v.x+m.zy*v.y+m.zz*v.z);
}
inline float3::float3(const float3x3& m) { // extract diagonal of matrix
	this->x = m.xx;
	this->y = m.yy;
	this->z = m.zz;
}
inline float3& float3::operator=(const float3x3& m) { // extract diagonal of matrix
	this->x = m.xx;
	this->y = m.yy;
	this->z = m.zz;
	return *this;
}

inline string to_string(float x); // forward-declare to_string(float x) for stringify()
struct floatNxN; // forward-declare floatNxN
struct floatN {
	uint N; // vector size is N
	float* V; // vector data
	floatN(const uint N, const float x=0.0f) { // create vector filled with zeros
		this->N = N;
		this->V = new float[N];
		for(uint i=0u; i<N; i++) this->V[i] = x;
	}
	floatN(const uint N, const float* V) {
		this->N = N;
		this->V = new float[N];
		for(uint i=0u; i<N; i++) this->V[i] = V[i];
	}
	floatN(const uint N, const floatNxN& m); // forward-declare floatNxN constructor
	floatN() = default;
	~floatN() {
		if(N==0u) delete[] V;
		N = 0u;
	}
	inline float& operator[](const uint i) {
		return V[i];
	}
	inline const float& operator[](const uint i) const {
		return V[i];
	}
	inline float operator()(const uint i) {
		return V[i];
	}
	inline const float operator()(const uint i) const {
		return V[i];
	}
	inline float* const operator()() {
		return V;
	}
	inline const float* const operator()() const {
		return V;
	}
	inline floatN& operator=(const floatN& v) {
		delete[] this->V;
		this->N = v.N;
		this->V = new float[v.N];
		for(uint i=0u; i<v.N; i++) this->V[i] = v.V[i];
		return *this;
	}
	inline floatN& operator=(const uint N) {
		delete[] this->V;
		this->N = N;
		this->V = new float[N];
		for(uint i=0u; i<N; i++) this->V[i] = 0.0f;
		return *this;
	}
	inline floatN& operator=(const floatNxN& m); // forward-declare floatNxN copy
	inline floatN& operator+=(const floatN& v) { // elementwise addition
		for(uint i=0u; i<N; i++) this->V[i] += v.V[i];
		return *this;
	}
	inline floatN& operator-=(const floatN& v) { // elementwise subtraction
		for(uint i=0u; i<N; i++) this->V[i] -= v.V[i];
		return *this;
	}
	inline floatN& operator+=(const float x) { // elementwise addition
		for(uint i=0u; i<N; i++) this->V[i] += x;
		return *this;
	}
	inline floatN& operator-=(const float x) { // elementwise subtraction
		for(uint i=0u; i<N; i++) this->V[i] -= x;
		return *this;
	}
	inline floatN& operator*=(const float x) { // elementwise multiplication
		for(uint i=0u; i<N; i++) this->V[i] *= x;
		return *this;
	}
	inline floatN& operator/=(const float x) { // elementwise division
		for(uint i=0u; i<N; i++) this->V[i] /= x;
		return *this;
	}
	inline floatN& operator+() {
		return *this;
	}
	inline const floatN& operator+() const {
		return *this;
	}
	inline floatN operator-() const { // elementwise negation
		floatN r = floatN(N);
		for(uint i=0u; i<N; i++) r.V[i] = -this->V[i];
		return r;
	}
	inline floatN operator+(const floatN& v) const { // elementwise addition
		floatN r = floatN(N);
		for(uint i=0u; i<N; i++) r.V[i] = this->V[i]+v.V[i];
		return r;
	}
	inline floatN operator-(const floatN& v) const { // elementwise subtraction
		floatN r = floatN(N);
		for(uint i=0u; i<N; i++) r.V[i] = this->V[i]-v.V[i];
		return r;
	}
	inline float operator*(const floatN& v) const { // dot product
		double r = 0.0;
		for(uint i=0u; i<N; i++) r += (double)this->V[i]*(double)v.V[i];
		return (float)r;
	}
	inline string stringify() const { // converts vector into string without spaces or newlines
		string s = "{"+to_string(V[0]);
		for(uint i=1u; i<N; i++) s += ","+to_string(V[i]);
		return s+"}";
	}
};
inline floatN operator+(const floatN& v, const float x) { // elementwise addition
	floatN r = floatN(v.N);
	for(uint i=0u; i<v.N; i++) r.V[i] = v.V[i]+x;
	return r;
}
inline floatN operator+(const float x, const floatN& v) { // elementwise addition
	return v+x;
}
inline floatN operator-(const floatN& v, const float x) { // elementwise subtraction
	return v+(-x);
}
inline floatN operator-(const float x, const floatN& v) { // elementwise subtraction
	floatN r = floatN(v.N);
	for(uint i=0u; i<v.N; i++) r.V[i] = x-v.V[i];
	return r;
}
inline floatN operator*(const floatN& v, const float x) { // elementwise multiplication
	floatN r = floatN(v.N);
	for(uint i=0u; i<v.N; i++) r.V[i] = v.V[i]*x;
	return r;
}
inline floatN operator*(const float x, const floatN& v) { // elementwise multiplication
	return v*x;
}
inline floatN operator/(const floatN& v, const float x) { // elementwise division
	floatN r = floatN(v.N);
	for(uint i=0u; i<v.N; i++) r.V[i] = v.V[i]/x;
	return r;
}

struct floatNxN {
	uint N; // matrix size is NxN
	float* M; // matrix data
	floatNxN(const uint N, const float x=0.0f) { // create matrix filled with zeros
		this->N = N;
		this->M = new float[N*N];
		for(uint i=0u; i<N*N; i++) this->M[i] = x;
	}
	floatNxN(const uint N, const float* M) {
		this->N = N;
		this->M = new float[N*N];
		for(uint i=0u; i<N*N; i++) this->M[i] = M[i];
	}
	floatNxN(const floatN& v) { // create diagonal matrix from vector
		this->N = v.N;
		this->M = new float[v.N*v.N];
		for(uint i=0u; i<v.N*v.N; i++) this->M[i] = 0.0f;
		for(uint i=0u; i<v.N; i++) this->M[v.N*i+i] = v.V[i];
	}
	floatNxN() = default;
	~floatNxN() {
		if(N==0u) delete[] M;
		N = 0u;
	}
	inline float& operator[](const uint i) {
		return M[i];
	}
	inline const float& operator[](const uint i) const {
		return M[i];
	}
	inline float operator()(const uint i) {
		return M[i];
	}
	inline const float operator()(const uint i) const {
		return M[i];
	}
	inline float operator()(const uint i, const uint j) {
		return M[N*i+j];
	}
	inline const float operator()(const uint i, const uint j) const {
		return M[N*i+j];
	}
	inline float* const operator()() {
		return M;
	}
	inline const float* const operator()() const {
		return M;
	}
	inline floatNxN& operator=(const floatNxN& m) {
		delete[] this->M;
		this->N = m.N;
		this->M = new float[m.N*m.N];
		for(uint i=0u; i<m.N*m.N; i++) this->M[i] = m.M[i];
		return *this;
	}
	inline floatNxN& operator=(const floatN& v) { // create diagonal matrix from vector
		delete[] this->M;
		this->N = v.N;
		this->M = new float[v.N*v.N];
		for(uint i=0u; i<v.N*v.N; i++) this->M[i] = 0.0f;
		for(uint i=0u; i<v.N; i++) this->M[v.N*i+i] = v.V[i];
		return *this;
	}
	inline floatNxN& operator=(const uint N) {
		delete[] this->M;
		this->N = N;
		this->M = new float[N*N];
		for(uint i=0u; i<N*N; i++) this->M[i] = 0.0f;
		return *this;
	}
	inline floatNxN& operator+=(const floatNxN& m) { // elementwise addition
		for(uint i=0u; i<N*N; i++) this->M[i] += m.M[i];
		return *this;
	}
	inline floatNxN& operator-=(const floatNxN& m) { // elementwise subtraction
		for(uint i=0u; i<N*N; i++) this->M[i] -= m.M[i];
		return *this;
	}
	inline floatNxN& operator*=(const floatNxN& m) { // matrix multiplication
		floatNxN B = floatNxN(N, this->M);
		for(uint i=0u; i<N; i++) {
			for(uint j=0u; j<N; j++) {
				double t = 0.0;
				for(uint k=0u; k<N; k++) t += (double)B.M[N*i+k]*(double)m.M[N*k+j];
				this->M[N*i+j] = (float)t;
			}
		}
		return *this;
	}
	inline floatNxN& operator+=(const float x) { // elementwise addition
		for(uint i=0u; i<N*N; i++) this->M[i] += x;
		return *this;
	}
	inline floatNxN& operator-=(const float x) { // elementwise subtraction
		for(uint i=0u; i<N*N; i++) this->M[i] -= x;
		return *this;
	}
	inline floatNxN& operator*=(const float x) { // elementwise multiplication
		for(uint i=0u; i<N*N; i++) this->M[i] *= x;
		return *this;
	}
	inline floatNxN& operator/=(const float x) { // elementwise division
		for(uint i=0u; i<N*N; i++) this->M[i] /= x;
		return *this;
	}
	inline floatNxN& operator+() {
		return *this;
	}
	inline const floatNxN& operator+() const {
		return *this;
	}
	inline floatNxN operator-() const { // elementwise negation
		floatNxN r = floatNxN(N);
		for(uint i=0u; i<N*N; i++) r.M[i] = -this->M[i];
		return r;
	}
	inline floatNxN operator+(const floatNxN& m) const { // elementwise addition
		floatNxN r = floatNxN(N);
		for(uint i=0u; i<N*N; i++) r.M[i] = this->M[i]+m.M[i];
		return r;
	}
	inline floatNxN operator-(const floatNxN& m) const { // elementwise subtraction
		floatNxN r = floatNxN(N);
		for(uint i=0u; i<N*N; i++) r.M[i] = this->M[i]-m.M[i];
		return r;
	}
	inline floatNxN operator*(const floatNxN& m) const { // matrix multiplication
		floatNxN r = floatNxN(N);
		for(uint i=0u; i<N; i++) {
			for(uint j=0u; j<N; j++) {
				double t = 0.0;
				for(uint k=0u; k<N; k++) t += (double)this->M[N*i+k]*(double)m.M[N*k+j];
				r.M[N*i+j] = (float)t;
			}
		}
		return r;
	}
	inline floatNxN operator^(const uint n) const { // matrix power
		floatNxN r = floatNxN(floatN(this->N, 1.0f)); // create unit matrix
		for(uint i=0u; i<n; i++) r = r*(*this);
		return r;
	}
	inline floatNxN invert() const { // returns inverse matrix
		const double tol = 10.0*epsilon_double;
		double* A = new double[2*N*N]; // calculating intermediate values as double is strictly necessary
		for(uint i=0u; i<N; i++) {
			for(uint j=0u; j<   N; j++) A[2u*N*i+j] = (double)M[N*i+j];
			for(uint j=N ; j<2u*N; j++) A[2u*N*i+j] = (double)(i+N==j);
		}
		for(uint k=0u; k<N-1u; k++) {
			if(fabs(A[2u*N*k+k])<=tol) {
				for(uint i=k+1u; i<N; i++) {
					if(fabs(A[2u*N*i+k])>tol) {
						for(uint j=0u; j<2u*N; j++) {
							const double t = A[2u*N*k+j];
							A[2u*N*k+j] = A[2u*N*i+j];
							A[2u*N*i+j] = t;
						}
						break;
					} else if(i+1u==N) {
						delete[] A;
						return floatNxN(N);
					}
				}
			}
			for(uint i=k+1u; i<N; i++) {
				const double t = A[2u*N*i+k]/A[2u*N*k+k];
				for(uint j=k; j<2u*N; j++) A[2u*N*i+j] -= A[2u*N*k+j]*t;
			}
		}
		double det = 1.0;
		for(uint k=0u; k<N; k++) det *= A[2u*N*k+k];
		if(fabs(det)<=tol) {
			delete[] A;
			return floatNxN(N);
		}
		for(uint k=N-1u; k>0u; k--) {
			for(int i=(int)k-1; i>=0; i--) {
				const double t = A[2u*N*i+k]/A[2u*N*k+k];
				for(uint j=k; j<2u*N; j++) A[2u*N*i+j] -= A[2u*N*k+j]*t;
			}
		}
		floatNxN r = floatNxN(N);
		for(uint i=0u; i<N; i++) {
			const double t = A[2u*N*i+i];
			for(uint j=0u; j<N; j++) r.M[N*i+j] = (float)(A[2u*N*i+N+j]/t);
		}
		delete[] A;
		return r;
	}
	inline floatNxN transpose() const { // returns transpose matrix
		floatNxN r = floatNxN(N);
		for(uint i=0u; i<N; i++) for(uint j=0u; j<N; j++) r.M[N*i+j] = M[N*j+i];
		return r;
	}
	inline string stringify() const { // converts matrix into string without spaces or newlines
		string s = "{"+to_string(M[0]);
		for(uint i=1u; i<N*N; i++) s += ","+to_string(M[i]);
		return s+"}";
	}
};
inline floatNxN operator+(const floatNxN& m, const float x) { // elementwise addition
	floatNxN r = floatNxN(m.N);
	for(uint i=0u; i<m.N*m.N; i++) r.M[i] = m.M[i]+x;
	return r;
}
inline floatNxN operator+(const float x, const floatNxN& m) { // elementwise addition
	return m+x;
}
inline floatNxN operator-(const floatNxN& m, const float x) { // elementwise subtraction
	return m+(-x);
}
inline floatNxN operator-(const float x, const floatNxN& m) { // elementwise subtraction
	floatNxN r = floatNxN(m.N);
	for(uint i=0u; i<m.N*m.N; i++) r.M[i] = x-m.M[i];
	return r;
}
inline floatNxN operator*(const floatNxN& m, const float x) { // elementwise multiplication
	floatNxN r = floatNxN(m.N);
	for(uint i=0u; i<m.N*m.N; i++) r.M[i] = m.M[i]*x;
	return r;
}
inline floatNxN operator*(const float x, const floatNxN& m) { // elementwise multiplication
	return m*x;
}
inline floatNxN operator/(const floatNxN& m, const float x) { // elementwise division
	floatNxN r = floatNxN(m.N);
	for(uint i=0u; i<m.N*m.N; i++) r.M[i] = m.M[i]/x;
	return r;
}
inline floatN operator*(const floatN& v, const floatNxN& m) { // multiply matrix with vector
	floatN r = floatN(v.N);
	for(uint i=0u; i<v.N; i++) {
		double t = 0.0;
		for(uint j=0u; j<v.N; j++) t += (double)m.M[m.N*i+j]*(double)v.V[j];
		r.V[i] = (float)t;
	}
	return r;
}
inline floatN operator*(const floatNxN& m, const floatN& v) { // multiply vector with matrix
	floatN r = floatN(m.N);
	for(uint i=0u; i<m.N; i++) {
		double t = 0.0;
		for(uint j=0u; j<m.N; j++) t += (double)v.V[j]*(double)m.M[m.N*j+i];
		r.V[i] = (float)t;
	}
	return r;
}
inline floatN::floatN(const uint N, const floatNxN& m) { // extract diagonal of matrix
	this->N = m.N;
	this->V = new float[m.N];
	for(uint i=0u; i<m.N; i++) this->V[i] = m.M[m.N*i+i];
}
inline floatN& floatN::operator=(const floatNxN& m) { // extract diagonal of matrix
	delete[] this->V;
	this->N = m.N;
	this->V = new float[m.N];
	for(uint i=0u; i<m.N; i++) this->V[i] = m.M[m.N*i+i];
	return *this;
}

struct double3x3; // forward-declare double3x3
struct double3 { // 3D vector type
	double x, y, z;
	double3(const double x, const double y, const double z) {
		this->x = x;
		this->y = y;
		this->z = z;
	}
	double3(const float x, const float y, const float z) {
		this->x = (double)x;
		this->y = (double)y;
		this->z = (double)z;
	}
	double3(const int x, const int y, const int z) {
		this->x = (double)x;
		this->y = (double)y;
		this->z = (double)z;
	}
	double3(const uint x, const uint y, const uint z) {
		this->x = (double)x;
		this->y = (double)y;
		this->z = (double)z;
	}
	double3(const double x) {
		this->x = this->y = this->z = x;
	}
	double3(const float x) {
		this->x = this->y = this->z = (double)x;
	}
	double3(const int x) {
		this->x = this->y = this->z = (double)x;
	}
	double3(const uint x) {
		this->x = this->y = this->z = (double)x;
	}
	double3(const double3x3& m); // forward-declare double3x3 constructor
	double3() = default;
	double3(const double3& v) = default; // copy constructor (elementwise copy)
	double3(double3&& v) = default; // move constructor
	~double3() = default;
	inline double3& operator=(const double3& v) = default; // copy assignment (elementwise copy)
	inline double3& operator=(double3&& v) = default; // move assignment
	inline double3& operator=(const double3x3& m); // forward-declare double3x3 copy
	inline double3& operator+=(const double3& v) { // elementwise addition
		this->x += v.x;
		this->y += v.y;
		this->z += v.z;
		return *this;
	}
	inline double3& operator-=(const double3& v) { // elementwise subtraction
		this->x -= v.x;
		this->y -= v.y;
		this->z -= v.z;
		return *this;
	}
	inline double3& operator+=(const double x) { // elementwise addition
		this->x += x;
		this->y += x;
		this->z += x;
		return *this;
	}
	inline double3& operator-=(const double x) { // elementwise subtraction
		this->x -= x;
		this->y -= x;
		this->z -= x;
		return *this;
	}
	inline double3& operator*=(const double x) { // elementwise multiplication
		this->x *= x;
		this->y *= x;
		this->z *= x;
		return *this;
	}
	inline double3& operator/=(const double x) { // elementwise division
		this->x /= x;
		this->y /= x;
		this->z /= x;
		return *this;
	}
	inline double3& operator+() {
		return *this;
	}
	inline const double3& operator+() const {
		return *this;
	}
	inline double3 operator-() const { // elementwise negation
		return double3(-this->x, -this->y, -this->z);
	}
	inline double3 operator+(const double3& v) const { // elementwise addition
		return double3(this->x+v.x, this->y+v.y, this->z+v.z);
	}
	inline double3 operator-(const double3& v) const { // elementwise subtraction
		return double3(this->x-v.x, this->y-v.y, this->z-v.z);
	}
	inline double operator*(const double3& v) const { // dot product
		return this->x*v.x+this->y*v.y+this->z*v.z;
	}
	inline bool operator==(const double3& v) const { // check for equality
		return this->x==v.x&&this->y==v.y&&this->z==v.z;
	}
	inline bool operator!=(const double3& v) const { // check for inequality
		return this->x!=v.x||this->y!=v.y||this->z!=v.z;
	}
	inline bool operator>(const double3& v) const { // compare length
		return sq(this->x)+sq(this->y)+sq(this->z)>sq(v.x)+sq(v.y)+sq(v.z);
	}
	inline bool operator<(const double3& v) const { // compare length
		return sq(this->x)+sq(this->y)+sq(this->z)<sq(v.x)+sq(v.y)+sq(v.z);
	}
	inline bool operator>=(const double3& v) const { // compare length
		return sq(this->x)+sq(this->y)+sq(this->z)>=sq(v.x)+sq(v.y)+sq(v.z);
	}
	inline bool operator<=(const double3& v) const { // compare length
		return sq(this->x)+sq(this->y)+sq(this->z)<=sq(v.x)+sq(v.y)+sq(v.z);
	}
};
inline double3 operator+(const double3& v, const double x) { // elementwise addition
	return double3(v.x+x, v.y+x, v.z+x);
}
inline double3 operator+(const double x, const double3& v) { // elementwise addition
	return v+x;
}
inline double3 operator-(const double3& v, const double x) { // elementwise subtraction
	return v+(-x);
}
inline double3 operator-(const double x, const double3& v) { // elementwise subtraction
	return double3(x-v.x, x-v.y, x-v.z);
}
inline double3 operator*(const double3& v, const double x) { // elementwise multiplication
	return double3(v.x*x, v.y*x, v.z*x);
}
inline double3 operator*(const double x, const double3& v) { // elementwise multiplication
	return v*x;
}
inline double3 operator/(const double3& v, const double x) { // elementwise division
	return double3(v.x/x, v.y/x, v.z/x);
}
inline double3 cross(const double3& v1, const double3& v2) {
	return double3(v1.y*v2.z-v1.z*v2.y, v1.z*v2.x-v1.x*v2.z, v1.x*v2.y-v1.y*v2.x);
}
inline double dot(const double3& v1, const double3& v2) {
	return v1.x*v2.x+v1.y*v2.y+v1.z*v2.z;
}
inline double length(const double3& v) {
	return sqrt(sq(v.x)+sq(v.y)+sq(v.z));
}
inline double3 normalize(const double3& v, const double s=1.0) {
	const double n = s/length(v);
	return double3(v.x*n, v.y*n, v.z*n);
}
inline double angle(const double3& v1, const double3& v2) {
	return acos(dot(v1, v2)/(length(v1)*length(v2)));
}
inline double3 xyz(const double3& rtp) { // spherical coordinates (r, theta, phi) -> cartesian coordiantes (x, y, z)
	const double rsint=rtp.x*sin(rtp.y);
	return double3(rsint*cos(rtp.z), rsint*sin(rtp.z), rtp.x*cos(rtp.y));
}
inline double3 rtp(const double3& xyz) { // cartesian coordiantes (x, y, z) -> spherical coordinates (r, theta, phi)
	const double r = length(xyz);
	return double3(r, acos(xyz.z/r), atan2(xyz.y, xyz.x));
}
inline double plane_distance(const double3& p, const double3& plane_point, const double3& plane_normal) { // returns distance of point p from plane
	return dot(plane_normal, p-plane_point);
}

struct double3x3 { // 3D matrix type
	double xx, xy, xz, yx, yy, yz, zx, zy, zz;
	double3x3(double x) { // create diagonal matrix
		this->xx = this->yy = this->zz = x;
		this->xy = this->xz = this->yx = this->yz = this->zx = this->zy = 0.0;
	}
	double3x3(const double xx, const double yy, const double zz) { // create diagonal matrix
		this->xx = xx;
		this->yy = yy;
		this->zz = zz;
		this->xy = this->xz = this->yx = this->yz = this->zx = this->zy = 0.0;
	}
	double3x3(const double xx, const double xy, const double xz, const double yx, const double yy, const double yz, const double zx, const double zy, const double zz) { // create matrix
		this->xx = xx; this->xy = xy; this->xz = xz;
		this->yx = yx; this->yy = yy; this->yz = yz;
		this->zx = zx; this->zy = zy; this->zz = zz;
	}
	double3x3(const double3& v) { // create diagonal matrix from vector
		this->xx = v.x;
		this->yy = v.y;
		this->zz = v.z;
		this->xy = this->xz = this->yx = this->yz = this->zx = this->zy = 0.0;
	}
	double3x3(const double3& v1, const double3& v2) { // create matrix from dyadic product of two vectors
		this->xx = v1.x*v2.x; this->xy = v1.x*v2.y; this->xz = v1.x*v2.z;
		this->yx = v1.y*v2.x; this->yy = v1.y*v2.y; this->yz = v1.y*v2.z;
		this->zx = v1.z*v2.x; this->zy = v1.z*v2.y; this->zz = v1.z*v2.z;
	}
	double3x3(const double3& v1, const double3& v2, const double3& v3) { // create matrix from three column vectors
		this->xx = v1.x; this->xy = v2.x; this->xz = v3.x;
		this->yx = v1.y; this->yy = v2.y; this->yz = v3.y;
		this->zx = v1.z; this->zy = v2.z; this->zz = v3.z;
	}
	double3x3(const double3& v, const double r) { // create rotation matrix for a rotation around the (normalized) axis v with angle r in radians
		const double sr=sin(r), cr=cos(r);
		this->xx = sq(v.x)+(1.0-sq(v.x))*cr; this->xy = v.x*v.y*(1.0-cr)-v.z*sr; this->xz = v.x*v.z*(1.0-cr)+v.y*sr;
		this->yx = v.x*v.y*(1.0-cr)+v.z*sr; this->yy = sq(v.y)+(1.0-sq(v.y))*cr; this->yz = v.y*v.z*(1.0-cr)-v.x*sr;
		this->zx = v.x*v.z*(1.0-cr)-v.y*sr; this->zy = v.y*v.z*(1.0-cr)+v.x*sr; this->zz = sq(v.z)+(1.0-sq(v.z))*cr;
	}
	double3x3() = default;
	double3x3(const double3x3& m) = default; // copy constructor (elementwise copy)
	double3x3(double3x3&& x) = default; // move constructor
	~double3x3() = default;
	inline double3x3& operator=(const double3x3& m) = default; // copy assignment (elementwise copy)
	inline double3x3& operator=(double3x3&& m) = default; // move assignment
	inline double3x3& operator+=(const double3x3& m) { // elementwise addition
		this->xx += m.xx; this->xy += m.xy; this->xz += m.xz;
		this->yx += m.yx; this->yy += m.yy; this->yz += m.yz;
		this->zx += m.zx; this->zy += m.zy; this->zz += m.zz;
		return *this;
	}
	inline double3x3& operator-=(const double3x3& m) { // elementwise subtraction
		this->xx -= m.xx; this->xy -= m.xy; this->xz -= m.xz;
		this->yx -= m.yx; this->yy -= m.yy; this->yz -= m.yz;
		this->zx -= m.zx; this->zy -= m.zy; this->zz -= m.zz;
		return *this;
	}
	inline double3x3& operator*=(const double3x3& m) { // matrix product
		this->xx = this->xx*m.xx+this->xy*m.yx+this->xz*m.zx; this->xy = this->xx*m.xy+this->xy*m.yy+this->xz*m.zy; this->xz = this->xx*m.xz+this->xy*m.yz+this->xz*m.zz;
		this->yx = this->yx*m.xx+this->yy*m.yx+this->yz*m.zx; this->yy = this->yx*m.xy+this->yy*m.yy+this->yz*m.zy; this->yz = this->yx*m.xz+this->yy*m.yz+this->yz*m.zz;
		this->zx = this->zx*m.xx+this->zy*m.yx+this->zz*m.zx; this->zy = this->zx*m.xy+this->zy*m.yy+this->zz*m.zy; this->zz = this->zx*m.xz+this->zy*m.yz+this->zz*m.zz;
		return *this;
	}
	inline double3x3& operator+=(const double x) { // elementwise addition
		this->xx += x; this->xy += x; this->xz += x;
		this->yx += x; this->yy += x; this->yz += x;
		this->zx += x; this->zy += x; this->zz += x;
		return *this;
	}
	inline double3x3& operator-=(const double x) { // elementwise subtraction
		this->xx -= x; this->xy -= x; this->xz -= x;
		this->yx -= x; this->yy -= x; this->yz -= x;
		this->zx -= x; this->zy -= x; this->zz -= x;
		return *this;
	}
	inline double3x3& operator*=(const double x) { // elementwise multiplication
		this->xx *= x; this->xy *= x; this->xz *= x;
		this->yx *= x; this->yy *= x; this->yz *= x;
		this->zx *= x; this->zy *= x; this->zz *= x;
		return *this;
	}
	inline double3x3& operator/=(const double x) { // elementwise division
		this->xx /= x; this->xy /= x; this->xz /= x;
		this->yx /= x; this->yy /= x; this->yz /= x;
		this->zx /= x; this->zy /= x; this->zz /= x;
		return *this;
	}
	inline double3x3& operator+() {
		return *this;
	}
	inline const double3x3& operator+() const {
		return *this;
	}
	inline double3x3 operator-() const { // elementwise negation
		return double3x3(-this->xx, -this->xy, -this->xz, -this->yx, -this->yy, -this->yz, -this->zx, -this->zy, -this->zz);
	}
	inline double3x3 operator+(const double3x3& m) const { // elementwise addition
		return double3x3(
			this->xx+m.xx, this->xy+m.xy, this->xz+m.xz,
			this->yx+m.yx, this->yy+m.yy, this->yz+m.yz,
			this->zx+m.zx, this->zy+m.zy, this->zz+m.zz
		);
	}
	inline double3x3 operator-(const double3x3& m) const { // elementwise subtraction
		return double3x3(
			this->xx-m.xx, this->xy-m.xy, this->xz-m.xz,
			this->yx-m.yx, this->yy-m.yy, this->yz-m.yz,
			this->zx-m.zx, this->zy-m.zy, this->zz-m.zz
		);
	}
	inline double3x3 operator*(const double3x3& m) const { // matrix product
		return double3x3(
			this->xx*m.xx+this->xy*m.yx+this->xz*m.zx, this->xx*m.xy+this->xy*m.yy+this->xz*m.zy, this->xx*m.xz+this->xy*m.yz+this->xz*m.zz,
			this->yx*m.xx+this->yy*m.yx+this->yz*m.zx, this->yx*m.xy+this->yy*m.yy+this->yz*m.zy, this->yx*m.xz+this->yy*m.yz+this->yz*m.zz,
			this->zx*m.xx+this->zy*m.yx+this->zz*m.zx, this->zx*m.xy+this->zy*m.yy+this->zz*m.zy, this->zx*m.xz+this->zy*m.yz+this->zz*m.zz
		);
	}
	inline double3x3 operator^(const uint n) const { // matrix power
		double3x3 r = double3x3(double3(1.0)); // create unit matrix
		for(uint i=0u; i<n; i++) r = r*(*this);
		return r;
	}
};
inline double3x3 operator+(const double3x3& m, const double x) { // elementwise addition
	return double3x3(
		m.xx+x, m.xy+x, m.xz+x,
		m.yx+x, m.yy+x, m.yz+x,
		m.zx+x, m.zy+x, m.zz+x
	);
}
inline double3x3 operator+(const double x, const double3x3& m) { // elementwise addition
	return m+x;
}
inline double3x3 operator-(const double3x3& m, const double x) { // elementwise subtraction
	return m+(-x);
}
inline double3x3 operator-(const double x, const double3x3& m) { // elementwise subtraction
	return double3x3(
		x-m.xx, x-m.xy, x-m.xz,
		x-m.yx, x-m.yy, x-m.yz,
		x-m.zx, x-m.zy, x-m.zz
	);
}
inline double3x3 operator*(const double3x3& m, const double x) { // elementwise multiplication
	return double3x3(
		m.xx*x, m.xy*x, m.xz*x,
		m.yx*x, m.yy*x, m.yz*x,
		m.zx*x, m.zy*x, m.zz*x
	);
}
inline double3x3 operator*(const double x, const double3x3& m) { // elementwise multiplication
	return m*x;
}
inline double3x3 operator/(const double3x3& m, const double x) { // elementwise division
	return double3x3(
		m.xx/x, m.xy/x, m.xz/x,
		m.yx/x, m.yy/x, m.yz/x,
		m.zx/x, m.zy/x, m.zz/x
	);
}
inline double3 operator*(const double3& v, const double3x3& m) { // multiply vector with matrix
	return double3(v.x*m.xx+v.y*m.yx+v.z*m.zx, v.x*m.xy+v.y*m.yy+v.z*m.zy, v.x*m.xz+v.y*m.yz+v.z*m.zz);
}
inline double3 operator*(const double3x3& m, const double3& v) { // multiply matrix with vector
		return double3(m.xx*v.x+m.xy*v.y+m.xz*v.z, m.yx*v.x+m.yy*v.y+m.yz*v.z, m.zx*v.x+m.zy*v.y+m.zz*v.z);
}
inline double3::double3(const double3x3& m) { // extract diagonal of matrix
	this->x = m.xx;
	this->y = m.yy;
	this->z = m.zz;
}
inline double3& double3::operator=(const double3x3& m) { // extract diagonal of matrix
	this->x = m.xx;
	this->y = m.yy;
	this->z = m.zz;
	return *this;
}

inline string to_string(double x); // forward-declare to_string(double x) for stringify()
struct doubleNxN; // forward-declare doubleNxN
struct doubleN {
	uint N; // vector size is N
	double* V; // vector data
	doubleN(const uint N, const double x=0.0) { // create vector filled with zeros
		this->N = N;
		this->V = new double[N];
		for(uint i=0u; i<N; i++) this->V[i] = x;
	}
	doubleN(const uint N, const double* V) {
		this->N = N;
		this->V = new double[N];
		for(uint i=0u; i<N; i++) this->V[i] = V[i];
	}
	doubleN(const uint N, const doubleNxN& m); // forward-declare doubleNxN constructor
	doubleN() = default;
	~doubleN() {
		if(N==0u) delete[] V;
		N = 0u;
	}
	inline double& operator[](const uint i) {
		return V[i];
	}
	inline const double& operator[](const uint i) const {
		return V[i];
	}
	inline double operator()(const uint i) {
		return V[i];
	}
	inline const double operator()(const uint i) const {
		return V[i];
	}
	inline double* const operator()() {
		return V;
	}
	inline const double* const operator()() const {
		return V;
	}
	inline doubleN& operator=(const doubleN& v) {
		delete[] this->V;
		this->N = v.N;
		this->V = new double[v.N];
		for(uint i=0u; i<v.N; i++) this->V[i] = v.V[i];
		return *this;
	}
	inline doubleN& operator=(const uint N) {
		delete[] this->V;
		this->N = N;
		this->V = new double[N];
		for(uint i=0u; i<N; i++) this->V[i] = 0.0;
		return *this;
	}
	inline doubleN& operator=(const doubleNxN& m); // forward-declare doubleNxN copy
	inline doubleN& operator+=(const doubleN& v) { // elementwise addition
		for(uint i=0u; i<N; i++) this->V[i] += v.V[i];
		return *this;
	}
	inline doubleN& operator-=(const doubleN& v) { // elementwise subtraction
		for(uint i=0u; i<N; i++) this->V[i] -= v.V[i];
		return *this;
	}
	inline doubleN& operator+=(const double x) { // elementwise addition
		for(uint i=0u; i<N; i++) this->V[i] += x;
		return *this;
	}
	inline doubleN& operator-=(const double x) { // elementwise subtraction
		for(uint i=0u; i<N; i++) this->V[i] -= x;
		return *this;
	}
	inline doubleN& operator*=(const double x) { // elementwise multiplication
		for(uint i=0u; i<N; i++) this->V[i] *= x;
		return *this;
	}
	inline doubleN& operator/=(const double x) { // elementwise division
		for(uint i=0u; i<N; i++) this->V[i] /= x;
		return *this;
	}
	inline doubleN& operator+() {
		return *this;
	}
	inline const doubleN& operator+() const {
		return *this;
	}
	inline doubleN operator-() const { // elementwise negation
		doubleN r = doubleN(N);
		for(uint i=0u; i<N; i++) r.V[i] = -this->V[i];
		return r;
	}
	inline doubleN operator+(const doubleN& v) const { // elementwise addition
		doubleN r = doubleN(N);
		for(uint i=0u; i<N; i++) r.V[i] = this->V[i]+v.V[i];
		return r;
	}
	inline doubleN operator-(const doubleN& v) const { // elementwise subtraction
		doubleN r = doubleN(N);
		for(uint i=0u; i<N; i++) r.V[i] = this->V[i]-v.V[i];
		return r;
	}
	inline double operator*(const doubleN& v) const { // dot product
		double r = 0.0;
		for(uint i=0u; i<N; i++) r += (double)this->V[i]*(double)v.V[i];
		return (double)r;
	}
	inline string stringify() const { // converts vector into string without spaces or newlines
		string s = "{"+to_string(V[0]);
		for(uint i=1u; i<N; i++) s += ","+to_string(V[i]);
		return s+"}";
	}
};
inline doubleN operator+(const doubleN& v, const double x) { // elementwise addition
	doubleN r = doubleN(v.N);
	for(uint i=0u; i<v.N; i++) r.V[i] = v.V[i]+x;
	return r;
}
inline doubleN operator+(const double x, const doubleN& v) { // elementwise addition
	return v+x;
}
inline doubleN operator-(const doubleN& v, const double x) { // elementwise subtraction
	return v+(-x);
}
inline doubleN operator-(const double x, const doubleN& v) { // elementwise subtraction
	doubleN r = doubleN(v.N);
	for(uint i=0u; i<v.N; i++) r.V[i] = x-v.V[i];
	return r;
}
inline doubleN operator*(const doubleN& v, const double x) { // elementwise multiplication
	doubleN r = doubleN(v.N);
	for(uint i=0u; i<v.N; i++) r.V[i] = v.V[i]*x;
	return r;
}
inline doubleN operator*(const double x, const doubleN& v) { // elementwise multiplication
	return v*x;
}
inline doubleN operator/(const doubleN& v, const double x) { // elementwise division
	doubleN r = doubleN(v.N);
	for(uint i=0u; i<v.N; i++) r.V[i] = v.V[i]/x;
	return r;
}

struct doubleNxN {
	uint N; // matrix size is NxN
	double* M; // matrix data
	doubleNxN(const uint N, const double x=0.0) { // create matrix filled with zeros
		this->N = N;
		this->M = new double[N*N];
		for(uint i=0u; i<N*N; i++) this->M[i] = x;
	}
	doubleNxN(const uint N, const double* M) {
		this->N = N;
		this->M = new double[N*N];
		for(uint i=0u; i<N*N; i++) this->M[i] = M[i];
	}
	doubleNxN(const doubleN& v) { // create diagonal matrix from vector
		this->N = v.N;
		this->M = new double[v.N*v.N];
		for(uint i=0u; i<v.N*v.N; i++) this->M[i] = 0.0;
		for(uint i=0u; i<v.N; i++) this->M[v.N*i+i] = v.V[i];
	}
	doubleNxN() = default;
	~doubleNxN() {
		if(N==0u) delete[] M;
		N = 0u;
	}
	inline double& operator[](const uint i) {
		return M[i];
	}
	inline const double& operator[](const uint i) const {
		return M[i];
	}
	inline double operator()(const uint i) {
		return M[i];
	}
	inline const double operator()(const uint i) const {
		return M[i];
	}
	inline double operator()(const uint i, const uint j) {
		return M[N*i+j];
	}
	inline const double operator()(const uint i, const uint j) const {
		return M[N*i+j];
	}
	inline double* const operator()() {
		return M;
	}
	inline const double* const operator()() const {
		return M;
	}
	inline doubleNxN& operator=(const doubleNxN& m) {
		delete[] this->M;
		this->N = m.N;
		this->M = new double[m.N*m.N];
		for(uint i=0u; i<m.N*m.N; i++) this->M[i] = m.M[i];
		return *this;
	}
	inline doubleNxN& operator=(const doubleN& v) { // create diagonal matrix from vector
		delete[] this->M;
		this->N = v.N;
		this->M = new double[v.N*v.N];
		for(uint i=0u; i<v.N*v.N; i++) this->M[i] = 0.0f;
		for(uint i=0u; i<v.N; i++) this->M[v.N*i+i] = v.V[i];
		return *this;
	}
	inline doubleNxN& operator=(const uint N) {
		delete[] this->M;
		this->N = N;
		this->M = new double[N*N];
		for(uint i=0u; i<N*N; i++) this->M[i] = 0.0f;
		return *this;
	}
	inline doubleNxN& operator+=(const doubleNxN& m) { // elementwise addition
		for(uint i=0u; i<N*N; i++) this->M[i] += m.M[i];
		return *this;
	}
	inline doubleNxN& operator-=(const doubleNxN& m) { // elementwise subtraction
		for(uint i=0u; i<N*N; i++) this->M[i] -= m.M[i];
		return *this;
	}
	inline doubleNxN& operator*=(const doubleNxN& m) { // matrix multiplication
		doubleNxN B = doubleNxN(N, this->M);
		for(uint i=0u; i<N; i++) {
			for(uint j=0u; j<N; j++) {
				double t = 0.0;
				for(uint k=0u; k<N; k++) t += (double)B.M[N*i+k]*(double)m.M[N*k+j];
				this->M[N*i+j] = (double)t;
			}
		}
		return *this;
	}
	inline doubleNxN& operator+=(const double x) { // elementwise addition
		for(uint i=0u; i<N*N; i++) this->M[i] += x;
		return *this;
	}
	inline doubleNxN& operator-=(const double x) { // elementwise subtraction
		for(uint i=0u; i<N*N; i++) this->M[i] -= x;
		return *this;
	}
	inline doubleNxN& operator*=(const double x) { // elementwise multiplication
		for(uint i=0u; i<N*N; i++) this->M[i] *= x;
		return *this;
	}
	inline doubleNxN& operator/=(const double x) { // elementwise division
		for(uint i=0u; i<N*N; i++) this->M[i] /= x;
		return *this;
	}
	inline doubleNxN& operator+() {
		return *this;
	}
	inline const doubleNxN& operator+() const {
		return *this;
	}
	inline doubleNxN operator-() const { // elementwise negation
		doubleNxN r = doubleNxN(N);
		for(uint i=0u; i<N*N; i++) r.M[i] = -this->M[i];
		return r;
	}
	inline doubleNxN operator+(const doubleNxN& m) const { // elementwise addition
		doubleNxN r = doubleNxN(N);
		for(uint i=0u; i<N*N; i++) r.M[i] = this->M[i]+m.M[i];
		return r;
	}
	inline doubleNxN operator-(const doubleNxN& m) const { // elementwise subtraction
		doubleNxN r = doubleNxN(N);
		for(uint i=0u; i<N*N; i++) r.M[i] = this->M[i]-m.M[i];
		return r;
	}
	inline doubleNxN operator*(const doubleNxN& m) const { // matrix multiplication
		doubleNxN r = doubleNxN(N);
		for(uint i=0u; i<N; i++) {
			for(uint j=0u; j<N; j++) {
				double t = 0.0;
				for(uint k=0u; k<N; k++) t += (double)this->M[N*i+k]*(double)m.M[N*k+j];
				r.M[N*i+j] = (double)t;
			}
		}
		return r;
	}
	inline doubleNxN operator^(const uint n) const { // matrix power
		doubleNxN r = doubleNxN(doubleN(this->N, 1.0f)); // create unit matrix
		for(uint i=0u; i<n; i++) r = r*(*this);
		return r;
	}
	inline doubleNxN invert() const { // returns inverse matrix
		const double tol = 10.0*(double)epsilon_float;
		double* A = new double[2*N*N]; // calculating intermediate values as double is strictly necessary
		for(uint i=0u; i<N; i++) {
			for(uint j=0u; j<   N; j++) A[2u*N*i+j] = M[N*i+j];
			for(uint j=N ; j<2u*N; j++) A[2u*N*i+j] = (double)(i+N==j);
		}
		for(uint k=0u; k<N-1u; k++) {
			if(fabs(A[2u*N*k+k])<=tol) {
				for(uint i=k+1u; i<N; i++) {
					if(fabs(A[2u*N*i+k])>tol) {
						for(uint j=0u; j<2u*N; j++) {
							const double t = A[2u*N*k+j];
							A[2u*N*k+j] = A[2u*N*i+j];
							A[2u*N*i+j] = t;
						}
						break;
					} else if(i+1u==N) {
						delete[] A;
						return doubleNxN(N);
					}
				}
			}
			for(uint i=k+1u; i<N; i++) {
				const double t = A[2u*N*i+k]/A[2u*N*k+k];
				for(uint j=k; j<2u*N; j++) A[2u*N*i+j] -= A[2u*N*k+j]*t;
			}
		}
		double det = 1.0;
		for(uint k=0u; k<N; k++) det *= A[2u*N*k+k];
		if(fabs(det)<=tol) {
			delete[] A;
			return doubleNxN(N);
		}
		for(uint k=N-1u; k>0u; k--) {
			for(int i=(int)k-1; i>=0; i--) {
				const double t = A[2u*N*i+k]/A[2u*N*k+k];
				for(uint j=k; j<2u*N; j++) A[2u*N*i+j] -= A[2u*N*k+j]*t;
			}
		}
		doubleNxN r = doubleNxN(N);
		for(uint i=0u; i<N; i++) {
			const double t = A[2u*N*i+i];
			for(uint j=0u; j<N; j++) r.M[N*i+j] = A[2u*N*i+N+j]/t;
		}
		delete[] A;
		return r;
	}
	inline doubleNxN transpose() const { // returns transpose matrix
		doubleNxN r = doubleNxN(N);
		for(uint i=0u; i<N; i++) for(uint j=0u; j<N; j++) r.M[N*i+j] = M[N*j+i];
		return r;
	}
	inline string stringify() const { // converts matrix into string without spaces or newlines
		string s = "{"+to_string(M[0]);
		for(uint i=1u; i<N*N; i++) s += ","+to_string(M[i]);
		return s+"}";
	}
};
inline doubleNxN operator+(const doubleNxN& m, const double x) { // elementwise addition
	doubleNxN r = doubleNxN(m.N);
	for(uint i=0u; i<m.N*m.N; i++) r.M[i] = m.M[i]+x;
	return r;
}
inline doubleNxN operator+(const double x, const doubleNxN& m) { // elementwise addition
	return m+x;
}
inline doubleNxN operator-(const doubleNxN& m, const double x) { // elementwise subtraction
	return m+(-x);
}
inline doubleNxN operator-(const double x, const doubleNxN& m) { // elementwise subtraction
	doubleNxN r = doubleNxN(m.N);
	for(uint i=0u; i<m.N*m.N; i++) r.M[i] = x-m.M[i];
	return r;
}
inline doubleNxN operator*(const doubleNxN& m, const double x) { // elementwise multiplication
	doubleNxN r = doubleNxN(m.N);
	for(uint i=0u; i<m.N*m.N; i++) r.M[i] = m.M[i]*x;
	return r;
}
inline doubleNxN operator*(const double x, const doubleNxN& m) { // elementwise multiplication
	return m*x;
}
inline doubleNxN operator/(const doubleNxN& m, const double x) { // elementwise division
	doubleNxN r = doubleNxN(m.N);
	for(uint i=0u; i<m.N*m.N; i++) r.M[i] = m.M[i]/x;
	return r;
}
inline doubleN operator*(const doubleN& v, const doubleNxN& m) { // multiply matrix with vector
	doubleN r = doubleN(v.N);
	for(uint i=0u; i<v.N; i++) {
		double t = 0.0;
		for(uint j=0u; j<v.N; j++) t += (double)m.M[m.N*i+j]*(double)v.V[j];
		r.V[i] = (double)t;
	}
	return r;
}
inline doubleN operator*(const doubleNxN& m, const doubleN& v) { // multiply vector with matrix
	doubleN r = doubleN(m.N);
	for(uint i=0u; i<m.N; i++) {
		double t = 0.0;
		for(uint j=0u; j<m.N; j++) t += (double)v.V[j]*(double)m.M[m.N*j+i];
		r.V[i] = (double)t;
	}
	return r;
}
inline doubleN::doubleN(const uint N, const doubleNxN& m) { // extract diagonal of matrix
	this->N = m.N;
	this->V = new double[m.N];
	for(uint i=0u; i<m.N; i++) this->V[i] = m.M[m.N*i+i];
}
inline doubleN& doubleN::operator=(const doubleNxN& m) { // extract diagonal of matrix
	delete[] this->V;
	this->N = m.N;
	this->V = new double[m.N];
	for(uint i=0u; i<m.N; i++) this->V[i] = m.M[m.N*i+i];
	return *this;
}

inline float plic_cube_reduced(const float V, const float n1, const float n2, const float n3) { // optimized solution from SZ and Kawano
	const float n12=n1+n2, n3V=n3*V;
	if(n12<=2.0f*n3V) return n3V+0.5f*n12; // case (5)
	const float sqn1=sq(n1), n26=6.0f*n2, v1=sqn1/n26; // after case (5) check n2>0 is true
	if(v1<=n3V && n3V<v1+0.5f*(n2-n1)) return 0.5f*(n1+sqrt(sqn1+8.0f*n2*(n3V-v1))); // case (2)
	const float V6 = n1*n26*n3V;
	if(n3V<v1) return cbrt(V6); // case (1)
	const float v3 = n3<n12 ? (sq(n3)*(3.0f*n12-n3)+sqn1*(n1-3.0f*n3)+sq(n2)*(n2-3.0f*n3))/(n1*n26) : 0.5f*n12; // after case (2) check n1>0 is true
	const float sqn12=sqn1+sq(n2), V6cbn12=V6-cb(n1)-cb(n2);
	const bool case34 = n3V<v3; // true: case (3), false: case (4)
	const float a = case34 ? V6cbn12 : 0.5f*(V6cbn12-cb(n3));
	const float b = case34 ?   sqn12 : 0.5f*(sqn12+sq(n3));
	const float c = case34 ?     n12 : 0.5f;
	const float t = sqrt(sq(c)-b);
	return c-2.0f*t*sin(0.33333334f*asin((cb(c)-0.5f*a-1.5f*b*c)/cb(t)));
}
inline float plic_cube(const float V0, const float3 n) { // unit cube - plane intersection: volume V0 in [0,1], normal vector n -> plane offset d0
	const float ax=fabs(n.x), ay=fabs(n.y), az=fabs(n.z), V=0.5f-fabs(V0-0.5f), l=ax+ay+az; // eliminate symmetry cases, normalize n using L1 norm
	const float n1 = fmin(fmin(ax, ay), az)/l;
	const float n3 = fmax(fmax(ax, ay), az)/l;
	const float n2 = fdim(1.0f, n1+n3); // ensure n2>=0
	const float d = plic_cube_reduced(V, n1, n2, n3);
	return l*copysignf(0.5f-d, V0-0.5f); // rescale result and apply symmetry for V0>0.5
}
inline float plic_cube_inverse(const float d0, const float3 n) { // unit cube - plane intersection: plane offset d0, normal vector n -> volume V0 in [0,1]
	const float n1 = fmin(fmin(fabs(n.x), fabs(n.y)), fabs(n.z)); // eliminate most cases due to symmetry
	const float n3 = fmax(fmax(fabs(n.x), fabs(n.y)), fabs(n.z));
	const float n2 = fabs(n.x)-n1+fabs(n.y)+fabs(n.z)-n3;
	const float d = 0.5f*(n1+n2+n3)-fabs(d0); // calculate PLIC with reduced symmetry, shift origin from (0.0,0.0,0.0) -> (0.5,0.5,0.5)
	float V; // 0.0<=V<=0.5
	if(fmin(n1+n2, n3)<=d && d<=n3) { // case (5)
		V = (d-0.5f*(n1+n2))/n3; // avoid division by n1 and n2
	} else if(d<n1) { // case (1)
		V = cb(d)/(6.0f*n1*n2*n3); // condition d<n1==0 is impossible if d==0.0f
	} else if(d<=n2) { // case (2)
		V = (3.0f*d*(d-n1)+sq(n1))/(6.0f*n2*n3); // avoid division by n1
	} else { // case (3) or (4)
		V = (cb(d)-cb(d-n1)-cb(d-n2)-cb(fdim(d, n3)))/(6.0f*n1*n2*n3);
	}
	return copysignf(0.5f-V, d0)+0.5f; // apply symmetry for V0>0.5
}
inline float plic_sphere(const float V0) { // intersection of a sphere with volume 1 and a plane: truncaterd volume V0 -> plane offset d0
	return 1.240701f*sin(0.5235988f-0.33333334f*atan2(2.0f*sqrt(V0-V0*V0), 2.0f*V0-1.0f)); // d0 = cbrt(6/pi)*sin(pi/6-atan2(2*sqrt(V0-V0^2),2*V0-1)/3)
}
inline float plic_sphere_inverse(const float d0) { // intersection of a sphere with volume 1 and a plane: plane offset d0 -> truncaterd volume V0
	return 1.0471976f*sq(0.6203505f+d0)*(1.240701f-d0); // r := cbrt(3/(4*pi)), v0 = pi/3*(r+d0)^2*(2*r-d0)
}

class SimplexNoise { // simplex noise in 2D/3D/4D, source: Stefan Gustavson, https://weber.itn.liu.se/~stegu/simplexnoise/SimplexNoise.java
private:
	struct float4 {
		float x, y, z, w;
		float4(const int x, const int y, const int z, const int w=0) {
			this->x = (float)x;
			this->y = (float)y;
			this->z = (float)z;
			this->w = (float)w;
		}
	};
	const float4 grad3[12] = {
		float4( 1, 1, 0), float4(-1, 1, 0), float4( 1,-1, 0), float4(-1,-1, 0),
		float4( 1, 0, 1), float4(-1, 0, 1), float4( 1, 0,-1), float4(-1, 0,-1),
		float4( 0, 1, 1), float4( 0,-1, 1), float4( 0, 1,-1), float4( 0,-1,-1)
	};
	const float4 grad4[32] = {
		float4( 0, 1, 1, 1), float4( 0, 1, 1,-1), float4( 0, 1,-1, 1), float4( 0, 1,-1,-1),
		float4( 0,-1, 1, 1), float4( 0,-1, 1,-1), float4( 0,-1,-1, 1), float4( 0,-1,-1,-1),
		float4( 1, 0, 1, 1), float4( 1, 0, 1,-1), float4( 1, 0,-1, 1), float4( 1, 0,-1,-1),
		float4(-1, 0, 1, 1), float4(-1, 0, 1,-1), float4(-1, 0,-1, 1), float4(-1, 0,-1,-1),
		float4( 1, 1, 0, 1), float4( 1, 1, 0,-1), float4( 1,-1, 0, 1), float4( 1,-1, 0,-1),
		float4(-1, 1, 0, 1), float4(-1, 1, 0,-1), float4(-1,-1, 0, 1), float4(-1,-1, 0,-1),
		float4( 1, 1, 1, 0), float4( 1, 1,-1, 0), float4( 1,-1, 1, 0), float4( 1,-1,-1, 0),
		float4(-1, 1, 1, 0), float4(-1, 1,-1, 0), float4(-1,-1, 1, 0), float4(-1,-1,-1, 0)
	};
	const uchar p[256] = {
		151,160,137, 91, 90, 15,131, 13,201, 95, 96, 53,194,233,  7,225,140, 36,103, 30, 69,142,  8, 99, 37,240, 21, 10, 23,190,  6,148,
		247,120,234, 75,  0, 26,197, 62, 94,252,219,203,117, 35, 11, 32, 57,177, 33, 88,237,149, 56, 87,174, 20,125,136,171,168, 68,175,
		 74,165, 71,134,139, 48, 27,166, 77,146,158,231, 83,111,229,122, 60,211,133,230,220,105, 92, 41, 55, 46,245, 40,244,102,143, 54,
		 65, 25, 63,161,  1,216, 80, 73,209, 76,132,187,208, 89, 18,169,200,196,135,130,116,188,159, 86,164,100,109,198,173,186,  3, 64,
		 52,217,226,250,124,123,  5,202, 38,147,118,126,255, 82, 85,212,207,206, 59,227, 47, 16, 58, 17,182,189, 28, 42,223,183,170,213,
		119,248,152,  2, 44,154,163, 70,221,153,101,155,167, 43,172,  9,129, 22, 39,253, 19, 98,108,110,79,113,224,232,178,185, 112,104,
		218,246, 97,228,251, 34,242,193,238,210,144, 12,191,179,162,241, 81, 51,145,235,249, 14,239,107, 49,192,214, 31,181,199,106,157,
		184, 84,204,176,115,121, 50, 45,127,  4,150,254,138,236,205, 93,222,114, 67, 29, 24, 72,243,141,128,195, 78, 66,215, 61,156,180
	};
	const float F2=0.5f*(sqrt(3.0f)-1.0f), G2=(3.0f-sqrt(3.0f))/6.0f;
	const float F3=1.0f/3.0f, G3=1.0f/6.0f;
	const float F4=(sqrt(5.0f)-1.0f)*0.25f, G4=(5.0f-sqrt(5.0f))*0.05f;
	uchar perm[512];
	uchar perm12[512];
	int floor(const float x) const { return (int)x-(x<=0.0f); }
	float dot(const float4& g, const float x, const float y) const { return g.x*x+g.y*y; }
	float dot(const float4& g, const float x, const float y, const float z) const { return g.x*x+g.y*y+g.z*z; }
	float dot(const float4& g, const float x, const float y, const float z, const float w) const { return g.x*x+g.y*y+g.z*z+g.w*w; }
public:
	SimplexNoise() {
		for(int i=0; i<512; i++) {
			perm[i] = p[i&255];
			perm12[i] = (uchar)(perm[i]%12);
		}
	}
	float noise(float x, float y) const { // 2D simplex noise
		float n0, n1, n2;
		float s = (x+y)*F2;
		int i=floor(x+s), j=floor(y+s);
		float t = (i+j)*G2;
		float X0=i-t, Y0=j-t;
		float x0=x-X0, y0=y-Y0;
		int i1, j1;
		if(x0>y0) { i1=1; j1=0; }
		else { i1=0; j1=1; }
		float x1=x0-  i1+     G2, y1=y0-  j1+     G2;
		float x2=x0-1.0f+2.0f*G2, y2=y0-1.0f+2.0f*G2;
		int ii=i&255, jj=j&255;
		int gi0 = perm12[ii   +perm[jj   ]];
		int gi1 = perm12[ii+i1+perm[jj+j1]];
		int gi2 = perm12[ii+ 1+perm[jj+ 1]];
		float t0 = 0.5f-x0*x0-y0*y0;
		if(t0<0) n0 = 0.0f; else { t0 *= t0; n0 = t0*t0*dot(grad3[gi0], x0, y0); }
		float t1 = 0.5f-x1*x1-y1*y1;
		if(t1<0) n1 = 0.0f; else { t1 *= t1; n1 = t1*t1*dot(grad3[gi1], x1, y1); }
		float t2 = 0.5f-x2*x2-y2*y2;
		if(t2<0) n2 = 0.0f; else { t2 *= t2; n2 = t2*t2*dot(grad3[gi2], x2, y2); }
		return 70.0f*(n0+n1+n2);
	}
	float noise(float x, float y, float z) const { // 3D simplex noise
		float n0, n1, n2, n3;
		float s = (x+y+z)*F3;
		int i=floor(x+s), j=floor(y+s), k=floor(z+s);
		float t = (i+j+k)*G3;
		float X0=i-t, Y0=j-t, Z0=k-t;
		float x0=x-X0, y0=y-Y0, z0=z-Z0;
		int i1, j1, k1, i2, j2, k2;
		if(x0>=y0) {
			if(y0>=z0)      { i1=1; j1=0; k1=0; i2=1; j2=1; k2=0; }
			else if(x0>=z0) { i1=1; j1=0; k1=0; i2=1; j2=0; k2=1; }
			else            { i1=0; j1=0; k1=1; i2=1; j2=0; k2=1; }
		} else {
			if(y0<z0)      { i1=0; j1=0; k1=1; i2=0; j2=1; k2=1; }
			else if(x0<z0) { i1=0; j1=1; k1=0; i2=0; j2=1; k2=1; }
			else           { i1=0; j1=1; k1=0; i2=1; j2=1; k2=0; }
		}
		float x1=x0-  i1+     G3, y1=y0-  j1+     G3, z1=z0-  k1+     G3;
		float x2=x0-  i2+2.0f*G3, y2=y0-  j2+2.0f*G3, z2=z0-  k2+2.0f*G3;
		float x3=x0-1.0f+3.0f*G3, y3=y0-1.0f+3.0f*G3, z3=z0-1.0f+3.0f*G3;
		int ii=i&255, jj=j&255, kk=k&255;
		int gi0 = perm12[ii   +perm[jj   +perm[kk   ]]];
		int gi1 = perm12[ii+i1+perm[jj+j1+perm[kk+k1]]];
		int gi2 = perm12[ii+i2+perm[jj+j2+perm[kk+k2]]];
		int gi3 = perm12[ii+ 1+perm[jj+ 1+perm[kk+ 1]]];
		float t0 = 0.6f-x0*x0-y0*y0-z0*z0;
		if(t0<0) n0 = 0.0f; else { t0 *= t0; n0 = t0*t0*dot(grad3[gi0], x0, y0, z0); }
		float t1 = 0.6f-x1*x1-y1*y1-z1*z1;
		if(t1<0) n1 = 0.0f; else { t1 *= t1; n1 = t1*t1*dot(grad3[gi1], x1, y1, z1); }
		float t2 = 0.6f-x2*x2-y2*y2-z2*z2;
		if(t2<0) n2 = 0.0f; else { t2 *= t2; n2 = t2*t2*dot(grad3[gi2], x2, y2, z2); }
		float t3 = 0.6f-x3*x3-y3*y3-z3*z3;
		if(t3<0) n3 = 0.0f; else { t3 *= t3; n3 = t3*t3*dot(grad3[gi3], x3, y3, z3); }
		return 32.0f*(n0+n1+n2+n3);
	}
	float noise(float x, float y, float z, float w) const { // 4D simplex noise
		float n0, n1, n2, n3, n4;
		float s = (x+y+z+w)*F4;
		int i=floor(x+s), j=floor(y+s), k=floor(z+s), l=floor(w+s);
		float t = (i+j+k+l)*G4;
		float X0=i-t, Y0=j-t, Z0=k-t, W0=l-t;
		float x0=x-X0, y0=y-Y0, z0=z-Z0, w0=w-W0;
		int rankx=0, ranky=0, rankz=0, rankw=0;
		if(x0>y0) rankx++; else ranky++;
		if(x0>z0) rankx++; else rankz++;
		if(x0>w0) rankx++; else rankw++;
		if(y0>z0) ranky++; else rankz++;
		if(y0>w0) ranky++; else rankw++;
		if(z0>w0) rankz++; else rankw++;
		int i1, j1, k1, l1, i2, j2, k2, l2, i3, j3, k3, l3;
		i1 = rankx>=3; j1 = ranky>=3; k1 = rankz>=3; l1 = rankw>=3;
		i2 = rankx>=2; j2 = ranky>=2; k2 = rankz>=2; l2 = rankw>=2;
		i3 = rankx>=1; j3 = ranky>=1; k3 = rankz>=1; l3 = rankw>=1;
		float x1=x0-  i1+     G4, y1=y0-  j1+     G4, z1=z0-  k1+     G4, w1=w0-  l1+     G4;
		float x2=x0-  i2+2.0f*G4, y2=y0-  j2+2.0f*G4, z2=z0-  k2+2.0f*G4, w2=w0-  l2+2.0f*G4;
		float x3=x0-  i3+3.0f*G4, y3=y0-  j3+3.0f*G4, z3=z0-  k3+3.0f*G4, w3=w0-  l3+3.0f*G4;
		float x4=x0-1.0f+4.0f*G4, y4=y0-1.0f+4.0f*G4, z4=z0-1.0f+4.0f*G4, w4=w0-1.0f+4.0f*G4;
		int ii=i&255, jj=j&255, kk=k&255, ll=l&255;
		int gi0 = perm[ii   +perm[jj   +perm[kk   +perm[ll   ]]]]%32;
		int gi1 = perm[ii+i1+perm[jj+j1+perm[kk+k1+perm[ll+l1]]]]%32;
		int gi2 = perm[ii+i2+perm[jj+j2+perm[kk+k2+perm[ll+l2]]]]%32;
		int gi3 = perm[ii+i3+perm[jj+j3+perm[kk+k3+perm[ll+l3]]]]%32;
		int gi4 = perm[ii+1 +perm[jj+ 1+perm[kk +1+perm[ll+ 1]]]]%32;
		float t0 = 0.6f-x0*x0-y0*y0-z0*z0-w0*w0;
		if(t0<0) n0 = 0.0f; else { t0 *= t0; n0 = t0*t0*dot(grad4[gi0], x0, y0, z0, w0); }
		float t1 = 0.6f-x1*x1-y1*y1-z1*z1-w1*w1;
		if(t1<0) n1 = 0.0f; else { t1 *= t1; n1 = t1*t1*dot(grad4[gi1], x1, y1, z1, w1); }
		float t2 = 0.6f-x2*x2-y2*y2-z2*z2-w2*w2;
		if(t2<0) n2 = 0.0f; else { t2 *= t2; n2 = t2*t2*dot(grad4[gi2], x2, y2, z2, w2); }
		float t3 = 0.6f-x3*x3-y3*y3-z3*z3-w3*w3;
		if(t3<0) n3 = 0.0f; else { t3 *= t3; n3 = t3*t3*dot(grad4[gi3], x3, y3, z3, w3); }
		float t4 = 0.6f-x4*x4-y4*y4-z4*z4-w4*w4;
		if(t4<0) n4 = 0.0f; else { t4 *= t4; n4 = t4*t4*dot(grad4[gi4], x4, y4, z4, w4); }
		return 27.0f*(n0+n1+n2+n3+n4);
	}
};

inline void split_float(float x, uint& integral, uint& decimal, int& exponent) {
	if(x>=10.0f) { // convert to base 10
		if(x>=1E32f) { x *= 1E-32f; exponent += 32; }
		if(x>=1E16f) { x *= 1E-16f; exponent += 16; }
		if(x>= 1E8f) { x *=  1E-8f; exponent +=  8; }
		if(x>= 1E4f) { x *=  1E-4f; exponent +=  4; }
		if(x>= 1E2f) { x *=  1E-2f; exponent +=  2; }
		if(x>= 1E1f) { x *=  1E-1f; exponent +=  1; }
	}
	if(x>0.0f && x<=1.0f) {
		if(x<1E-31f) { x *=  1E32f; exponent -= 32; }
		if(x<1E-15f) { x *=  1E16f; exponent -= 16; }
		if(x< 1E-7f) { x *=   1E8f; exponent -=  8; }
		if(x< 1E-3f) { x *=   1E4f; exponent -=  4; }
		if(x< 1E-1f) { x *=   1E2f; exponent -=  2; }
		if(x<  1E0f) { x *=   1E1f; exponent -=  1; }
	}
	integral = (uint)x;
	float remainder = (x-integral)*1E8f; // 8 decimal digits
	decimal = (uint)remainder;
	if(remainder-(float)decimal>=0.5f) { // correct rounding of last decimal digit
		decimal++;
		if(decimal>=100000000u) { // decimal overflow
			decimal = 0u;
			integral++;
			if(integral>=10u) { // decimal overflow causes integral overflow
				integral = 1u;
				exponent++;
			}
		}
	}
}
inline void split_double(double x, uint& integral, ulong& decimal, int& exponent) {
	if(x>=10.0) { // convert to base 10
		if(x>=1E256) { x *= 1E-256; exponent += 256; }
		if(x>=1E128) { x *= 1E-128; exponent += 128; }
		if(x>= 1E64) { x *=  1E-64; exponent +=  64; }
		if(x>= 1E32) { x *=  1E-32; exponent +=  32; }
		if(x>= 1E16) { x *=  1E-16; exponent +=  16; }
		if(x>=  1E8) { x *=   1E-8; exponent +=   8; }
		if(x>=  1E4) { x *=   1E-4; exponent +=   4; }
		if(x>=  1E2) { x *=   1E-2; exponent +=   2; }
		if(x>=  1E1) { x *=   1E-1; exponent +=   1; }
	}
	if(x>0.0 && x<=1.0) {
		if(x<1E-255) { x *=  1E256; exponent -= 256; }
		if(x<1E-127) { x *=  1E128; exponent -= 128; }
		if(x< 1E-63) { x *=   1E64; exponent -=  64; }
		if(x< 1E-31) { x *=   1E32; exponent -=  32; }
		if(x< 1E-15) { x *=   1E16; exponent -=  16; }
		if(x<  1E-7) { x *=    1E8; exponent -=   8; }
		if(x<  1E-3) { x *=    1E4; exponent -=   4; }
		if(x<  1E-1) { x *=    1E2; exponent -=   2; }
		if(x<   1E0) { x *=    1E1; exponent -=   1; }
	}
	integral = (uint)x;
	double remainder = (x-integral)*1E16; // 16 decimal digits
	decimal = (ulong)remainder;
	if(remainder-(double)decimal>=0.5) { // correct rounding of last decimal digit
		decimal++;
		if(decimal>=10000000000000000ull) { // decimal overflow
			decimal = 0ull;
			integral++;
			if(integral>=10u) { // decimal overflow causes integral overflow
				integral = 1u;
				exponent++;
			}
		}
	}
}
inline string decimal_to_string_float(uint x, int digits) {
	string r = "";
	while((digits--)>0) {
		r = (char)(x%10u+48u)+r;
		x /= 10u;
	}
	return r;
}
inline string decimal_to_string_double(ulong x, int digits) {
	string r = "";
	while((digits--)>0) {
		r = (char)(x%10ull+48ull)+r;
		x /= 10ull;
	}
	return r;
}

inline vector<string> get_main_arguments(int argc, char* argv[]) {
	return argc>1 ? vector<string>(argv+1, argv+argc) : vector<string>();
}

inline string to_string(const string& s){
	return s;
}
inline string to_string(char c) {
	return string(1, c);
}
inline string to_string(uchar c) {
	return string(1, c);
}
inline string to_string(ulong x) {
	string r = "";
	do {
		r = (char)(x%10ull+48ull)+r;
		x /= 10ull;
	} while(x);
	return r;
}
inline string to_string(slong x) {
	return x>=0ll ? to_string((ulong)x) : "-"+to_string((ulong)(-x));
}
inline string to_string(uint x) {
	string r = "";
	do {
		r = (char)(x%10u+48u)+r;
		x /= 10u;
	} while(x);
	return r;
}
inline string to_string(int x) {
	return x>=0 ? to_string((uint)x) : "-"+to_string((uint)(-x));
}
inline string to_string(float x) { // convert float to string with full precision (<string> to_string() prints only 6 decimals)
	string s = "";
	if(x<0.0f) { s += "-"; x = -x; }
	if(std::isnan(x)) return s+"NaN";
	if(std::isinf(x)) return s+"Inf";
	uint integral, decimal;
	int exponent = 0;
	split_float(x, integral, decimal, exponent);
	return s+to_string(integral)+"."+decimal_to_string_float(decimal, 8)+(exponent!=0?"E"+to_string(exponent):"");
}
inline string to_string(double x) { // convert double to string with full precision (<string> to_string() prints only 6 decimals)
	string s = "";
	if(x<0.0) { s += "-"; x = -x; }
	if(std::isnan(x)) return s+"NaN";
	if(std::isinf(x)) return s+"Inf";
	uint integral;
	ulong decimal;
	int exponent = 0;
	split_double(x, integral, decimal, exponent);
	return s+to_string(integral)+"."+decimal_to_string_double(decimal, 16)+(exponent!=0?"E"+to_string(exponent):"");
}
inline string to_string(float x, const uint decimals) { // convert float to string with specified number of decimals
	string s = "";
	if(x<0.0f) { s += "-"; x = -x; }
	if(std::isnan(x)) return s+"NaN";
	if(std::isinf(x)||x>(float)max_ulong) return s+"Inf";
	const float power = pow(10.0f, min(decimals, 8u));
	x += 0.5f/power; // rounding
	const ulong integral = (ulong)x;
	const uint decimal = (uint)((x-(float)integral)*power);
	return s+to_string(integral)+(decimals==0u ? "" : "."+decimal_to_string_float(decimal, min((int)decimals, 8)));
}
inline string to_string(double x, const uint decimals) { // convert float to string with specified number of decimals
	string s = "";
	if(x<0.0) { s += "-"; x = -x; }
	if(std::isnan(x)) return s+"NaN";
	if(std::isinf(x)||x>(double)max_ulong) return s+"Inf";
	const double power = pow(10.0, min(decimals, 16u));
	x += 0.5/power; // rounding
	const ulong integral = (ulong)x;
	const ulong decimal = (ulong)((x-(double)integral)*power);
	return s+to_string(integral)+(decimals==0u ? "" : "."+decimal_to_string_double(decimal, min((int)decimals, 16)));
}
template<class T> inline string to_string(const vector<T>& v) { // this vector notation is needed when working with Configuration_File
	string s="[";
	for(uint i=0u; i<(uint)v.size()-1u; i++) {
		s += to_string(v[i]) + ", ";
	}
	s += (uint)v.size()>0u ? to_string(v[(uint)v.size()-1u]) : "";
	return s+"]";
}

inline uint length(const string& s) {
	return (uint)s.length();
}
inline bool contains(const string& s, const string& match) {
	return s.find(match)!=string::npos;
}
inline bool contains_any(const string& s, const vector<string>& matches) {
	for(uint i=0u; i<(uint)matches.size(); i++) if(contains(s, matches[i])) return true;
	return false;
}
inline string to_lower(const string& s) {
	string r = "";
	for(uint i=0u; i<(uint)s.length(); i++) {
		const uchar c = s.at(i);
		r += c>64u&&c<91u ? c+32u : c;
	}
	return r;
}
inline string to_upper(const string& s) {
	string r = "";
	for(uint i=0u; i<(uint)s.length(); i++) {
		const uchar c = s.at(i);
		r += c>96u&&c<123u ? c-32u : c;
	}
	return r;
}
inline bool equals(const string& a, const string& b) {
	return to_lower(a)==to_lower(b);
}
inline string replace(const string& s, const string& from, const string& to) {
	string r = s;
	int p = 0;
	while((p=(int)r.find(from, p))!=string::npos) {
		r.replace(p, from.length(), to);
		p += (int)to.length();
	}
	return r;
}
inline string substring(const string& s, const uint start, uint length=max_uint) {
	return s.substr(start, min(length, (uint)s.length()-start));
}
inline string trim(const string& s) { // removes whitespace characters from beginnig and end of string s
	const int l = (int)s.length();
	int a=0, b=l-1;
	char c;
	while(a<l && ((c=s[a])==' '||c=='\t'||c=='\n'||c=='\v'||c=='\f'||c=='\r'||c=='\0')) a++;
	while(b>a && ((c=s[b])==' '||c=='\t'||c=='\n'||c=='\v'||c=='\f'||c=='\r'||c=='\0')) b--;
	return s.substr(a, 1+b-a);
}
inline bool begins_with(const string& s, const string& match) {
	if(match.size()>s.size()) return false;
	else return equal(match.begin(), match.end(), s.begin());
}
inline bool ends_with(const string& s, const string& match) {
	if(match.size()>s.size()) return false;
	else return equal(match.rbegin(), match.rend(), s.rbegin());
}
template<class T> inline bool contains(const vector<T>& v, const T& match) {
	return find(v.begin(), v.end(), match)!=v.end();
}

inline string alignl(const uint n, const string& x="") { // converts x to string with spaces behind such that length is n if x is not longer than n
	string s = x;
	for(uint i=0u; i<n; i++) s += " ";
	return s.substr(0, max(n, (uint)x.length()));
}
inline string alignr(const uint n, const string& x="") { // converts x to string with spaces in front such that length is n if x is not longer than n
	string s = "";
	for(uint i=0u; i<n; i++) s += " ";
	s += x;
	return s.substr((uint)min((int)s.length()-(int)n, (int)n), s.length());
}
template<typename T> inline string alignl(const uint n, const T x) { // converts x to string with spaces behind such that length is n if x does not have more digits than n
	return alignl(n, to_string(x));
}
template<typename T> inline string alignr(const uint n, const T x) { // converts x to string with spaces in front such that length is n if x does not have more digits than n
	return alignr(n, to_string(x));
}

inline string print_time(const double t) { // input: time in seconds, output: formatted string ___y___d__h__m__s
	uint s = (uint)(t+0.5);
	const uint y = s/31556925; s %= 31556925;
	const uint d = s/   86400; s %=    86400;
	const uint h = s/    3600; s %=     3600;
	const uint m = s/      60; s %=       60;
	const bool dy=t>31556924.5, dd=t>86399.5, dh=t>3599.5, dm=t>59.5; // display years, days, hours, minutes?
	return (dy?                                 to_string(y)+"y":"")+
	       (dd?(dy&&d<10?"00":dy&&d<100?"0":"")+to_string(d)+"d":"")+
	       (dh?(dd&&h<10?"0"               :"")+to_string(h)+"h":"")+
	       (dm?(dh&&m<10?"0"               :"")+to_string(m)+"m":"")+
	           (dm&&s<10?"0"               :"")+to_string(s)+"s"    ;
}
inline string print_percentage(const float x) {
	return alignr(3u, to_uint(100.0f*x))+"%";
}
inline string print_progress(const float x, const int n=10) {
	const int p = to_int((float)(n-7)*x);
	string s = "[";
	for(int i=0; i<p; i++) s += "=";
	for(int i=p; i<n-7; i++) s += " ";
	return s+"] "+print_percentage(x);
}

inline void print(const string& s="") {
	std::cout << s;
}
inline void println(const string& s="") {
	std::cout << s+"\n";
}
inline void reprint(const string& s="") {
	std::cout << "\r"+s;
}
inline void wait() {
	std::cin.get();
}
template<typename T> inline void println(const T& x) {
	println(to_string(x));
}

class Image {
private:
	uint w=0u, h=0u; // width, height
	int* d = nullptr; // pixel data
	bool external_pointer = false;
public:
	inline Image(const uint width, const uint height) {
		this->w = width;
		this->h = height;
		this->d = new int[width*height];
	}
	inline Image(const uint width, const uint height, int* const data) {
		this->w = width;
		this->h = height;
		this->d = data;
		this->external_pointer = true;
	}
	inline ~Image() {
		if(!external_pointer) delete[] d;
	}
	inline const uint width() const {
		return this->w;
	}
	inline const uint height() const {
		return this->h;
	}
	inline const uint length() const {
		return this->w*this->h;
	}
	inline int* data() {
		return this->d;
	}
	inline const int color(const uint x, const uint y) const {
		return this->d[x+y*this->w];
	}
	inline const int color(const uint i) const {
		return this->d[i];
	}
	inline void set_color(const uint x, const uint y, const int color) {
		this->d[x+y*this->w] = color;
	}
	inline void set_color(const uint i, const int color) {
		this->d[i] = color;
	}
};
inline Image* rescale(const Image* image, const uint newwidth, const uint newheight, Image* rescaled=nullptr) {
	if(rescaled==nullptr||rescaled->width()!=newwidth||rescaled->height()!=newheight) {
		delete rescaled;
		rescaled = new Image(newwidth, newheight);
	}
	if(newwidth<image->width()&&newheight<image->height()) { // scale down (average downscaling)
		for(uint ny=0u; ny<newheight; ny++) {
			for(uint nx=0u; nx<newwidth; nx++) {
				int r=0, g=0, b=0;
				for(uint y=ny*image->height()/newheight; y<(ny+1u)*image->height()/newheight; y++) {
					for(uint x=nx*image->width()/newwidth; x<(nx+1u)*image->width()/newwidth; x++) {
						const int color = image->color(x, y);
						r += (color>>16)&255;
						g += (color>> 8)&255;
						b +=  color     &255;
					}
				}
				const uint k = ((nx+1u)*image->width()/newwidth-nx*image->width()/newwidth)*((ny+1u)*image->height()/newheight-ny*image->height()/newheight);
				rescaled->set_color(nx, ny, ((r/k)&255)<<16|((g/k)&255)<<8|((b/k)&255));
			}
		}
	} else if(newwidth<image->width()&&newheight>=image->height()) { // x: average downscaling, y: integer upscaling
		for(uint ny=0u; ny<newheight; ny++) {
			for(uint nx=0u; nx<newwidth; nx++) {
				int r=0, g=0, b=0;
				const uint y = ny*image->height()/newheight;
				for(uint x=nx*image->width()/newwidth; x<(nx+1u)*image->width()/newwidth; x++) {
					const int color = image->color(x, y);
					r += (color>>16)&255;
					g += (color>> 8)&255;
					b +=  color     &255;
				}
				const uint k = (nx+1u)*image->width()/newwidth-nx*image->width()/newwidth;
				rescaled->set_color(nx, ny, ((r/k)&255)<<16|((g/k)&255)<<8|((b/k)&255));
			}
		}
	} else if(newwidth>=image->width()&&newheight<image->height()) { // x: integer upscaling, y: average downscaling
		for(uint ny=0u; ny<newheight; ny++) {
			for(uint nx=0u; nx<newwidth; nx++) {
				int r=0, g=0, b=0;
				const uint x = nx*image->width()/newwidth;
				for(uint y=ny*image->height()/newheight; y<(ny+1u)*image->height()/newheight; y++) {
					const int color = image->color(x, y);
					r += (color>>16)&255;
					g += (color>> 8)&255;
					b +=  color     &255;
				}
				const uint k = (ny+1u)*image->height()/newheight-ny*image->height()/newheight;
				rescaled->set_color(nx, ny, ((r/k)&255)<<16|((g/k)&255)<<8|((b/k)&255));
			}
		}
	} else { // scale up (integer upscaling)
		for(uint ny=0u; ny<newheight; ny++) {
			for(uint nx=0u; nx<newwidth; nx++) {
				const int color = image->color(nx*image->width()/newwidth, ny*image->height()/newheight);
				rescaled->set_color(nx, ny, color);
			}
		}
	}
	return rescaled;
}
inline int color(const int red, const int green, const int blue) {
	return (red&255)<<16|(green&255)<<8|(blue&255);
}
inline int color(const int red, const int green, const int blue, const int alpha) {
	return (alpha&255)<<24|(red&255)<<16|(green&255)<<8|(blue&255);
}
inline int color(const float red, const float green, const float blue) {
	return clamp((int)(255.0f*red+0.5f), 0, 255)<<16|clamp((int)(255.0f*green+0.5f), 0, 255)<<8|clamp((int)(255.0f*blue+0.5f), 0, 255);
}
inline int color(const float red, const float green, const float blue, const float alpha) {
	return clamp((int)(255.0f*alpha+0.5f), 0, 255)<<24|clamp((int)(255.0f*red+0.5f), 0, 255)<<16|clamp((int)(255.0f*green+0.5f), 0, 255)<<8|clamp((int)(255.0f*blue+0.5f), 0, 255);
}
inline int color(const float3 rgb) {
	return color(rgb.x, rgb.y, rgb.z);
}
inline int red(const int color) {
	return (color>>16)&255;
}
inline int green(const int color) {
	return (color>>8)&255;
}
inline int blue(const int color) {
	return color&255;
}
inline int alpha(const int color) {
	return (color>>24)&255;
}
inline int brightness(const int color) {
	return (red(color)+green(color)+blue(color))/3;
}
inline int grayscale(const int color) {
	const int b = brightness(color);
	return b<<16|b<<8|b;
}
inline int invert(const int color) { // invert color
	return (255-red(color))<<16|(255-green(color))<<8|(255-blue(color));
}
inline int invert_brightness(const int color) { // invert brightness, but retain color
	const int r = red(color), g=green(color), b=blue(color);
	return ::color(255-(g+b)/2, 255-(r+b)/2, 255-(r+g)/2);
}
inline int color_mul(const int c, const float x) { // c*x
	const int r = min((int)fma((float)red  (c), x, 0.5f), 255);
	const int g = min((int)fma((float)green(c), x, 0.5f), 255);
	const int b = min((int)fma((float)blue (c), x, 0.5f), 255);
	return r<<16|g<<8|b; // values are already clamped
}
inline int color_add(const int c1, const int c2) {
	const int r1=red(c1), g1=green(c1), b1=blue(c1);
	const int r2=red(c2), g2=green(c2), b2=blue(c2);
	return min(r1+r2, 255)<<16|min(g1+g2, 255)<<8|min(b1+b2, 255); // values are already clamped
}
inline int color_average(const int c1, const int c2) {
	const int r1=red(c1), g1=green(c1), b1=blue(c1);
	const int r2=red(c2), g2=green(c2), b2=blue(c2);
	return ((r1+r2)/2)<<16|((g1+g2)/2)<<8|((b1+b2)/2);
}
inline int color_mix(const int c1, const int c2, const float w) { // w*c1+(1-w)*c2
	const float3 fc1=float3((float)red(c1), (float)green(c1), (float)blue(c1)), fc2=float3((float)red(c2), (float)green(c2), (float)blue(c2));
	const float3 fcm = w*fc1+(1.0f-w)*fc2+0.5f;
	return color((int)fcm.x, (int)fcm.y, (int)fcm.z);
}
inline int color_mix_3(const int c0, const int c1, const int c2, const float w0, const float w1, const float w2) { // w1*c1+w2*c2+w3*c3, w0+w1+w2 = 1
	const float3 fc0=float3((float)red(c0), (float)green(c0), (float)blue(c0)),  fc1=float3((float)red(c1), (float)green(c1), (float)blue(c1)), fc2=float3((float)red(c2), (float)green(c2), (float)blue(c2));
	const float3 fcm = w0*fc0+w1*fc1+w2*fc2+0.5f;
	return color((int)fcm.x, (int)fcm.y, (int)fcm.z);
}
inline float3 rgb_to_hsv(const int red, const int green, const int blue) {
	const int cmax = max(max(red, green), blue);
	const int cmin = min(min(red, green), blue);
	const int d = cmax-cmin;
	const float fd = (float)d;
	const float h = 60.0f*( d==0 ? 0.0f : cmax==red ? fmod((float)(green-blue)/fd+6.0f, 6.0f) : cmax==green ? (float)(blue-red)/fd+2.0f : (float)(red-green)/fd+4.0f);
	const float s = cmax==0 ? 0.0f : fd/(float)cmax;
	const float v = (float)cmax/255.0f;
	return float3(h, s, v);
}
inline float3 rgb_to_hsv(const int color) {
	return rgb_to_hsv(red(color), green(color), blue(color));
}
inline int hsv_to_rgb(const float h, const float s, const float v) {
	const float c = v*s;
	const float x = c*(1.0f-fabs(fmod(h/60.0f, 2.0f)-1.0f));
	const float m = v-c;
	float r=0.0f, g=0.0f, b=0.0f;
	if(0.0f<=h&&h<60.0f) { r = c; g = x; }
	else if(h<120.0f) { r = x; g = c; }
	else if(h<180.0f) { g = c; b = x; }
	else if(h<240.0f) { g = x; b = c; }
	else if(h<300.0f) { r = x; b = c; }
	else if(h<360.0f) { r = c; b = x; }
	return color(r+m, g+m, b+m);
}
inline int hsv_to_rgb(const float3& hsv) {
	return hsv_to_rgb(hsv.x, hsv.y, hsv.z);
}
inline int colorscale_rainbow(float x) { // coloring scheme (float [0, 1] -> int color)
	x = clamp(6.0f*(1.0f-x), 0.0f, 6.0f);
	float r=0.0f, g=0.0f, b=0.0f; // black
	if(x<1.2f) { // red - yellow
		r = 1.0f;
		g = x*0.83333333f;
	} else if(x>=1.2f&&x<2.0f) { // yellow - green
		r = 2.5f-x*1.25f;
		g = 1.0f;
	} else if(x>=2.0f&&x<3.0f) { // green - cyan
		g = 1.0f;
		b = x-2.0f;
	} else if(x>=3.0f&&x<4.0f) { // cyan - blue
		g = 4.0f-x;
		b = 1.0f;
	} else if(x>=4.0f&&x<5.0f) { // blue - violet
		r = x*0.4f-1.6f;
		b = 3.0f-x*0.5f;
	} else { // violet - black
		r = 2.4f-x*0.4f;
		b = 3.0f-x*0.5f;
	}
	return color(r, g, b);
}
inline int colorscale_iron(float x) { // coloring scheme (float [0, 1] -> int color)
	x = clamp(4.0f*(1.0f-x), 0.0f, 4.0f);
	float r=1.0f, g=0.0f, b=0.0f;
	if(x<0.66666667f) { // white - yellow
		g = 1.0f;
		b = 1.0f-x*1.5f;
	} else if(x<2.0f) { // yellow - red
		g = 1.5f-x*0.75f;
	} else if(x<3.0f) { // red - violet
		r = 2.0f-x*0.5f;
		b = x-2.0f;
	} else { // violet - black
		r = 2.0f-x*0.5f;
		b = 4.0f-x;
	}
	return color(r, g, b);
}
inline int colorscale_twocolor(float x) { // coloring scheme (float [0, 1] -> int color)
	return x>0.5f ? color_mix(0xFFAA00, 0x181818, clamp(2.0f*x-1.0f, 0.0f, 1.0f)) : color_mix(0x181818, 0x0080FF, clamp(2.0f*x, 0.0f, 1.0f)); // red - gray - blue
}

#define color_black      0
#define color_dark_blue  1
#define color_dark_green 2
#define color_light_blue 3
#define color_dark_red   4
#define color_magenta    5
#define color_orange     6
#define color_light_gray 7
#define color_gray       8
#define color_blue       9
#define color_green     10
#define color_cyan      11
#define color_red       12
#define color_pink      13
#define color_yellow    14
#define color_white     15

#ifdef UTILITIES_CONSOLE_INPUT
#define UTILITIES_CONSOLE_COLOR
#endif // UTILITIES_CONSOLE_INPUT
#ifdef UTILITIES_CONSOLE_COLOR
#if defined(_WIN32)
#ifndef UTILITIES_REGEX
#include <algorithm> // for transform() in get_exe_path()
#endif // UTILITIES_REGEX
#define WIN32_LEAN_AND_MEAN
#define VC_EXTRALEAN
#include <Windows.h> // for displaying colors and getting console size
#undef min
#undef max
#elif defined(__linux__)||defined(__APPLE__)
#include <sys/ioctl.h> // for getting console size
#include <unistd.h> // for getting path of executable
#if defined(__APPLE__)
#include <mach-o/dyld.h> // _NSGetExecutablePath
#endif // Apple
#else // Linux or Apple
#undef UTILITIES_CONSOLE_COLOR
#endif // Windows/Linux
#endif // UTILITIES_CONSOLE_COLOR
#ifdef UTILITIES_CONSOLE_COLOR
inline string get_exe_path() { // returns path where executable is located, ends with a "/"
	string path = "";
#if defined(_WIN32)
	wchar_t wc[260] = {0};
	GetModuleFileNameW(NULL, wc, 260);
	std::wstring ws(wc);
	transform(ws.begin(), ws.end(), back_inserter(path), [](wchar_t c) { return (char)c; });
	path = replace(path, "\\", "/");
#elif defined(__APPLE__)
	uint length = 0u;
	_NSGetExecutablePath(nullptr, &length);
	path.resize(length+1u, 0);
	_NSGetExecutablePath(path.data(), &length);
#else // Linux
	char c[260];
	int length = (int)readlink("/proc/self/exe", c, 260);
	path = string(c, length>0 ? length : 0);
#endif // Windows/Linux
	return path.substr(0, path.rfind('/')+1);
}
inline void get_console_size(uint& width, uint& height) {
#if defined(_WIN32)
	static const HANDLE handle = GetStdHandle(STD_OUTPUT_HANDLE);
	CONSOLE_SCREEN_BUFFER_INFO csbi;
	GetConsoleScreenBufferInfo(handle, &csbi);
	width = (uint)(csbi.srWindow.Right-csbi.srWindow.Left+1); // (uint)(csbi.dwSize.X); gives size of screen buffer
	height = (uint)(csbi.srWindow.Bottom-csbi.srWindow.Top+1); // (uint)(csbi.dwSize.Y); gives size of screen buffer
#else // Linux
	struct winsize w;
	ioctl(fileno(stdout), TIOCGWINSZ, &w);
	width = (uint)(w.ws_col);
	height = (uint)(w.ws_row);
#endif // Windows/Linux
}
inline void get_console_font_size(uint& width, uint& height) {
#if defined(_WIN32)
	static const HANDLE handle = GetStdHandle(STD_OUTPUT_HANDLE);
	CONSOLE_FONT_INFO cfi;
	GetCurrentConsoleFont(handle, false, &cfi);
	width = (uint)(cfi.dwFontSize.X);
	height = (uint)(cfi.dwFontSize.Y);
#else // Linux
	//struct winsize w;
	//ioctl(fileno(stdout), TIOCGWINSZ, &w);
	width = 8u;//(uint)(w.ws_xpixel/w.ws_col);
	height = 16u;//(uint)(w.ws_ypixel/w.ws_row);
#endif // Windows/Linux
}
inline void set_console_cursor(const uint x, const uint y) {
#if defined(_WIN32)
	static const HANDLE handle = GetStdHandle(STD_OUTPUT_HANDLE);
	SetConsoleCursorPosition(handle, {(short)x, (short)y});
#else // Linux
	std::cout << "\033["+to_string(y+1u)+";"+to_string(x+1u)+"f";
#endif // Windows/Linux
}
inline void show_console_cursor(const bool show) {
#if defined(_WIN32)
	static const HANDLE handle = GetStdHandle(STD_OUTPUT_HANDLE);
	CONSOLE_CURSOR_INFO cci;
	GetConsoleCursorInfo(handle, &cci);
	cci.bVisible = show; // show/hide cursor
	SetConsoleCursorInfo(handle, &cci);
#else // Linux
	std::cout << (show ? "\033[?25h" : "\033[?25l"); // show/hide cursor
#endif // Windows/Linux
}
inline void clear_console() {
#if defined(_WIN32)
	static const HANDLE handle = GetStdHandle(STD_OUTPUT_HANDLE);
	CONSOLE_SCREEN_BUFFER_INFO csbi;
	const COORD topLeft = { 0, 0 };
	std::cout.flush(); // cout uses a buffer to batch writes to the underlying console, we don't want stale buffered text to randomly be written out
	GetConsoleScreenBufferInfo(handle, &csbi); // figure out the current width and height of the console window
	const DWORD length = csbi.dwSize.X*csbi.dwSize.Y;
	DWORD written;
	FillConsoleOutputCharacter(handle, TEXT(' '), length, topLeft, &written); // flood-fill the console with spaces to clear it
	FillConsoleOutputAttribute(handle, csbi.wAttributes, length, topLeft, &written); // reset attributes of every character to default, this clears all background colour formatting
	SetConsoleCursorPosition(handle, topLeft); // move the cursor back to the top left for the next sequence of writes
#else // Linux
	std::cout << "\033[2J";
#endif // Windows/Linux
}
inline int get_console_color(const int color) {
	const int r=(color>>16)&255, g=(color>>8)&255, b=color&255;
	const int matches[16] = {
		sq(r-  0)+sq(g-  0)+sq(b-  0), // color_black      0   0   0   0
		sq(r-  0)+sq(g- 55)+sq(b-218), // color_dark_blue  1   0  55 218
		sq(r- 19)+sq(g-161)+sq(b- 14), // color_dark_green 2  19 161  14
		sq(r- 58)+sq(g-150)+sq(b-221), // color_light_blue 3  58 150 221
		sq(r-197)+sq(g- 15)+sq(b- 31), // color_dark_red   4 197  15  31
		sq(r-136)+sq(g- 23)+sq(b-152), // color_magenta    5 136  23 152
		sq(r-193)+sq(g-156)+sq(b-  0), // color_orange     6 193 156   0
		sq(r-204)+sq(g-204)+sq(b-204), // color_light_gray 7 204 204 204
		sq(r-118)+sq(g-118)+sq(b-118), // color_gray       8 118 118 118
		sq(r- 59)+sq(g-120)+sq(b-255), // color_blue       9  59 120 255
		sq(r- 22)+sq(g-198)+sq(b- 12), // color_green     10  22 198  12
		sq(r- 97)+sq(g-214)+sq(b-214), // color_cyan      11  97 214 214
		sq(r-231)+sq(g- 72)+sq(b- 86), // color_red       12 231  72  86
		sq(r-180)+sq(g-  0)+sq(b-158), // color_pink      13 180   0 158
		sq(r-249)+sq(g-241)+sq(b-165), // color_yellow    14 249 241 165
		sq(r-242)+sq(g-242)+sq(b-242)  // color_white     15 242 242 242
	};
	int m=195075, k=0;
	for(int i=0; i<16; i++) if(matches[i]<m) m = matches[k=i];
	return k;
}
inline ushort get_console_color_dither(const int color) {
	const int r=(color>>16)&255, g=(color>>8)&255, b=color&255;
	const int red  [16] = {  0,  0, 19, 58,197,136,193,204,118, 59, 22, 97,231,180,249,242};
	const int green[16] = {  0, 55,161,150, 15, 23,156,204,118,120,198,214, 72,  0,241,242};
	const int blue [16] = {  0,218, 14,221, 31,152,  0,204,118,255, 12,214, 86,158,165,242};
	int m=195075, k=0;
	for(int i=0; i<16; i++) {
		for(int j=0; j<16; j++) {
			const int mixred=(red[i]+6*red[j])/7, mixgreen=(green[i]+6*green[j])/7, mixblue=(blue[i]+6*blue[j])/7; // (char)176: pixel ratio 1:6
			const int match = sq(r-mixred)+sq(g-mixgreen)+sq(b-mixblue);
			if(match<m) {
				m = match;
				k = i<<4|j;
			}
		}
	}
	for(int i=0; i<16; i++) {
		for(int j=0; j<16; j++) {
			const int mixred=(2*red[i]+5*red[j])/7, mixgreen=(2*green[i]+5*green[j])/7, mixblue=(2*blue[i]+5*blue[j])/7; // (char)177: pixel ratio 2:5
			const int match = sq(r-mixred)+sq(g-mixgreen)+sq(b-mixblue);
			if(match<m) {
				m = match;
				k = 1<<8|i<<4|j;
			}
		}
	}
	for(int i=0; i<16; i++) {
		for(int j=0; j<i; j++) {
			const int mixred=(red[i]+red[j])/2, mixgreen=(green[i]+green[j])/2, mixblue=(blue[i]+blue[j])/2; // (char)178: pixel ratio 1:2
			const int match = sq(r-mixred)+sq(g-mixgreen)+sq(b-mixblue);
			if(match<m) {
				m = match;
				k = 2<<8|i<<4|j;
			}
		}
	}
	return (ushort)k;
}
inline string get_textcolor_code(const int textcolor) { // Linux only
	switch(textcolor) {
		case  0: return "30"; // color_black      0
		case  1: return "34"; // color_dark_blue  1
		case  2: return "32"; // color_dark_green 2
		case  3: return "36"; // color_light_blue 3
		case  4: return "31"; // color_dark_red   4
		case  5: return "35"; // color_magenta    5
		case  6: return "33"; // color_orange     6
		case  7: return "37"; // color_light_gray 7
		case  8: return "90"; // color_gray       8
		case  9: return "94"; // color_blue       9
		case 10: return "92"; // color_green     10
		case 11: return "96"; // color_cyan      11
		case 12: return "91"; // color_red       12
		case 13: return "95"; // color_pink      13
		case 14: return "93"; // color_yellow    14
		case 15: return "97"; // color_white     15
		default: return "37";
	}
}
inline string get_backgroundcolor_code(const int backgroundcolor) { // Linux only
	switch(backgroundcolor) {
		case  0: return  "40"; // color_black      0
		case  1: return  "44"; // color_dark_blue  1
		case  2: return  "42"; // color_dark_green 2
		case  3: return  "46"; // color_light_blue 3
		case  4: return  "41"; // color_dark_red   4
		case  5: return  "45"; // color_magenta    5
		case  6: return  "43"; // color_orange     6
		case  7: return  "47"; // color_light_gray 7
		case  8: return "100"; // color_gray       8
		case  9: return "104"; // color_blue       9
		case 10: return "102"; // color_green     10
		case 11: return "106"; // color_cyan      11
		case 12: return "101"; // color_red       12
		case 13: return "105"; // color_pink      13
		case 14: return "103"; // color_yellow    14
		case 15: return "107"; // color_white     15
		default: return  "40";
	}
}
inline string get_print_color(const int textcolor) { // Linux only
	return "\033["+get_textcolor_code(textcolor)+"m";
}
inline string get_print_color(const int textcolor, const int backgroundcolor) { // Linux only
	return "\033["+get_textcolor_code(textcolor)+";"+get_backgroundcolor_code(backgroundcolor)+"m";
}
inline void print_color(const int textcolor) {
#if defined(_WIN32)
	static const HANDLE handle = GetStdHandle(STD_OUTPUT_HANDLE);
	SetConsoleTextAttribute(handle, textcolor);
#else // Linux
	std::cout << get_print_color(textcolor);
#endif // Windows/Linux
}
inline void print_color(const int textcolor, const int backgroundcolor) {
#if defined(_WIN32)
	static const HANDLE handle = GetStdHandle(STD_OUTPUT_HANDLE);
	SetConsoleTextAttribute(handle, backgroundcolor<<4|textcolor);
#else // Linux
	std::cout << get_print_color(textcolor, backgroundcolor);
#endif // Windows/Linux
}
inline void print_color_reset() {
#if defined(_WIN32)
	static const HANDLE handle = GetStdHandle(STD_OUTPUT_HANDLE);
	SetConsoleTextAttribute(handle, 7); // reset color
#else // Linux
	std::cout << "\033[0m"; // reset color
#endif // Windows/Linux
}
inline void print(const string& s, const int textcolor) {
	print_color(textcolor);
	std::cout << s;
	print_color_reset();
}
inline void print(const string& s, const int textcolor, const int backgroundcolor) {
	print_color(textcolor, backgroundcolor);
	std::cout << s;
	print_color_reset();
}
inline void print_no_reset(const string& s, const int textcolor) { // print with color, but don't reset color afterwards (faster)
	print_color(textcolor);
	std::cout << s;
}
inline void print_no_reset(const string& s, const int textcolor, const int backgroundcolor) { // print with color, but don't reset color afterwards (faster)
	print_color(textcolor, backgroundcolor);
	std::cout << s;
}
static Image* print_image_rescaled = nullptr;
inline void print_image(const Image* image, const uint textwidth=0u, const uint textheight=0u) {
	uint newwidth=(uint)(CONSOLE_WIDTH), newheight=(uint)(CONSOLE_WIDTH)*9u/16u;
	if(textwidth==0u&&textheight==0u) {
		uint fontwidth=8u, fontheight=16u;
		get_console_size(newwidth, newheight);
		get_console_font_size(fontwidth, fontheight);
		newwidth--;
		newheight = newwidth*image->height()*2u*fontwidth/(image->width()*fontheight);
	} else {
		newwidth = textwidth;
		newheight = 2u*textheight;
	}
	print_image_rescaled = rescale(image, newwidth, newheight, print_image_rescaled);
#if defined(_WIN32)
	const string s = string("")+(char)223; // trick to double vertical resolution: use graphic character
	for(uint y=0u; y<newheight-1u; y+=2u) {
		int ltc = get_console_color(print_image_rescaled->color(0u, y   ));
		int lbc = get_console_color(print_image_rescaled->color(0u, y+1u));
		string segment = s;
		for(uint x=1; x<newwidth; x++) {
			const int tc = get_console_color(print_image_rescaled->color(x, y   ));
			const int bc = get_console_color(print_image_rescaled->color(x, y+1u));
			if(tc==ltc&&bc==lbc) { // still the same color
				segment += s; // skip color change command if color is the same for consecutive pixels
			} else { // color has changed
				print_no_reset(segment, ltc, lbc); // only switch color once before each segment
				segment = s; // reset segment
				ltc = tc;
				lbc = bc;
			}
		}
		print(segment, ltc, lbc); // print last segment, then reset color
		println();
	}
#else // Linux
	const string s = "\u2580"; // trick to double vertical resolution: use graphic character
	string r = ""; // append color changes to a string, print string in the end (much faster)
	for(uint y=0u; y<newheight-1u; y+=2u) {
		int ltc=-1, lbc=-1;
		for(uint x=0u; x<newwidth; x++) {
			const int tc = get_console_color(print_image_rescaled->color(x, y   ));
			const int bc = get_console_color(print_image_rescaled->color(x, y+1u));
			if(tc==ltc&&bc==lbc) { // still the same color
				r += s; // skip color change command if color is the same for consecutive pixels
			} else { // color has changed
				r += get_print_color(tc, bc)+s; // only switch color once before each segment
				ltc = tc;
				lbc = bc;
			}
		}
		r += "\033[0m\n"; // reset color after every line and add new line
	}
	print(r);
#endif // Windows/Linux
}
inline void print_image_bw(const Image* image, const uint textwidth=0u, const uint textheight=0u) { // black and white mode
	uint newwidth=(uint)(CONSOLE_WIDTH), newheight=(uint)(CONSOLE_WIDTH)*9u/16u;
	if(textwidth==0u&&textheight==0u) {
		uint fontwidth=8u, fontheight=16u;
		get_console_size(newwidth, newheight);
		get_console_font_size(fontwidth, fontheight);
		newwidth--;
		newheight = newwidth*image->height()*2u*fontwidth/(image->width()*fontheight);
	} else {
		newwidth = textwidth;
		newheight = 2u*textheight;
	}
	print_image_rescaled = rescale(image, newwidth, newheight, print_image_rescaled);
	const string bb = " ";
#if defined(_WIN32)
	const string ww = string("")+(char)219; // trick to double vertical resolution: use graphic characters
	const string bw = string("")+(char)220;
	const string wb = string("")+(char)223;
#else // Linux
	const string ww = "\u2588"; // trick to double vertical resolution: use graphic characters
	const string bw = "\u2584";
	const string wb = "\u2580";
#endif // Windows/Linux
	string r = ""; // append to a string, print string in the end (much faster)
	for(uint y=0u; y<newheight-1u; y+=2u) {
		for(uint x=0u; x<newwidth; x++) {
			const bool t = brightness(print_image_rescaled->color(x, y   ))>127;
			const bool b = brightness(print_image_rescaled->color(x, y+1u))>127;
			r += t ? b ? ww : wb : b ? bw : bb;
		}
		r += "\n"; // add new line
	}
	print(r);
}
#ifdef UTILITIES_CONSOLE_DITHER_LOOKUP
static bool print_image_dither_lookup_initialized = false;
static ushort* print_image_dither_lookup = nullptr;
inline void print_image_dither_initialize_lookup() { // initialize lookup table parallelized (much faster)
	if(!print_image_dither_lookup_initialized) {
		print_image_dither_lookup = new ushort[16777216];
		parallel_for(16777216u, [&](uint n) {
			print_image_dither_lookup[n] = get_console_color_dither((int)n);
		});
		print_image_dither_lookup_initialized = true;
	}
}
#endif // UTILITIES_CONSOLE_DITHER_LOOKUP
inline void print_image_dither(const Image* image, const uint textwidth=0u, const uint textheight=0u) {
#ifdef UTILITIES_CONSOLE_DITHER_LOOKUP
	print_image_dither_initialize_lookup();
#endif // UTILITIES_CONSOLE_DITHER_LOOKUP
	uint newwidth=(uint)(CONSOLE_WIDTH), newheight=(uint)(CONSOLE_WIDTH)*9u/16u;
	if(textwidth==0u&&textheight==0u) {
		uint fontwidth=8u, fontheight=16u;
		get_console_size(newwidth, newheight);
		get_console_font_size(fontwidth, fontheight);
		newwidth--;
		newheight = newwidth*image->height()*fontwidth/(image->width()*fontheight);
	} else {
		newwidth = textwidth;
		newheight = textheight;
	}
	print_image_rescaled = rescale(image, newwidth, newheight, print_image_rescaled);
#if defined(_WIN32)
	const string s[3] = { string("")+(char)176, string("")+(char)177, string("")+(char)178 };
	for(uint y=0u; y<newheight; y++) {
#ifdef UTILITIES_CONSOLE_DITHER_LOOKUP
		const int dither = (int)print_image_dither_lookup[print_image_rescaled->color(0u, y)];
#else // UTILITIES_CONSOLE_DITHER_LOOKUP
		const int dither = (int)get_console_color_dither(print_image_rescaled->color(0u, y));
#endif // UTILITIES_CONSOLE_DITHER_LOOKUP
		const int shade = dither>>8;
		int ltc=(dither>>4)&0xF, lbc=dither&0xF;
		string segment = s[shade];
		for(uint x=1u; x<newwidth; x++) {
#ifdef UTILITIES_CONSOLE_DITHER_LOOKUP
			const int dither = (int)print_image_dither_lookup[print_image_rescaled->color(x, y)];
#else // UTILITIES_CONSOLE_DITHER_LOOKUP
			const int dither = (int)get_console_color_dither(print_image_rescaled->color(x, y));
#endif // UTILITIES_CONSOLE_DITHER_LOOKUP
			const int shade = dither>>8;
			const int tc=(dither>>4)&0xF, bc=dither&0xF;
			if(tc!=ltc||bc!=lbc) { // color has changed
				print_no_reset(segment, ltc, lbc); // only switch color once before each segment
				segment = ""; // reset segment
				ltc = tc;
				lbc = bc;
			}
			segment += s[shade];
		}
		print(segment, ltc, lbc); // print last segment, then reset color
		println();
	}
#else // Linux
	const string s[3] = { "\u2591", "\u2592", "\u2593" };
	string r = ""; // append color changes to a string, print string in the end (much faster)
	for(uint y=0u; y<newheight; y++) {
		int ltc=-1, lbc=-1;
		for(uint x=0u; x<newwidth; x++) {
#ifdef UTILITIES_CONSOLE_DITHER_LOOKUP
			const int dither = (int)print_image_dither_lookup[print_image_rescaled->color(x, y)];
#else // UTILITIES_CONSOLE_DITHER_LOOKUP
			const int dither = (int)get_console_color_dither(print_image_rescaled->color(x, y));
#endif // UTILITIES_CONSOLE_DITHER_LOOKUP
			const int shade = dither>>8;
			const int tc=(dither>>4)&0xF, bc=dither&0xF;
			if(tc!=ltc||bc!=lbc) { // color has changed
				r += get_print_color(tc, bc); // only switch color once before each segment
				ltc = tc;
				lbc = bc;
			}
			r += s[shade];
		}
		r += "\033[0m\n"; // reset color after every line and add new line
	}
	print(r);
#endif // Windows/Linux
}
static Image* print_video_lastframe = nullptr;
inline void print_video_dither(const Image* image, const uint textwidth=0u, const uint textheight=0u) { // optimized for video (only draws pixels that differ from last frame)
#ifdef UTILITIES_CONSOLE_DITHER_LOOKUP
	print_image_dither_initialize_lookup();
#endif // UTILITIES_CONSOLE_DITHER_LOOKUP
	uint newwidth=(uint)(CONSOLE_WIDTH), newheight=(uint)(CONSOLE_WIDTH)*9u/16u;
	if(textwidth==0u&&textheight==0u) {
		uint fontwidth=8u, fontheight=16u;
		get_console_size(newwidth, newheight);
		get_console_font_size(fontwidth, fontheight);
		newwidth--;
		newheight = newwidth*image->height()*fontwidth/(image->width()*fontheight);
	} else {
		newwidth = textwidth;
		newheight = textheight;
	}
	print_image_rescaled = rescale(image, newwidth, newheight, print_image_rescaled); // do resolution downsampling
	for(uint i=0u; i<newwidth*newheight; i++) { // do color subsampling
#ifdef UTILITIES_CONSOLE_DITHER_LOOKUP
		print_image_rescaled->set_color(i, (int)print_image_dither_lookup[print_image_rescaled->color(i)]);
#else // UTILITIES_CONSOLE_DITHER_LOOKUP
		print_image_rescaled->set_color(i, (int)get_console_color_dither(print_image_rescaled->color(i)));
#endif // UTILITIES_CONSOLE_DITHER_LOOKUP
	}
	bool print_video_lastframe_is_new = false;
	if(print_video_lastframe==nullptr||print_video_lastframe->width()!=newwidth||print_video_lastframe->height()!=newheight) {
		delete print_video_lastframe;
		print_video_lastframe = new Image(newwidth, newheight);
		print_video_lastframe_is_new = true;
	}
	uint changed_pixels = 0u; // reuse print_video_lastframe as list of all changed pixels
	for(uint i=0u; i<newwidth*newheight; i++) {
		if(print_video_lastframe_is_new || print_image_rescaled->color(i)!=print_video_lastframe->color(i)) print_video_lastframe->set_color(changed_pixels++, (int)i);
	}
	if(changed_pixels>0u) { // at least one pixel has changed
#if defined(_WIN32)
		const string s[3] = { string("")+(char)176, string("")+(char)177, string("")+(char)178 };
		const uint n = (uint)print_video_lastframe->color(0); // do first entry manually
		uint lx=n%newwidth, ly=n/newwidth;
		const int dither = print_image_rescaled->color(lx, ly);
		const int shade = dither>>8;
		int ltc=(dither>>4)&0xF, lbc=dither&0xF;
		set_console_cursor(lx, ly);
		string segment = s[shade]; // reset segment
		for(uint i=1u; i<changed_pixels; i++) {
			const uint n = (uint)print_video_lastframe->color(i); // get index of changed pixel
			const uint x=n%newwidth, y=n/newwidth; // get coordinates of changed pixel
			const int dither = print_image_rescaled->color(x, y); // get color of changed pixel
			const int shade = dither>>8;
			const int tc=(dither>>4)&0xF, bc=dither&0xF;
			if(tc!=ltc||bc!=lbc||x!=lx+1u||y!=ly) { // color has changed or cursor has moved
				print_no_reset(segment, ltc, lbc); // print old segment, then change color
				segment = ""; // reset segment
				ltc = tc;
				lbc = bc;
			}
			if(x!=lx+1u||y!=ly) { // cursor has moved
				set_console_cursor(x, y); // set new cursor position
			}
			segment += s[shade];
			lx = x;
			ly = y;
		}
		print(segment, ltc, lbc); // print last segment, then reset color
#else // Linux
		const string s[3] = { "\u2591", "\u2592", "\u2593" };
		string r = ""; // append color changes to a string, print string in the end (much faster)
		uint lx=max_uint, ly=max_uint;
		int ltc=-1, lbc=-1;
		for(uint i=0u; i<changed_pixels; i++) {
			const uint n = (uint)print_video_lastframe->color(i); // get index of changed pixel
			const uint x=n%newwidth, y=n/newwidth; // get coordinates of changed pixel
			const int dither = print_image_rescaled->color(x, y); // get color of changed pixel
			const int shade = dither>>8;
			const int tc=(dither>>4)&0xF, bc=dither&0xF;
			if(x!=lx+1u||y!=ly) { // cursor has moved
				r += "\033["+to_string(y+1u)+";"+to_string(x+1u)+"f";
			}
			if(tc!=ltc||bc!=lbc||x!=lx+1u||y!=ly) { // color has changed or cursor has moved
				r += get_print_color(tc, bc); // print old segment, then change color
				ltc = tc;
				lbc = bc;
			}
			r += s[shade];
			lx = x;
			ly = y;
		}
		r += "\033[0m"; // reset color
		print(r);
#endif // Windows/Linux
	}
	set_console_cursor(0u, newheight);
	Image* swap = print_video_lastframe; // swap currentframe and lastframe
	print_video_lastframe = print_image_rescaled;
	print_image_rescaled = swap;
}
inline Image* screenshot(Image* image=nullptr) {
#if defined(_WIN32)
	HDC hdc = GetDC(NULL); // get the desktop device context
	HDC cdc = CreateCompatibleDC(hdc); // create a device context to use yourself
	const uint height = (uint)GetSystemMetrics(SM_CYVIRTUALSCREEN); // get the width and height of the screen
	const uint width  = 16u*height/9u; // only capture left monitor for dual screen setups, for both screens use (uint)GetSystemMetrics(SM_CXVIRTUALSCREEN);
	HBITMAP hbitmap = CreateCompatibleBitmap(hdc, width, height); // create a bitmap
	SelectObject(cdc, hbitmap); // use the previously created device context with the bitmap
	BITMAPINFOHEADER bmi = { 0 };
	bmi.biSize = sizeof(BITMAPINFOHEADER);
	bmi.biPlanes = 1;
	bmi.biBitCount = 32;
	bmi.biWidth = (int)width;
	bmi.biHeight = -(int)height; // flip image upright
	bmi.biCompression = BI_RGB;
	bmi.biSizeImage = 3*width*height;
	BitBlt(cdc, 0, 0, width, height, hdc, 0, 0, SRCCOPY); // copy from desktop device context to bitmap device context
	ReleaseDC(NULL, hdc);
	if(image==nullptr||image->width()!=width||image->height()!=height) {
		delete image;
		image = new Image(width, height);
	}
	GetDIBits(cdc, hbitmap, 0, height, image->data(), (BITMAPINFO*)&bmi, DIB_RGB_COLORS);
	DeleteObject(hbitmap);
	DeleteDC(cdc);
#endif // Windows
	return image;
}
inline void print_color_test() {
#ifdef _WIN32
	const string s = string("")+(char)223; // trick to double vertical resolution: use graphic character
#else // Linux
	const string s = "\u2580"; // trick to double vertical resolution: use graphic character
#endif // Windows/Linux
	print(s, color_magenta   , color_black     );
	print(s, color_blue      , color_gray      );
	print(s, color_light_blue, color_light_gray);
	print(s, color_cyan      , color_white     );
	print(s, color_green     , color_pink      );
	print(s, color_yellow    , color_dark_blue );
	print(s, color_orange    , color_dark_green);
	print(s, color_red       , color_dark_red  );
	println();
}
#ifdef UTILITIES_CONSOLE_INPUT
// ASCII codes (key>0): 8 backspace, 9 tab, 10 newline, 27 escape, 127 delete, !"#$%&'()*+,-./0-9:;<=>?@A-Z[]^_`a-z{|}~üäÄöÖÜßµ´§°¹³²
// control key codes (key<0): -38/-40/-37/-39 up/down/left/right arrow, -33/-34 page up/down, -36/-35 pos1/end
// other key codes (key<0): -45 insert, -144 num lock, -20 caps lock, -91 windows key, -93 kontext menu key, -112 to -123 F1 to F12
// not working: ¹ (251), num lock (-144), caps lock (-20), windows key (-91), kontext menu key (-93), F11 (-122)
#if defined(_WIN32)
inline int key_press() { // not working: F11 (-122, toggles fullscreen)
	KEY_EVENT_RECORD keyevent;
	INPUT_RECORD irec;
	DWORD events;
	while(true) {
		ReadConsoleInput(GetStdHandle(STD_INPUT_HANDLE), &irec, 1, &events);
		if(irec.EventType==KEY_EVENT&&((KEY_EVENT_RECORD&)irec.Event).bKeyDown) {
			keyevent = (KEY_EVENT_RECORD&)irec.Event;
			const int ca = (int)keyevent.uChar.AsciiChar;
			const int cv = (int)keyevent.wVirtualKeyCode;
			const int key = ca==0 ? -cv : ca+(ca>0?0:256);
			switch(key) {
				case  -16: continue; // disable Shift
				case  -17: continue; // disable Ctrl / AltGr
				case  -18: continue; // disable Alt / AltGr
				case -220: continue; // disable first detection of "^" key (not "^" symbol)
				case -221: continue; // disable first detection of "`" key (not "`" symbol)
				case -191: continue; // disable AltGr + "#"
				case  -52: continue; // disable AltGr + "4"
				case  -53: continue; // disable AltGr + "5"
				case  -54: continue; // disable AltGr + "6"
				case  -12: continue; // disable num block 5 with num lock deactivated
				case   13: return  10; // enter
				case  -46: return 127; // delete
				case  -49: return 251; // ¹
				case    0: continue;
				case    1: continue; // disable Ctrl + a (selects all text)
				case    2: continue; // disable Ctrl + b
				case    3: continue; // disable Ctrl + c (terminates program)
				case    4: continue; // disable Ctrl + d
				case    5: continue; // disable Ctrl + e
				case    6: continue; // disable Ctrl + f (opens search)
				case    7: continue; // disable Ctrl + g
				//case    8: continue; // disable Ctrl + h (ascii for backspace)
				//case    9: continue; // disable Ctrl + i (ascii for tab)
				case   10: continue; // disable Ctrl + j
				case   11: continue; // disable Ctrl + k
				case   12: continue; // disable Ctrl + l
				//case   13: continue; // disable Ctrl + m (breaks console, ascii for new line)
				case   14: continue; // disable Ctrl + n
				case   15: continue; // disable Ctrl + o
				case   16: continue; // disable Ctrl + p
				case   17: continue; // disable Ctrl + q
				case   18: continue; // disable Ctrl + r
				case   19: continue; // disable Ctrl + s
				case   20: continue; // disable Ctrl + t
				case   21: continue; // disable Ctrl + u
				case   22: continue; // disable Ctrl + v (inserts clipboard)
				case   23: continue; // disable Ctrl + w
				case   24: continue; // disable Ctrl + x
				case   25: continue; // disable Ctrl + y
				case   26: continue; // disable Ctrl + z
				default: return key; // any other ASCII/virtual character
			}
		}
	}
}
#else // Linux
#include <termios.h>
inline int key_press() { // not working: ¹ (251), num lock (-144), caps lock (-20), windows key (-91), kontext menu key (-93)
	struct termios term;
	tcgetattr(0, &term);
	while(true) {
		term.c_lflag &= ~(ICANON|ECHO); // turn off line buffering and echoing
		tcsetattr(0, TCSANOW, &term);
		int nbbytes;
		ioctl(0, FIONREAD, &nbbytes); // 0 is STDIN
		while(!nbbytes) {
			sleep(0.01);
			fflush(stdout);
			ioctl(0, FIONREAD, &nbbytes); // 0 is STDIN
		}
		int key = (int)getchar();
		if(key==27||key==194||key==195) { // escape, 194/195 is escape for °ß´äöüÄÖÜ
			key = (int)getchar();
			if(key==91) { // [ following escape
				key = (int)getchar(); // get code of next char after \e[
				if(key==49) { // F5-F8
					key = 62+(int)getchar(); // 53, 55-57
					if(key==115) key++; // F5 code is too low by 1
					getchar(); // take in following ~ (126), but discard code
				} else if(key==50) { // insert or F9-F12
					key = (int)getchar();
					if(key==126) { // insert
						key = 45;
					} else { // F9-F12
						key += 71; // 48, 49, 51, 52
						if(key<121) key++; // F11 and F12 are too low by 1
						getchar(); // take in following ~ (126), but discard code
					}
				} else if(key==51||key==53||key==54) { // delete, page up/down
					getchar(); // take in following ~ (126), but discard code
				}
			} else if(key==79) { // F1-F4
				key = 32+(int)getchar(); // 80-83
			}
			key = -key; // use negative numbers for escaped keys
		}
		term.c_lflag |= (ICANON|ECHO); // turn on line buffering and echoing
		tcsetattr(0, TCSANOW, &term);
		switch(key) {
			case  127: return   8; // backspace
			case  -27: return  27; // escape
			case  -51: return 127; // delete
			case -164: return 132; // ä
			case -182: return 148; // ö
			case -188: return 129; // ü
			case -132: return 142; // Ä
			case -150: return 153; // Ö
			case -156: return 154; // Ü
			case -159: return 225; // ß
			case -181: return 230; // µ
			case -167: return 245; // §
			case -176: return 248; // °
			case -178: return 253; // ²
			case -179: return 252; // ³
			case -180: return 239; // ´
			case  -65: return -38; // up arrow
			case  -66: return -40; // down arrow
			case  -68: return -37; // left arrow
			case  -67: return -39; // right arrow
			case  -53: return -33; // page up
			case  -54: return -34; // page down
			case  -72: return -36; // pos1
			case  -70: return -35; // end
			case    0: continue;
			case    1: continue; // disable Ctrl + a
			case    2: continue; // disable Ctrl + b
			case    3: continue; // disable Ctrl + c (terminates program)
			case    4: continue; // disable Ctrl + d
			case    5: continue; // disable Ctrl + e
			case    6: continue; // disable Ctrl + f
			case    7: continue; // disable Ctrl + g
			case    8: continue; // disable Ctrl + h
			//case    9: continue; // disable Ctrl + i (ascii for tab)
			//case   10: continue; // disable Ctrl + j (ascii for new line)
			case   11: continue; // disable Ctrl + k
			case   12: continue; // disable Ctrl + l
			case   13: continue; // disable Ctrl + m
			case   14: continue; // disable Ctrl + n
			case   15: continue; // disable Ctrl + o
			case   16: continue; // disable Ctrl + p
			case   17: continue; // disable Ctrl + q
			case   18: continue; // disable Ctrl + r
			case   19: continue; // disable Ctrl + s
			case   20: continue; // disable Ctrl + t
			case   21: continue; // disable Ctrl + u
			case   22: continue; // disable Ctrl + v
			case   23: continue; // disable Ctrl + w
			case   24: continue; // disable Ctrl + x
			case   25: continue; // disable Ctrl + y
			case   26: continue; // disable Ctrl + z (terminates program)
			default: return key; // any other ASCII character
		}
	}
}
#endif // Windows/Linux
#endif // UTILITIES_CONSOLE_INPUT
#else // UTILITIES_CONSOLE_COLOR
inline void print(const string& s, const int textcolor) {
	std::cout << s;
}
inline void print(const string& s, const int textcolor, const int backgroundcolor) {
	std::cout << s;
}
#endif // UTILITIES_CONSOLE_COLOR

#ifdef UTILITIES_REGEX
inline vector<string> split_regex(const string& s, const string& separator="\\s+") {
	vector<string> r;
	const std::regex rgx(separator);
	std::sregex_token_iterator token(s.begin(), s.end()+1, rgx, -1), end;
	while(token!=end) {
		r.push_back(*token);
		token++;
	}
	return r;
}
inline bool equals_regex(const string& s, const string& match) { // returns true if string exactly matches regex
	return regex_match(s.begin(), s.end(), std::regex(match));
}
inline uint matches_regex(const string& s, const string& match) { // counts number of matches
	std::regex words_regex(match);
	auto words_begin = std::sregex_iterator(s.begin(), s.end(), words_regex);
	auto words_end = std::sregex_iterator();
	return (uint)std::distance(words_begin, words_end);
}
inline bool contains_regex(const string& s, const string& match) {
	return matches_regex(s, match)>=1;
}
inline string replace_regex(const string& s, const string& from, const string& to) {
	return regex_replace(s, std::regex(from), to);
}
inline bool is_number(const string& s) {
	return equals_regex(s, "\\d+(u|l|ul|ll|ull)?")||equals_regex(s, "0x(\\d|[a-fA-F])+(u|l|ul|ll|ull)?")||equals_regex(s, "0b[01]+(u|l|ul|ll|ull)?")||equals_regex(s, "(((\\d+\\.?\\d*|\\.\\d+)([eE][+-]?\\d+[fF]?)?)|(\\d+\\.\\d*|\\.\\d+)[fF]?)");
}
inline void print_message(const string& message, const string& keyword="", const int keyword_color=-1, const int colons=true) { // print formatted message
	const uint k=length(keyword)+2u, w=CONSOLE_WIDTH-4u-k;
	string p=colons?": ":"  ", f="";
	for(uint j=0u; j<k; j++) f += " ";
	vector<string> v = split_regex(message, "[\\s\\0]+");
	uint l = 0u; // length of current line of words
	for(uint i=0u; i<(uint)v.size(); i++) {
		const string word = v.at(i);
		const uint wordlength = length(word);
		l += wordlength+1u; // word + space
		if(l<=w+1u) { // word fits -> append word and space
			p += word+" ";
		} else if(wordlength>w) { // word overflows -> split word into next line
			p += substring(word, 0, w-(l-wordlength-1u))+" |\n| "+f;
			v[i] = substring(v[i], w-(l-wordlength-1u)); i--; // reuse same vector element for overflowing part, decrement i to start next line with this overflowing part
			l = 0u; // reset line length
		} else { // word does not fit -> fill remaining line with spaces
			l = l-length(v.at(i--))-1u; // remove word from line, decrement i to start next line with this word
			for(uint j=l; j<=w; j++) p += " ";
			p += "|\n| "+f;
			l = 0u; // reset line length
		}
	}
	for(uint j=l; j<=w; j++) p += " ";
	if(keyword_color<0||keyword_color>15) {
		print("\r| "+keyword);
	} else {
		print("\r| ");
		print(keyword, keyword_color);
	}
	println(p+"|");
}
inline void print_error(const string& s) { // print formatted error message
	print_message(s, "Error", color_red);
#ifdef _WIN32
	print_message("Press Enter to exit.", "     ", -1, false);
#endif // _WIN32
	string b = "";
	for(int i=0; i<CONSOLE_WIDTH-2; i++) b += "-";
	println("'"+b+"'");
#ifdef _WIN32
	wait();
#endif //_WIN32
	exit(1);
}
inline void print_warning(const string& s) { // print formatted warning message
	print_message(s, "Warning", color_orange);
}
inline void print_info(const string& s) { // print formatted info message
	print_message(s, "Info", color_green);
}

inline void parse_sanity_check_error(const string& s, const string& regex, const string& type) {
	if(!equals_regex(s, regex)) print_error("\""+s+"\" cannot be parsed to "+type+".");
}
inline int to_int(const string& s) {
	const string t = trim(s);
	parse_sanity_check_error(t, "[+-]?\\d+", "int");
	return atoi(t.c_str());
}
inline uint to_uint(const string& s) {
	const string t = trim(s);
	parse_sanity_check_error(t, "\\+?\\d+", "uint");
	return (uint)atoi(t.c_str());
}
inline slong to_slong(const string& s) {
	const string t = trim(s);
	parse_sanity_check_error(t, "[+-]?\\d+", "slong");
	return (slong)atoll(t.c_str());
}
inline ulong to_ulong(const string& s) {
	const string t = trim(s);
	parse_sanity_check_error(t, "\\+?\\d+", "ulong");
	return (ulong)atoll(t.c_str());
}
inline float to_float(const string& s) {
	const string t = trim(s);
	parse_sanity_check_error(t, "[+-]?(((\\d+\\.?\\d*|\\.\\d+)([eE][+-]?\\d+[fF]?)?)|(\\d+\\.\\d*|\\.\\d+)[fF]?)", "float");
	return (float)atof(t.c_str());
}
inline double to_double(const string& s) {
	const string t = trim(s);
	parse_sanity_check_error(t, "[+-]?(((\\d+\\.?\\d*|\\.\\d+)([eE][+-]?\\d+[fF]?)?)|(\\d+\\.\\d*|\\.\\d+)[fF]?)", "double");
	return atof(t.c_str());
}

inline bool parse_sanity_check(const string& s, const string& regex) {
	return equals_regex(s, regex);
}
inline int to_int(const string& s, const int default_value) {
	const string t = trim(s);
	return parse_sanity_check(t, "[+-]?\\d+") ? atoi(t.c_str()) : default_value;
}
inline uint to_uint(const string& s, const uint default_value) {
	const string t = trim(s);
	return parse_sanity_check(t, "\\+?\\d+") ? (uint)atoi(t.c_str()) : default_value;
}
inline slong to_slong(const string& s, const slong default_value) {
	const string t = trim(s);
	return parse_sanity_check(t, "[+-]?\\d+") ? (slong)atoll(t.c_str()) : default_value;
}
inline ulong to_ulong(const string& s, const ulong default_value) {
	const string t = trim(s);
	return parse_sanity_check(t, "\\+?\\d+") ? (ulong)atoll(t.c_str()) : default_value;
}
inline float to_float(const string& s, const float default_value) {
	const string t = trim(s);
	return parse_sanity_check(t, "[+-]?(((\\d+\\.?\\d*|\\.\\d+)([eE][+-]?\\d+[fF]?)?)|(\\d+\\.\\d*|\\.\\d+)[fF]?)") ? (float)atof(t.c_str()) : default_value;
}
inline double to_double(const string& s, const double default_value) {
	const string t = trim(s);
	return parse_sanity_check(t, "[+-]?(((\\d+\\.?\\d*|\\.\\d+)([eE][+-]?\\d+[fF]?)?)|(\\d+\\.\\d*|\\.\\d+)[fF]?)") ? atof(t.c_str()) : default_value;
}
#else // UTILITIES_REGEX
inline void print_message(const string& message, const string& keyword="", const int keyword_color=-1, const int colons=true) { // print message
	println(keyword+": "+message);
}
inline void print_error(const string& s) { // print error message
	println("Error: "+s);
#ifdef _WIN32
	println("       Press Enter to exit.");
	wait();
#endif //_WIN32
	exit(1);
}
inline void print_warning(const string& s) { // print warning message
	println("Warning: "+s);
}
inline void print_info(const string& s) { // print info message
	println("Info: "+s);
}
#endif // UTILITIES_REGEX

#ifdef UTILITIES_FILE
#include <fstream> // read/write files
#ifndef UTILITIES_NO_CPP17
#include <filesystem> // automatically create directory before writing file, requires C++17
inline vector<string> find_files(const string& path, const string& extension=".*") {
	vector<string> files;
	if(std::filesystem::is_directory(path)&&std::filesystem::exists(path)) {
		for(const auto& entry : std::filesystem::directory_iterator(path)) {
			if(extension==".*"||entry.path().extension().string()==extension) files.push_back(entry.path().string());
		}
	}
	return files;
}
#endif // UTILITIES_NO_CPP17
inline void create_folder(const string& path) { // create folder if it not already exists
	const int slash_position = (int)path.rfind('/'); // find last slash dividing the path from the filename
	if(slash_position==(int)string::npos) return; // no slash found
	const string f = path.substr(0, slash_position); // cut off file name if there is any
#ifndef UTILITIES_NO_CPP17
	if(!std::filesystem::is_directory(f)||!std::filesystem::exists(f)) std::filesystem::create_directories(f); // create folder if it not already exists
#endif // UTILITIES_NO_CPP17
}
inline string create_file_extension(const string& filename, const string& extension) {
	return filename.substr(0, filename.rfind('.'))+(extension.at(0)!='.'?".":"")+extension; // remove existing file extension if existing and replace it with new one
}
inline string read_file(const string& filename) {
	std::ifstream file(filename, std::ios::in);
	if(file.fail()) print_error("File \""+filename+"\" does not exist!");
	const string r((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
	file.close();
	return r;
}
inline void write_file(const string& filename, const string& content="") {
	create_folder(filename);
	std::ofstream file(filename, std::ios::out);
	file.write(content.c_str(), content.length());
	file.close();
}
inline void write_line(const string& filename, const string& content="") {
	const string s = content+(ends_with(content, "\n")?"":"\n"); // add new line if needed
	std::ofstream file(filename, std::ios::ios_base::app);
	file.write(s.c_str(), s.length());
	file.close();
}
template<typename T> inline void write_file(const string& filename, const string& header, const uint n, const T* y) {
	string s = header;
	if(length(s)>0u && !ends_with(s, "\n")) s += "\n";
	for(uint i=0u; i<n; i++) s += to_string(i)+"\t"+to_string(y[i])+"\n";
	write_file(filename, s);
}
template<typename T, typename U> inline void write_file(const string& filename, const string& header, const uint n, const T* x, const U* y) {
	string s = header;
	if(length(s)>0u && !ends_with(s, "\n")) s += "\n";
	for(uint i=0u; i<n; i++) s += to_string(x[i])+"\t"+to_string(y[i])+"\n";
	write_file(filename, s);
}
#pragma warning(disable:6385)
inline Image* read_bmp(const string& filename, Image* image=nullptr) {
	std::ifstream file(create_file_extension(filename, ".bmp"), std::ios::in|std::ios::binary);
	if(file.fail()) print_error("File \""+filename+"\" does not exist!");
	uint width=0u, height=0u;
	file.seekg(0, std::ios::end);
	const uint filesize = (uint)file.tellg();
	file.seekg(0, std::ios::beg);
	uchar* data = new uchar[filesize];
	file.read((char*)data, filesize);
	file.close();
	if(filesize==0u) print_error("File \""+filename+"\" is corrupt!");
	for(uint i=0u; i<4u; i++) {
		width  |= data[18+i]<<(8u*i);
		height |= data[22+i]<<(8u*i);
	}
	const uint pad=(4u-(3u*width)%4u)%4u, imagesize=(3u*width+pad)*height;
	if(filesize!=54u+imagesize) print_error("File \""+filename+"\" is corrupt or unsupported!");
	if(image==nullptr||image->width()!=width||image->height()!=height) {
		delete image;
		image = new Image(width, height);
	}
	for(uint y=0u; y<height; y++) {
		for(uint x=0u; x<width; x++) {
			const uint i = 54+3u*x+y*(3u*width+pad);
			image->set_color(x, height-1-y, ((int)data[i]&255)|((int)data[i+1u]&255)<<8|((int)data[i+2u]&255)<<16);
		}
	}
	delete[] data;
	return image;
}
inline void write_bmp(const string& filename, const Image* image) {
	create_folder(filename);
	const uint pad=(4u-(3u*image->width())%4u)%4u, filesize=54u+(3u*image->width()+pad)*image->height(); // horizontal line must be a multiple of 4 bytes long, header is 54 bytes
	uchar header[54] = { 'B','M', 0,0,0,0, 0,0,0,0, 54,0,0,0, 40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0,24,0 };
	for(uint i=0u; i<4u; i++) {
		header[ 2u+i] = (uchar)((filesize       >>(8u*i))&255);
		header[18u+i] = (uchar)((image->width() >>(8u*i))&255);
		header[22u+i] = (uchar)((image->height()>>(8u*i))&255);
	}
	uchar* data = new uchar[filesize];
	for(uint i=0u; i<54u; i++) data[i] = header[i];
	for(uint y=0u; y<image->height(); y++) {
		for(uint x=0u; x<image->width(); x++) {
			const int color = image->color(x, image->height()-1u-y);
			const uint i = 54u+3u*x+y*(3u*image->width()+pad);
			data[i   ] = (uchar)( color     &255);
			data[i+1u] = (uchar)((color>> 8)&255);
			data[i+2u] = (uchar)((color>>16)&255);
		}
		for(uint p=0u; p<pad; p++) data[54+(3u*image->width()+p)+y*(3u*image->width()+pad)] = 0;
	}
	std::ofstream file(create_file_extension(filename, ".bmp"), std::ios::out|std::ios::binary);
	file.write((char*)data, filesize);
	file.close();
	delete[] data;
}
inline Image* read_qoi(const string& filename, Image* image=nullptr) { // 4-channel .qoi decoder, source: https://qoiformat.org/qoi-specification.pdf
	std::ifstream file(create_file_extension(filename, ".qoi"), std::ios::in|std::ios::binary);
	if(file.fail()) print_error("File \""+filename+"\" does not exist!");
	file.seekg(0, std::ios::end);
	const uint filesize = (uint)file.tellg();
	file.seekg(0, std::ios::beg);
	uchar* data = new uchar[filesize];
	file.read((char*)data, filesize);
	file.close();
	if(filesize<22u||((int*)data)[0]!=1718185841) print_error("File \""+filename+"\" is corrupt!");
	uint width=0u, height=0u;
	for(uint i=0u; i<4u; i++) {
		width  |= data[ 7u-i]<<(8u*i);
		height |= data[11u-i]<<(8u*i);
	}
	if(image==nullptr||image->width()!=width||image->height()!=height) {
		delete image;
		image = new Image(width, height);
	}
	uint counter = 14u; // header size
	int pixel=0xFF000000, runlength=0, lookup[64]={0};
	for(uint i=0u; i<image->length()&&counter<filesize-8u; i++) {
		if(runlength>0) {
			runlength--;
		} else {
			const int chunk = (int)data[counter++];
			const int tag = chunk&0b11000000;
			if(chunk==0b11111110) { // raw rgb encoding (4 bytes)
				const int r = (int)data[counter++];
				const int g = (int)data[counter++];
				const int b = (int)data[counter++];
				pixel = color(r, g, b, 255);
			} else if(chunk==0b11111111) { // raw rgba encoding (5 bytes)
				const int r = (int)data[counter++];
				const int g = (int)data[counter++];
				const int b = (int)data[counter++];
				const int a = (int)data[counter++];
				pixel = color(r, g, b, a);
			} else if(tag==0b00000000) { // lookup table encoding (1 byte)
				pixel = lookup[chunk]; // chunk&0b00111111
			} else if(tag==0b01000000) { // difference encoding (1 byte)
				const int r = red  (pixel)+((chunk>>4)&3)-2;
				const int g = green(pixel)+((chunk>>2)&3)-2;
				const int b = blue (pixel)+( chunk    &3)-2;
				pixel = color(r, g, b, alpha(pixel));
			} else if(tag==0b10000000) { // luma encoding (2 bytes)
				const int luma = (int)data[counter++];
				const int dg=(chunk&0b00111111)-32, drdg=((luma>>4)&0xF)-8, dbdg=(luma&0xF)-8;
				const int r = red  (pixel)+drdg+dg;
				const int g = green(pixel)+     dg;
				const int b = blue (pixel)+dbdg+dg;
				pixel = color(r, g, b, alpha(pixel));
			} else if(tag==0b11000000) { // runlength encoding (1 byte)
				runlength = chunk&0b00111111;
			}
			lookup[(3*red(pixel)+5*green(pixel)+7*blue(pixel)+11*alpha(pixel))%64] = pixel;
		}
		image->set_color(i, pixel);
	}
	delete[] data;
	return image;
}
inline void write_qoi(const string& filename, const Image* image) { // 3-channel .qoi encoder, source: https://qoiformat.org/qoi-specification.pdf
	create_folder(filename);
	uchar header[14] = { 'q','o','i','f', 0,0,0,0, 0,0,0,0, 3,0 }; // 3 channels for red, green, blue (no alpha)
	const uchar padding[8] = { 0,0,0,0,0,0,0,1 };
	for(uint i=0u; i<4u; i++) {
		header[ 7u-i] = (char)((image->width() >>(8u*i))&255);
		header[11u-i] = (char)((image->height()>>(8u*i))&255);
	}
	int current=0, previous=current, runlength=0, lookup[64]={0};
	uchar* data = new uchar[22u+4u*image->length()];
	uint filesize = 0u;
	for(uint i=0u; i<14u; i++) data[filesize++] = header[i];
	for(uint i=0u; i<image->length(); i++) {
		current = image->color(i)&0x00FFFFFF; // ignore alpha channel
		if(current==previous) { // runlength encoding (1 byte)
			runlength++;
			if(runlength>=62||i>=image->length()-1u) {
				data[filesize++] = (uchar)(0b11000000|(runlength-1));
				runlength = 0;
			}
		} else {
			if(runlength>0) { // runlength encoding (1 byte)
				data[filesize++] = (uchar)(0b11000000|(runlength-1));
				runlength = 0;
			}
			const int hash = (3*red(current)+5*green(current)+7*blue(current)+11*255)%64;
			if(lookup[hash]==current) { // lookup table encoding (1 byte)
				data[filesize++] = (uchar)hash; // 0b00000000|hash
			} else {
				lookup[hash] = current;
				const int dr=red(current)-red(previous), dg=green(current)-green(previous), db=blue(current)-blue(previous);
				const int drdg=dr-dg, dbdg=db-dg;
				if(dr>=-2&&dr<2 && dg>=-2&&dg<2 && db>=-2&&db<2) { // difference encoding (1 byte)
					data[filesize++] = (uchar)(0b01000000|(dr+2)<<4|(dg+2)<<2|(db+2));
				} else if(dg>=-32&&dg<32 && drdg>=-8&&drdg<8 && dbdg>=-8&&dbdg<8) { // luma encoding (2 bytes)
					data[filesize++] = (uchar)(0b10000000|(dg+32));
					data[filesize++] = (uchar)((drdg+8)<<4|(dbdg+8));
				} else {
					data[filesize++] = (uchar)0b11111110; // raw rgb encoding (4 bytes)
					data[filesize++] = (uchar)red(current);
					data[filesize++] = (uchar)green(current);
					data[filesize++] = (uchar)blue(current);
				}
			}
		}
		previous = current;
	}
	for(uint i=0u; i<8u; i++) data[filesize++] = padding[i];
	std::ofstream file(create_file_extension(filename, ".qoi"), std::ios::out|std::ios::binary);
	file.write((char*)data, filesize);
	file.close();
	delete[] data;
}
#ifdef UTILITIES_PNG
#include "lodepng.hpp"
inline Image* read_png(const string& filename, Image* image=nullptr) {
	uint width=0u, height=0u;
	vector<uchar> data;
	lodepng::decode(data, width, height, create_file_extension(filename, ".png"), LCT_RGB);
	if(image==nullptr||image->width()!=width||image->height()!=height) {
		delete image;
		image = new Image(width, height);
	}
	for(uint i=0u; i<width*height; i++) {
		image->set_color(i, data[3u*i]<<16|data[3u*i+1u]<<8|data[3u*i+2u]);
	}
	return image;
}
inline void write_png(const string& filename, const Image* image) {
	create_folder(filename);
	uchar* data = new uchar[3u*image->length()];
	for(uint i=0u; i<image->length(); i++) {
		const int color = image->color(i);
		data[3u*i   ] = (color>>16)&255;
		data[3u*i+1u] = (color>> 8)&255;
		data[3u*i+2u] =  color     &255;
	}
	lodepng::encode(create_file_extension(filename, ".png"), data, image->width(), image->height(), LCT_RGB);
	delete[] data;
}
#endif // UTILITIES_PNG

struct Mesh { // triangle mesh
	uint triangle_number = 0u;
	float3 center, pmin, pmax;
	float3* p0;
	float3* p1;
	float3* p2;
	inline Mesh(const uint triangle_number, const float3& center) {
		this->triangle_number = triangle_number;
		this->center = this->pmin = this->pmax = center;
		this->p0 = new float3[triangle_number];
		this->p1 = new float3[triangle_number];
		this->p2 = new float3[triangle_number];
	}
	inline ~Mesh() {
		delete[] p0;
		delete[] p1;
		delete[] p2;
	}
	inline void find_bounds() {
		pmin = pmax = p0[0];
		for(uint i=1u; i<triangle_number; i++) {
			const float3 p0i=p0[i], p1i=p1[i], p2i=p2[i];
			pmin.x = fmin(fmin(fmin(p0i.x, p1i.x), p2i.x), pmin.x);
			pmin.y = fmin(fmin(fmin(p0i.y, p1i.y), p2i.y), pmin.y);
			pmin.z = fmin(fmin(fmin(p0i.z, p1i.z), p2i.z), pmin.z);
			pmax.x = fmax(fmax(fmax(p0i.x, p1i.x), p2i.x), pmax.x);
			pmax.y = fmax(fmax(fmax(p0i.y, p1i.y), p2i.y), pmax.y);
			pmax.z = fmax(fmax(fmax(p0i.z, p1i.z), p2i.z), pmax.z);
		}
	}
	inline void scale(const float scale) {
		for(uint i=0u; i<triangle_number; i++) {
			p0[i] = scale*(p0[i]-center)+center;
			p1[i] = scale*(p1[i]-center)+center;
			p2[i] = scale*(p2[i]-center)+center;
		}
		pmin = scale*(pmin-center)+center;
		pmax = scale*(pmax-center)+center;
	}
	inline void translate(const float3& translation) {
		for(uint i=0u; i<triangle_number; i++) {
			p0[i] += translation;
			p1[i] += translation;
			p2[i] += translation;
		}
		center += translation;
		pmin += translation;
		pmax += translation;
	}
	inline void rotate(const float3x3& rotation) {
		for(uint i=0u; i<triangle_number; i++) {
			p0[i] = rotation*(p0[i]-center)+center;
			p1[i] = rotation*(p1[i]-center)+center;
			p2[i] = rotation*(p2[i]-center)+center;
		}
		find_bounds();
	}
	inline void set_center(const float3& center) {
		this->center = center;
	}
	inline const float3& get_center() const {
		return center;
	}
	inline const float3 get_bounding_box_size() const {
		return pmax-pmin;
	}
	inline const float3 get_bounding_box_center() const {
		return 0.5f*(pmin+pmax);
	}
	inline float get_min_size() const {
		return fmin(fmin(pmax.x-pmin.x, pmax.y-pmin.y), pmax.z-pmin.z);
	}
	inline float get_max_size() const {
		return fmax(fmax(pmax.x-pmin.x, pmax.y-pmin.y), pmax.z-pmin.z);
	}
	inline float get_scale_for_box_fit(const float3& box_size) const { // scale factor to exactly fit mesh bounding box in simulation box
		return fmin(fmin(box_size.x/(pmax.x-pmin.x), box_size.y/(pmax.y-pmin.y)), box_size.z/(pmax.z-pmin.z));
	}
};
inline Mesh* read_stl_raw(const string& path, const bool reposition, const float3& box_size, const float3& center, const float3x3& rotation, const float size) { // read binary .stl file
	const string filename = create_file_extension(path, ".stl");
	std::ifstream file(filename, std::ios::in|std::ios::binary);
	if(file.fail()) print_error("File \""+filename+"\" does not exist!");
	file.seekg(0, std::ios::end);
	const uint filesize = (uint)file.tellg();
	file.seekg(0, std::ios::beg);
	uchar* data = new uchar[filesize];
	file.read((char*)data, filesize);
	file.close();
	if(filesize==0u) print_error("File \""+filename+"\" is corrupt!");
	const uint triangle_number = ((uint*)data)[20];
	uint counter = 84u;
	if(triangle_number>0u&&filesize==84u+50u*triangle_number) print_info("Loading \""+filename+"\" with "+to_string(triangle_number)+" triangles.");
	else print_error("File \""+filename+"\" is corrupt or unsupported! Only binary .stl files are supported.");
	Mesh* mesh = new Mesh(triangle_number, center);
	mesh->p0[0] = float3(0.0f); // to fix warning C6001
	for(uint i=0u; i<triangle_number; i++) {
		const float* triangle_data = (float*)(data+counter);
		counter += 50u;
		mesh->p0[i] = rotation*float3(triangle_data[ 3], triangle_data[ 4], triangle_data[ 5]); // read positions of triangle vertices and rotate them
		mesh->p1[i] = rotation*float3(triangle_data[ 6], triangle_data[ 7], triangle_data[ 8]);
		mesh->p2[i] = rotation*float3(triangle_data[ 9], triangle_data[10], triangle_data[11]);
	}
	delete[] data;
	mesh->find_bounds();
	float scale = 1.0f;
	if(size==0.0f) { // auto-rescale to largest possible size
		scale = mesh->get_scale_for_box_fit(box_size);
	} else if(size>0.0f) { // rescale longest bounding box side length of mesh to specified size
		scale = size/mesh->get_max_size();
	} else { // rescale to specified size relative to original size (input size as negative number)
		scale = -size;
	}
	const float3 offset = reposition ? -0.5f*(mesh->pmin+mesh->pmax) : float3(0.0f); // auto-reposition mesh
	for(uint i=0u; i<triangle_number; i++) { // rescale mesh
		mesh->p0[i] = center+scale*(offset+mesh->p0[i]);
		mesh->p1[i] = center+scale*(offset+mesh->p1[i]);
		mesh->p2[i] = center+scale*(offset+mesh->p2[i]);
	}
	mesh->find_bounds();
	return mesh;
}
inline Mesh* read_stl(const string& path, const float3& box_size, const float3& center, const float3x3& rotation, const float size) { // read binary .stl file (rescale and reposition)
	return read_stl_raw(path, true, box_size, center, rotation, size);
}
inline Mesh* read_stl(const string& path, const float3& box_size, const float3& center, const float size) { // read binary .stl file (rescale and reposition, no rotation)
	return read_stl_raw(path, true, box_size, center, float3x3(1.0f), size);
}
inline Mesh* read_stl(const string& path, const float scale=1.0f, const float3x3& rotation=float3x3(1.0f), const float3& offset=float3(0.0f)) { // read binary .stl file (do not auto-rescale and auto-reposition)
	return read_stl_raw(path, false, float3(1.0f), offset, rotation, -fabs(scale));
}

class Configuration_File {
private:
	string filename;
	struct Entry {
		string name;
		string value;
		inline Entry(const string& name, const string& value) {
			this->name = name;
			this->value = value;
		}
	};
	vector<Entry> entries;
	inline vector<string> split_vector(string values) const {
		if((values.front()=='{'&&values.back()=='}') || (values.front()=='['&&values.back()==']') || (values.front()=='('&&values.back()==')') || (values.front()=='<'&&values.back()=='>')) values = values.substr(1, values.length()-2);
		vector<string> elements = split_regex(values, "\\s*[,;]\\s*");
		for(uint j=0u; j<(uint)elements.size(); j++) elements[j] = trim(elements[j]);
		return elements;
	}
	template<class T> inline T from_string(const string& s, const string& name) const {
		if constexpr(std::is_same<T, bool>::value) return contains(to_lower(s),"true") ? (T)1 : contains(to_lower(s),"false") ? (T)0 : (T)to_int(s);
		else if constexpr(std::is_same<T,  char >::value) return (T)to_int   (s);
		else if constexpr(std::is_same<T, uchar >::value) return (T)to_uint  (s);
		else if constexpr(std::is_same<T,  short>::value) return (T)to_int   (s);
		else if constexpr(std::is_same<T, ushort>::value) return (T)to_uint  (s);
		else if constexpr(std::is_same<T,  int  >::value) return (T)to_int   (s);
		else if constexpr(std::is_same<T, uint  >::value) return (T)to_uint  (s);
		else if constexpr(std::is_same<T, slong >::value) return (T)to_slong (s);
		else if constexpr(std::is_same<T, ulong >::value) return (T)to_ulong (s);
		else if constexpr(std::is_same<T, float >::value) return (T)to_float (s);
		else if constexpr(std::is_same<T, double>::value) return (T)to_double(s);
		else if constexpr(std::is_same<T, string>::value) return (s.front()=='"'&&s.back()=='"') ? s.substr(1, s.length()-2) : s;
		else print_error("Generic type in call value<T>(\""+name+"\") is not supported.");
		return T(); // otherwise compiler complains about missing return statement
	}
	template<typename T> inline bool extract_value(const string& name, T& v) const {
		for(uint i=0u; i<(uint)entries.size(); i++) {
			if(entries[i].name==name) {
				v = from_string<T>(entries[i].value, name);
				return true;
			}
		}
		return false;
	}
	template<typename T> inline bool extract_value(const string& name, vector<T>& v) const {
		for(uint i=0u; i<(uint)entries.size(); i++) {
			if(entries[i].name==name) {
				const vector<string> elements = split_vector(entries[i].value);
				for(uint j=0u; j<(uint)elements.size(); j++) v.push_back(from_string<T>(elements[j], name));
				return true;
			}
		}
		return false;
	}
public:
	inline Configuration_File(const string& filename) {
		this->filename = string(filename);
		const string file = read_file(filename);
		vector<string> lines = split_regex(file, "\\s*\\n\\s*");
		for(uint i=0u; i<(uint)lines.size(); i++) {
			lines[i] = trim(lines[i]); // remove whitespace
			if(length(lines[i])==0u || begins_with(lines[i], "//") || lines[i].front()=='#' || (lines[i].front()=='['&&lines[i].back()==']')) lines.erase(lines.begin()+i); // remove comment lines
		}
		for(uint i=0u; i<(uint)lines.size(); i++) {
			const vector<string> data = split_regex(lines[i], "\\s*=\\s*"); // split lines in "name = value"
			if((int)data.size()!=2) print_error("Entry \""+lines[i]+"\" in configuration file \""+filename+"\" is invalid.");
			entries.push_back(Entry(trim(data[0]), trim(data[1])));
		}
	}
	template<class T> inline T value(const string& name) const { // valid for bool, char, uchar, short, ushort, int, uint, slong, ulong, float, double, string and vectors of these types
		T t;
		if(extract_value(name, t)) return t;
		else print_error("There is no entry called \""+name+"\" in configuration file \""+filename+"\".");
	}
	template<class T> inline T value(const string& name, const T& default_value) const { // valid for bool, char, uchar, short, ushort, int, uint, slong, ulong, float, double, string and vectors of these types
		T t;
		if(extract_value(name, t)) return t;
		else return default_value;
	}
	inline void print_entries() const {
		for(uint i=0u; i<(uint)entries.size(); i++) {
			println(entries[i].name+" = "+entries[i].value);
		}
	}
};
#endif // UTILITIES_FILE