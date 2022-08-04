#pragma once

#include "utilities.hpp"
#define R(...) string(" "#__VA_ARGS__" ") // evil stringification macro, similar syntax to raw string R"(...)"

string opencl_c_container(); // outsourced to kernel.cpp
string get_opencl_c_code() {
	string r = opencl_c_container();
	r = replace(r, " ", "\n"); // replace all spaces by new lines
	r = replace(r, "#ifdef\n", "#ifdef "); // except for the arguments after some preprocessor options that need to be in the same line
	r = replace(r, "#ifndef\n", "#ifndef ");
	r = replace(r, "#define\n", "#define "); // #define with two arguments will not work
	r = replace(r, "#if\n", "#if "); // don't leave any spaces in arguments
	r = replace(r, "#elif\n", "#elif "); // don't leave any spaces in arguments
	r = replace(r, "#pragma\n", "#pragma ");
	return "\n"+r;
}

// everything below is just for syntax highlighting in the editor, this does not change any functionality
// full catalogue: https://www.khronos.org/files/opencl-1-2-quick-reference-card.pdf

// general
#define get_global_id(x) // global index, set x=0
#define get_global_size(x) // global range, set x=0
#define get_local_id(x) // local index within group, set x=0
#define get_local_size(x) // group size, set x=0
#define	get_num_groups(x) // number of groups, set x=0
#define get_group_id(x) // group ID, set x=0
#define get_global_offset(x) // global offset, set x=0
#define get_work_dim // number of dimensions in use
#define __attribute__(x)
#define always_inline // for inlining functions
#define barrier(x) // barrier for local work group, x is CLK_LOCAL_MEM_FENCE or CLK_GLOBAL_MEM_FENCE
#define mem_fence(x) // orders loads/stores, x is CLK_LOCAL_MEM_FENCE or CLK_GLOBAL_MEM_FENCE
#define read_mem_fence(x) // orders loads, x is CLK_LOCAL_MEM_FENCE or CLK_GLOBAL_MEM_FENCE
#define write_mem_fence(x) // orders stores, x is CLK_LOCAL_MEM_FENCE or CLK_GLOBAL_MEM_FENCE
#define CLK_LOCAL_MEM_FENCE
#define CLK_GLOBAL_MEM_FENCE
#define kernel
#define constant
#define global
#define local
#define __kernel
#define __constant
#define __global
#define __local
#define __private // private keyword already exists in C++

// 32-bit integer atomics
#define atomic_add(p,x) // (*p)+=x
#define atomic_sub(p,x) // (*p)-=x
#define atomic_xchg(p,x) // t=(*p);(*p)=x;x=t;
#define atomic_inc(p) // (*p)+=1
#define atomic_dec(p) // (*p)-=1
#define atomic_cmpxchg(p,c,x) // (*p)=((*p)==c?x:(*p))
#define atomic_max(p,x) // (*p)=max(*p,x)
#define atomic_min(p,x) // (*p)=min(*p,x)
#define atomic_and(p,x) // (*p)=(*p)&x
#define atomic_or(p,x) // (*p)=(*p)|x
#define atomic_xor(p,x) // (*p)=(*p)^x

// 64-bit integer atomics (cl_khr_int64_base_atomics extension must be supported by the device)
#define atom_add(p,x) // (*p)+=x
#define atom_sub(p,x) // (*p)-=x
#define atom_xchg(p,x) // t=(*p);(*p)=x;x=t;
#define atom_inc(p) // (*p)+=1
#define atom_dec(p) // (*p)-=1
#define atom_cmpxchg(p,c,x) // (*p)=((*p)==c?x:(*p))
#define atom_max(p,x) // (*p)=max(*p,x)
#define atom_min(p,x) // (*p)=min(*p,x)
#define atom_and(p,x) // (*p)=(*p)&x
#define atom_or(p,x) // (*p)=(*p)|x
#define atom_xor(p,x) // (*p)=(*p)^x

// integer functions
#define abs(x) // |x|
#define clz(x) // count leading 0 bits, slow, instead use as_uint((float)((x&0x07FF)<<12))>>23
#define mad_sat(a,b,c) // a*b+c
#define max(x,y)
#define min(x,y)

// integer and floating-point functions
#define clamp(x,a,b)
#define sign(x)

// floating-point functions
#define acos(x)
#define acosh(x)
#define acospi(x) // acos(x)/pi
#define asin(x)
#define asinh(x)
#define asinpi(x) // asin(x)/pi
#define atan(x)
#define atan2(x,y) // atan(x/y)
#define atanh(x)
#define atanpi(x) // atan(x)/pi
#define atan2pi(x,y) // atan(x/y)/pi
#define cbrt(x) // x^(1/3)
#define copysign(x,y) // x with sign changed to sign of y
#define cos(x)
#define cosh(x)
#define cospi(x) // cos(pi*x)
#define degrees(x) // x*180/pi
#define erfc(x) // complementary error function
#define erf(x) // error function
#define exp(x) // e^x
#define exp2(x) // 2^x
#define exp10(x) // 10^x
#define expm1(x) // e^x-1
#define fabs(x) // |x|
#define fdim(x,y) // max(x-y,0)
#define floor(x) // (float)((int)x)
#define fma(a,b,c) // a*b+c
#define fmax(x,y) // max(x,y)
#define fmin(x,y) // min(x,y)
#define fmod(x,y) // x%y
#define hypot(x,y) // (x^2+y^2)^(1/2)
#define isfinite(x) // test for finite value
#define isinf(x) // test for infinity
#define isnan(x) // test for NaN
#define isnormal(x) // test for normal value
#define ldexp(x,n) // x*2^n (n is integer)
#define lgamma(x) // log gamma function
#define log(x) // ln(x)
#define log2(x) // log_2(x)
#define log10(x) // log_10(x)
#define log1p(x) // ln(1+x)
#define mad(a,b,c) // a*b+c (approximation)
#define maxmag(x,y) // max(|x|,|y|)
#define minmag(x,y) // min(|x|,|y|)
#define native_rsqrt(x) // x^(-1/2)
#define native_sqrt(x) // x^(1/2)
#define pow(x,y) // x^y
#define pown(x,n) // x^n, where n is an integer
#define powr(x,y) // x^y, where x>=0
#define radians(x) // x*pi/180
#define rootn(x,y) // x^(1/y)
#define rsqrt(x) // x^(-1/2), slower, use native_rsqrt(x) instead
#define signbit(x) // test for sign bit
#define sin(x)
#define sinh(x)
#define sinpi(x) // sin(pi*x)
#define sqrt(x) // x^(1/2), slower, use native_sqrt(x) instead
#define step(x,y) // y<x ? 0 : 1
#define stepsmooth(a,b,x) // step and interpolate
#define tan(x)
#define tanh(x)
#define tanpi(x) // tan(pi*x)
#define tgamma(x) // gamma function
#define vload_half(o,p) // load half from global memory
#define vstore_half_rte(x,o,p) // store half in global memory

// vector functions
#define cross(x,y) // xxy
#define distance(x,y) // |y-x|
#define dot(x,y) // x*y
#define length(x) // |x|
#define normalize(x) // x/|x|
#define fast_distance(x,y) // |y-x|
#define fast_length(x) // |x|
#define fast_normalize(x) // x/|x|

// data types
#define half
#define half2
#define half3
#define half4
#define half8
#define half16
#define float2
#define float3
#define float4
#define float8
#define float16
#define double2
#define double3
#define double4
#define double8
#define double16
#define char2
#define char3
#define char4
#define char8
#define char16
#define short2
#define short3
#define short4
#define short8
#define short16
#define int2
#define int3
#define int4
#define int8
#define int16
#define long2
#define long3
#define long4
#define long8
#define long16
#define uchar
#define uchar2
#define uchar3
#define uchar4
#define uchar8
#define uchar16
#define ushort
#define ushort2
#define ushort3
#define ushort4
#define ushort8
#define ushort16
#define uint
#define uint2
#define uint3
#define uint4
#define uint8
#define uint16
#define ulong
#define ulong2
#define ulong3
#define ulong4
#define ulong8
#define ulong16

// interpret functions
#define as_half(x)
#define as_half2(x)
#define as_half3(x)
#define as_half4(x)
#define as_half8(x)
#define as_half16(x)
#define as_float(x)
#define as_float2(x)
#define as_float3(x)
#define as_float4(x)
#define as_float8(x)
#define as_float16(x)
#define as_double(x)
#define as_double2(x)
#define as_double3(x)
#define as_double4(x)
#define as_double8(x)
#define as_double16(x)
#define as_char(x)
#define as_char2(x)
#define as_char3(x)
#define as_char4(x)
#define as_char8(x)
#define as_char16(x)
#define as_short2(x)
#define as_short3(x)
#define as_short4(x)
#define as_short8(x)
#define as_short16(x)
#define as_int(x)
#define as_int2(x)
#define as_int3(x)
#define as_int4(x)
#define as_int8(x)
#define as_int16(x)
#define as_long(x)
#define as_long2(x)
#define as_long3(x)
#define as_long4(x)
#define as_long8(x)
#define as_long16(x)
#define as_uchar(x)
#define as_uchar2(x)
#define as_uchar3(x)
#define as_uchar4(x)
#define as_uchar8(x)
#define as_uchar16(x)
#define as_ushort(x)
#define as_ushort2(x)
#define as_ushort3(x)
#define as_ushort4(x)
#define as_ushort8(x)
#define as_ushort16(x)
#define as_uint(x)
#define as_uint2(x)
#define as_uint3(x)
#define as_uint4(x)
#define as_uint8(x)
#define as_uint16(x)
#define as_ulong(x)
#define as_ulong2(x)
#define as_ulong3(x)
#define as_ulong4(x)
#define as_ulong8(x)
#define as_ulong16(x)