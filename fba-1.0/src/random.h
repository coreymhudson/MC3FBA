#ifndef RNG_H
#define RNG_H

// __________________________________________________________________________
// random.h   - a Random Number Generator Class
// random.cpp - contains the non-inline class methods

// __________________________________________________________________________
// This C++ code uses the simple, very fast "KISS" (Keep It Simple
// Stupid) random number generator suggested by George Marsaglia in a
// Usenet posting from 1999.  He describes it as "one of my favorite
// generators".  It generates high-quality random numbers that
// apparently pass all commonly used tests for randomness.  In fact, it
// generates random numbers by combining the results of three other good
// random number generators that have different periods and are
// constructed from completely different algorithms.  It does not have
// the ultra-long period of some other generators - a "problem" that can
// be fixed fairly easily - but that seems to be its only potential
// problem.  The period is about 2^123.

// The ziggurat method of Marsaglia is used to generate exponential and
// normal variates.  The method as well as source code can be found in
// the article "The Ziggurat Method for Generating Random Variables" by
// Marsaglia and Tsang, Journal of Statistical Software 5, 2000.

// The method for generating gamma variables appears in "A Simple Method
// for Generating Gamma Variables" by Marsaglia and Tsang, ACM
// Transactions on Mathematical Software, Vol. 26, No 3, Sep 2000, pages
// 363-372.

// The code for Poisson and Binomial random numbers comes from
// Numerical Recipes in C.

// Some of this code is unlikely to work correctly as is on 64 bit
// machines.

#include <cstdlib>
#include <ctime>
#include <string.h>
#include <math.h>
#include <stdio.h>
#ifdef _WIN32
#include <process.h>
#define getpid _getpid
#else
#include <unistd.h>
#endif

//#ifdef _WIN32
  static const double PI   =  3.1415926535897932;
  static const double AD_l =  0.6931471805599453;
  static const double AD_a =  5.7133631526454228;
  static const double AD_b =  3.4142135623730950;
  static const double AD_c = -1.6734053240284925;
  static const double AD_p =  0.9802581434685472;
  static const double AD_A =  5.6005707569738080;
  static const double AD_B =  3.3468106480569850;
  static const double AD_H =  0.0026106723602095;
  static const double AD_D =  0.0857864376269050;
//#endif //_WIN32

namespace KW_RNG {

typedef signed int  sint;
typedef unsigned int uint;
typedef signed long  slong;
typedef unsigned long ulong;

class RNG
{
private:
  ulong z, w, jsr, jcong; // Seeds

  ulong kn[128], ke[256];
  double wn[128],fn[128], we[256],fe[256];

/*
#ifndef _WIN32
  static const double PI   =  3.1415926535897932;
  static const double AD_l =  0.6931471805599453;
  static const double AD_a =  5.7133631526454228;
  static const double AD_b =  3.4142135623730950;
  static const double AD_c = -1.6734053240284925;
  static const double AD_p =  0.9802581434685472;
  static const double AD_A =  5.6005707569738080;
  static const double AD_B =  3.3468106480569850;
  static const double AD_H =  0.0026106723602095;
  static const double AD_D =  0.0857864376269050;
#endif //_WIN32
*/

public:
  RNG() { init(); zigset(); }
  RNG(ulong z_, ulong w_, ulong jsr_, ulong jcong_ ) :
    z(z_), w(w_), jsr(jsr_), jcong(jcong_) { zigset(); }
  ~RNG() { }


  inline ulong znew() 
    { return (z = 36969 * (z & 65535) + (z >> 16)); }
  inline ulong wnew() 
    { return (w = 18000 * (w & 65535) + (w >> 16)); }
  inline ulong MWC()  
    { return (((znew() & 65535) << 16) + wnew()); }
  inline ulong SHR3()
    { jsr ^= ((jsr & 32767) << 17); jsr ^= (jsr >> 13); return (jsr ^= ((jsr << 5) & 0xFFFFFFFF)); }
  inline ulong CONG() 
    { return (jcong = (69069 * jcong + 1234567) & 0xFFFFFFFF); }
  inline double RNOR() {
    slong h = rand_int32();
    ulong i = h & 127;
    return (((ulong) abs((sint) h) < kn[i]) ? h * wn[i] : nfix(h, i));
  }
  inline double REXP() {
    ulong j = rand_int32();
    ulong i = j & 255;
    return ((j < ke[i]) ? j * we[i] : efix(j, i));
  }

  double nfix(slong h, ulong i);
  double efix(ulong j, ulong i);
  void zigset();

  inline void init()
    { ulong yo = time(0) + getpid();
      z = w = jsr = jcong = yo; }
  inline void init(ulong z_, ulong w_, ulong jsr_, ulong jcong_ )
    { z = z_; w = w_; jsr = jsr_; jcong = jcong_; }

  inline ulong rand_int32()         // [0,2^32-1]
    { return ((MWC() ^ CONG()) + SHR3()) & 0xFFFFFFFF; }
  inline long rand_int31()          // [0,2^31-1]
    { return long(rand_int32() >> 1);}
  inline double rand_closed01()     // [0,1]
    { return ((double) rand_int32() / 4294967295.0); }
  inline double rand_open01()       // (0,1)
    { return (((double) rand_int32() + 0.5) / 4294967296.0); }
  inline double rand_halfclosed01() // [0,1)
    { return ((double) rand_int32() / 4294967296.0); }
  inline double rand_halfopen01()   // (0,1]
    { return (((double) rand_int32() + 0.5) / 4294967295.5); }

  // Continuous Distributions
  inline double uniform(double x = 0.0, double y = 1.0)
    { return rand_closed01() * (y - x) + x; }
  inline double normal(double mu = 0.0, double sd = 1.0)
    { return RNOR() * sd + mu; }
  inline double exponential(double lambda = 1)
    { return REXP() / lambda; }
  double gamma(double shape = 1, double scale = 1);
  double chi_square(double df)
    { return gamma(df / 2.0, 0.5); }
  double beta(double a1, double a2)
    { double x1 = gamma(a1, 1); return (x1 / (x1 + gamma(a2, 1))); }

  // Discrete Distributions
  double poisson(double lambda);
  int binomial(double pp, int n);
	
}; // class RNG
	
} // namespace

static KW_RNG::RNG _my_random;
inline int my_random() { return _my_random.rand_int31(); }
inline void my_srandom(int x) { _my_random.init(x,!x*13,x*x+1,(x>>16)+(x<<16)); }
inline int my_binomial(double pp, int n) { return _my_random.binomial(pp,n); }
inline double my_random01() { return _my_random.rand_halfopen01(); }
	
	// Max line size in files
#define FBUFF_SIZE 1000000
	
#define MY_RAND_MAX 0x7FFFFFFF
	
	// Min & Max
#ifndef min
#define defmin(type) inline type min(type a, type b) { return a<b ? a : b; }
	defmin(int)
	defmin(double)
	defmin(unsigned long)
#endif //min
#ifndef max
#define defmax(type) inline type max(type a, type b) { return a>b ? a : b; }
	defmax(int)
	defmax(double)
	defmax(unsigned long)
#endif //max
	
	// Max Int
#ifndef MAX_INT
#define MAX_INT 0x7FFFFFFF
#endif //MAX_INT
	
	// Tag Int
#define TAG_INT 0x40000000
	
	// Oldies ....
#define S_VECTOR_RAW
	
	//Fast search or replace
inline int* fast_rpl(int *m, const int a, const int b) {
	while(*m!=a) m++;
	*m = b;
	return m;
}
inline int* fast_search(int *m, const int size, const int a) {
	int *p = m+size;
	while(m != p--) if(*p == a) return p;
	return NULL;
}
	
	// 1/(RANDMAX+1)
#define inv_RANDMAX (1.0/(1.0+double(MY_RAND_MAX)))
	
	// random number in ]0,1[, _very_ accurate around 0
inline double random_float() {
	int r=my_random();
	double mul=inv_RANDMAX;
	while(r<=0x7FFFFF) {
		r<<=8;
		r+=(my_random()&0xFF);
		mul*=(1.0/256.0);
	}
	return double(r)*mul;
}
	
	// Return true with probability p. Very accurate when p is small.
#define test_proba(p) (random_float()<(p))
	
// Random bit generator, sparwise.
static int _random_bits_stored = 0;
static int _random_bits = 0;

inline int random_bit() {
	register int a = _random_bits;
	_random_bits = a >> 1;
	if(_random_bits_stored--) return a&0x1;
	a = my_random();
	_random_bits = a >> 1;
	_random_bits_stored = 30;
	return a&0x1;
}

/*inline void r_choose_n(const std::vector<int> samples, *int results[], int n){
	int results[n];
	int prov[samples.size()];
	int count=0;
	for(int i = 0; i < samples.size(); i++){
		if(samples.at(i) > 0){
			prov[count] = i;
			count++;
		}
	}
	int j, k, tmp;
	for(j = count-1; j > 0; j--){
		k = my_random() % j;
		tmp = prov[k];
		prov[k] = prov[j];
		prov[j] = tmp;
	}
	for(int i = 0; i < n; i++){
		results[i] = prov[i];
	}
}*/

#endif // RNG_H

