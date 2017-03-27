#ifndef UTILITY_H
#define UTILITY_H

#include <iostream>
#include <ctime>
#undef near

using std::vector;

namespace Utility
{
///////////////////////////////////////////////////////////////////////////////////////////////
// Returns true if the value v is between min and max
///////////////////////////////////////////////////////////////////////////////////////////////
template<typename T>
inline bool InRange(const T &v, const T &min, const T &max)
{
  return v >= min && v <= max;
}
///////////////////////////////////////////////////////////////////////////////////////////////
// Returns true if the value v is between min and max within the tolerance TOL
///////////////////////////////////////////////////////////////////////////////////////////////
template<typename T>
inline bool InRange(const T &v, const T &min, const T &max, const T &TOL)
{
  return InRange(v, min-TOL, max+TOL);
}
///////////////////////////////////////////////////////////////////////////////////////////////
// Returns the value with the larger absolute value
///////////////////////////////////////////////////////////////////////////////////////////////
template<typename T>
inline T MaxAbs(const T &a, const T &b)
{
  return (a*a > b*b) ? a : b;
}
///////////////////////////////////////////////////////////////////////////////////////////////
// Returns the value with the smaller absolute value
///////////////////////////////////////////////////////////////////////////////////////////////
template<typename T>
inline T MinAbs(const T &a, const T &b)
{
  return (a*a < b*b) ? a : b;
}
///////////////////////////////////////////////////////////////////////////////////////////////
// A simple Matlab like timing function
///////////////////////////////////////////////////////////////////////////////////////////////
inline double tictoc()
{
  double t_new = clock() / (double) CLOCKS_PER_SEC;
  static double t = t_new;
  double dt = t_new - t;
  t = t_new;
  return dt; 
}
///////////////////////////////////////////////////////////////////////////////////////////////////
// Binary search between low and high for value (assumes array is sorted)
///////////////////////////////////////////////////////////////////////////////////////////////////
template<typename T>
inline int SearchSorted(const T* A, T value, unsigned low, unsigned high)
{
  unsigned mid;
  while ((low+1)<=(high+1))
  {
   mid = (low + high) >> 1;                 // set mid to mean of high and low
   if (A[mid] > value)       high = mid-1;  // mid was too high, reduce high
   else if (A[mid] < value)  low = mid+1;   // mid was too low, increase low
   else return mid;
  }
  return -1; // value not found in array, return -1
}
///////////////////////////////////////////////////////////////////////////////////////////////////
// Flat search between low and high for value (assumes array is NOT sorted)
///////////////////////////////////////////////////////////////////////////////////////////////////
template<typename T>
inline int SearchUnsorted(const T* A, T value, unsigned low, unsigned high)
{
  for (low; low<=high; low++) if (A[low] == value) return low;
  return -1;
}
///////////////////////////////////////////////////////////////////////////////////////////////////
// Regular swap
///////////////////////////////////////////////////////////////////////////////////////////////////
template<typename T>
inline void Swap(T &a, T &b)
{
  T temp = a;
  a = b;
  b = temp;
}
///////////////////////////////////////////////////////////////////////////////////////////////////
// Returns the sign of a double precision number
///////////////////////////////////////////////////////////////////////////////////////////////////
inline double Sign(double x)  { return (x >= 0.0) ? 1.0 : -1.0; }
///////////////////////////////////////////////////////////////////////////////////////////////////
// Returns the sign of a single precision number
///////////////////////////////////////////////////////////////////////////////////////////////////
inline float Sign(float x)  { return (x >= 0.0f) ? 1.0f : -1.0f; }
///////////////////////////////////////////////////////////////////////////////////////////////////
// Returns the pythagorean distance with the two lengths
///////////////////////////////////////////////////////////////////////////////////////////////////
inline double Norm(double x, double y) { return sqrt(x*x+y*y); }
///////////////////////////////////////////////////////////////////////////////////////////////////
// Returns square of a value
///////////////////////////////////////////////////////////////////////////////////////////////////
template<typename T>
inline T Square(T x) { return x*x; }
}

#endif
