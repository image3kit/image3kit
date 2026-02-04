// Modified by A.Q.Raeini, Jan 2018
/*
 * Copyright (c) 2006-2015, Visualization and Multimedia Lab,
 *						  University of Zurich <http://vmml.ifi.uzh.ch>,
 *						  Eyescale Software GmbH,
 *						  Blue Brain Project, EPFL
 *
 * This file is part of VMMLib <https://github.com/VMML/vmmlib/>
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.  Redistributions in binary
 * form must reproduce the above copyright notice, this list of conditions and
 * the following disclaimer in the documentation and/or other materials provided
 * with the distribution.  Neither the name of the Visualization and Multimedia
 * Lab, University of Zurich nor the names of its contributors may be used to
 * endorse or promote products derived from this software without specific prior
 * written permission.
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef VMMLIB__VECTOR__HPP
#define VMMLIB__VECTOR__HPP

//#include <vmmlib/vmmlib_config.hpp>
//#include <vmmlib/math.hpp>
//#include <vmmlib/std::enable_if.hpp>
//#include <vmmlib/exception.hpp>

#include <iostream>
#include <iomanip>
#include <cstring>
#include <limits>
#include <algorithm>
#include <cstddef>
#include <cassert>
#include <cmath>




		/**
		 *   heavily inspired by boost::std::enable_if
		 *   http://www.boost.org, file: boost/utility/std::enable_if.hpp,
		 *   Copyright 2003 Jaakko Järvi, Jeremiah Willcock, Andrew Lumsdaine
		 */
		template < bool condition, typename T = void >
		struct enable_iff { typedef T type; };

		template< typename T >
		struct enable_iff<false,T> {};



template<typename T> struct var3;


#define _MT_ template< size_t M, typename T  >


template< size_t M, typename T = double >
class Vctr
{
public:
	typedef T									   value_type;
	typedef T*									  iterator;
	typedef const T*								const_iterator;
	typedef std::reverse_iterator< iterator >	   reverse_iterator;
	typedef std::reverse_iterator< const_iterator > const_reverse_iterator;

	static const size_t DIMENSION = M;

	// constructors
	Vctr() : array() {} // http://stackoverflow.com/questions/5602030
	Vctr(const Vctr&) = default;
	explicit Vctr(const T& a); // sets all components to a;
	Vctr(const T& x, const T& y);
	Vctr(const T& x, const T& y, const T& z);
	Vctr(const var3<T>& v3);
	Vctr(const T& x, const T& y, const T& z, const T& w);

#ifndef SWIG
	// initializes the first M-1 values from vector_, the last from last_
	Vctr(const Vctr< M-1,T>& vector_, T last_);
#endif

	Vctr(const T* values);

#ifdef __OSG_MATH
	template< typename OSGVEC3 >
	explicit Vctr(const OSGVEC3& from,	typename std::enable_if< M==3, OSGVEC3 >::type* = 0);
#endif

	// vec< M > with homogeneous coordinates <-> vec< M-1 > conversion ctor
	// to-homogenous-coordinates ctor
	template< size_t N >
	Vctr(const Vctr<N,T>& source_, typename std::enable_if< N==M - 1 >::type* = 0);

	// from-homogenous-coordinates vector
	template< size_t N >
	Vctr(const Vctr<N,T>& source_, typename std::enable_if< N==M + 1 >::type* = 0 );

	template< typename U > Vctr(const Vctr<M, U >& source_);

	// iterators
	inline iterator begin();
	inline iterator end();
	inline const_iterator begin() const;
	inline const_iterator end() const;
	inline reverse_iterator rbegin();
	inline reverse_iterator rend();
	inline const_reverse_iterator rbegin() const;
	inline const_reverse_iterator rend() const;

#  ifndef VMMLIB_NO_CONVERSION_OPERATORS
	// conversion operators
	inline operator T*();
	inline operator const T*() const;
#  else
	inline T& operator[](size_t index);
	inline const T& operator[](size_t index) const;
#  endif

	// accessors
	inline T& operator()(size_t index);
	inline const T& operator()(size_t index) const;

	inline T& at(size_t index);
	inline const T& at(size_t index) const;

	// element accessors for M <= 4;
	inline T& x();
	inline T& y();
	inline T& z();
	inline T& w();
	inline const T& x() const;
	inline const T& y() const;
	inline const T& z() const;
	inline const T& w() const;

	// pixel color element accessors for M<= 4
	inline T& r();
	inline T& g();
	inline T& b();
	inline T& a();
	inline const T& r() const;
	inline const T& g() const;
	inline const T& b() const;
	inline const T& a() const;

	bool operator==(const Vctr& other) const;
	bool operator!=(const Vctr& other) const;
	bool equals(const Vctr& other,  T tolerance = std::numeric_limits<T>::epsilon()) const;
	bool operator<(const Vctr& other) const;

	// remember kids: c_arrays are dangerous and evil!
	Vctr& operator=(const T* c_array);
	T operator=(T filler);

	Vctr& operator=(const Vctr& other);

	// returns void to avoid 'silent' loss of precision when chaining
	template< typename U > void operator=(const Vctr<M, U >& other);

	// to-homogenous-coordinates assignment operator
	// non-chainable because of sfinae
	template< size_t N >
	typename std::enable_if< N==M - 1 >::type*
		operator=(const Vctr< N,T>& source_);

	// from-homogenous-coordinates assignment operator
	// non-chainable because of sfinae
	template< size_t N >
	typename std::enable_if< N==M + 1 >::type*
		operator=(const Vctr< N,T>& source_);

	Vctr operator*(const Vctr& other) const;
	Vctr operator/(const Vctr& other) const;
	Vctr operator+(const Vctr& other) const;
	Vctr operator-(const Vctr& other) const;

	void operator*=(const Vctr& other);
	void operator/=(const Vctr& other);
	void operator+=(const Vctr& other);
	void operator-=(const Vctr& other);

	Vctr operator*(const T other) const;
	Vctr operator/(const T other) const;
	Vctr operator+(const T other) const;
	Vctr operator-(const T other) const;

	void operator*=(const T other);
	void operator/=(const T other);
	void operator+=(const T other);
	void operator-=(const T other);

	Vctr operator-() const;

	const Vctr& negate();

	void set(T a); // sets all components to a;
#ifndef SWIG
	void set(const Vctr< M-1,T>& v, T a);
#endif
	template< size_t N >
	void set(const Vctr< N,T>& v);

	// sets the first few components to a certain value
	void set(T x, T y);
	void set(T x, T y, T z);
	void set(T x, T y, T z, T w);

	template< typename input_iterator_t >
	void iter_set(input_iterator_t begin_, input_iterator_t end_);

	// compute the cross product of two vectors
	// note: there's also a free function:
	// vector<> cross(const vector<>, const vector<>)

	// res = vec1.cross(vec2) => retval res = vec1 x vec2
	template< typename TT >
	Vctr cross(const Vctr<M, TT >& rhs,   typename std::enable_if< M==3, TT >::type* = 0) const;

	// res.cross(vec1, vec2) => (this) = vec1 x vec2
	template< typename TT >
	void cross(const Vctr<M, TT >& a, const Vctr<M, TT >& b, typename std::enable_if< M==3, TT >::type* = 0);


	// compute the dot product of two vectors
	// note: there's also a free function:
	// T dot(const vector<>, const vector<>);
	inline T dot(const Vctr& other) const;


	// normalize the vector
	// note: there's also a free function:
	// vector<> normalize(const vector<>);
	inline T normalize();

	//sets all vector components to random values
	//remember to set srand(seed);
	//if seed is set to -1, srand(seed) was set outside set_random
	//otherwise srand(seed) will be called with the given seed
	void set_random(int seed = -1);

	inline T mag() const;
	inline T magSqr() const;

	inline T distance(const Vctr& other) const;
	inline T squared_distance(const Vctr& other) const;

	/** @return the product of all elements of this vector */
	T product() const;

	template< typename TT >
	Vctr< 3,T> rotate(const T theta, Vctr<M, TT > axis,    typename std::enable_if< M==3, TT >::type* = 0) const;

	// right hand system, CCW triangle
	// (*this) = normal of v0,v1,v2
	void compute_normal(const Vctr& v0, const Vctr& v1, const Vctr& v2);
	// retval = normal of (this), v1, v2
	Vctr compute_normal(const Vctr& v1, const Vctr& v2) const;

	/** @return the sub vector at the given position and length. */
	template< size_t N >
	Vctr< N,T>& get_sub_vector(size_t offset = 0, typename std::enable_if< M >= N >::type* = 0);

	/** @return the sub vector at the given position and length. */
	template< size_t N >
	const Vctr< N,T>& get_sub_vector(size_t offset = 0, typename std::enable_if< M >= N >::type* = 0) const;

	// sphere functions - sphere layout: center xyz, radius w
	template< typename TT >
	inline Vctr< 3,T> project_point_onto_sphere(	const Vctr< 3, TT >& point, typename std::enable_if< M==4, TT >::type* = 0) const;

	// returns a negative distance if the point lies in the sphere
	template< typename TT >
	inline T distance_to_sphere(const Vctr< 3, TT >& point, typename std::enable_if< M==4, TT >::type* = 0) const;

	// plane functions - plane layout; normal xyz, distance w
	template< typename TT >
	inline T distance_to_plane(const Vctr< 3, TT >& point, typename std::enable_if< M==4, TT >::type* = 0) const;

	template< typename TT >
	inline Vctr< 3,T> project_point_onto_plane(	const Vctr< 3, TT >& point, typename std::enable_if< M==4, TT >::type* = 0) const;

	// returns the index of the minimal resp. maximal value in the vector
	size_t	  find_min_index() const;
	size_t	  find_max_index() const;

	// returns the index of the minimal resp. maximal value in the vector
	size_t	  find_abs_min_index() const;
	size_t	  find_abs_max_index() const;

	// returns minimal resp. maximal value in the vector
	T&		  find_min();
	T&		  find_max();
	const T&	find_min() const;
	const T&	find_max() const;

	void clamp(const T& min = 0., const T& max = 1.);

	template< typename TT >
	void scale_to(Vctr<M, TT >& scaled_vector, T min_value = -1., T max_value = 1.) const;

	inline static size_t size(); // returns M

	bool is_unit_vector() const;

	// perturbs each component by randomly + or - the perturbation parameter
	void perturb(T perturbation = 0.0001);

	void sqrt_elementwise();
	double norm() const; //L2 norm

	// computes the reciprocal value for each component, x = 1/x;
	// WARNING: might res in nans if division by 0!
	void reciprocal();
	// computes the reciprocal value for each component, x = 1/x;
	// checks every component for 0, sets to max value if zero.
	void reciprocal_safe();

	template< typename TT >
	void cast_from(const Vctr<M, TT >& other);

	size_t nnz() const;

	// test each component of the vector for isnan and isinf
	//  inline bool is_valid() const; -> moved to class validator

	void read_from_stream(std::istream& is);
	void write_to_stream(std::ostream& os) const;

	friend std::ostream& operator<< (std::ostream& os, const Vctr& vector_)
	{
		const std::ios::fmtflags flags = os.flags();
		const int				prec  = os.precision();

		os.setf(std::ios::right, std::ios::adjustfield);
		os.precision(5);
		os << "[ ";
		for(size_t index = 0; index < M; ++index)
			os << std::setw(10) << vector_.at(index) << " ";
		os << "]";
		os.precision(prec);
		os.setf(flags);
		return os;
	}

	T array[ M ];	//!< storage

#ifndef SWIG
	// Vctr3 defaults
	static const Vctr FORWARD;
	static const Vctr BACKWARD;
	static const Vctr UP;
	static const Vctr DOWN;
	static const Vctr LEFT;
	static const Vctr RIGHT;

	static const Vctr ONE;
	static const Vctr ZERO;

	// Unit vectors
	static const Vctr UNIT_X;
	static const Vctr UNIT_Y;
	static const Vctr UNIT_Z;
#endif

}; // class vector

//
// typedefs and statics
//
#ifndef SWIG
_MT_ const Vctr<M,T> Vctr<M,T>::FORWARD(0, 0, -1);
_MT_ const Vctr<M,T> Vctr<M,T>::BACKWARD(0, 0, 1);
_MT_ const Vctr<M,T> Vctr<M,T>::UP(0, 1, 0);
_MT_ const Vctr<M,T> Vctr<M,T>::DOWN(0, -1, 0);
_MT_ const Vctr<M,T> Vctr<M,T>::LEFT(-1, 0, 0);
_MT_ const Vctr<M,T> Vctr<M,T>::RIGHT(1, 0, 0);
_MT_
const Vctr<M,T> Vctr<M,T>::ONE(static_cast<T>(1));
_MT_ const Vctr<M,T> Vctr<M,T>::ZERO(static_cast<T>(0));
_MT_
const Vctr<M,T> Vctr<M,T>::UNIT_X(1, 0, 0);
_MT_ const Vctr<M,T> Vctr<M,T>::UNIT_Y(0, 1, 0);
_MT_ const Vctr<M,T> Vctr<M,T>::UNIT_Z(0, 0, 1);
#endif

#ifndef VMMLIB_NO_TYPEDEFS

typedef Vctr<2,int> Vctr2i;
typedef Vctr<3,int> Vctr3i;
typedef Vctr<4,int> Vctr4i;
typedef Vctr<2,unsigned> Vctr2ui;
typedef Vctr<3,unsigned> Vctr3ui;
typedef Vctr<4,unsigned> Vctr4ui;
typedef Vctr<3,double> Vctr3d;
typedef Vctr<4,double> Vctr4d;
typedef Vctr<2,float> Vctr2f;
typedef Vctr<3,float> Vctr3f;
typedef Vctr<4,float> Vctr4f;
typedef Vctr<3,uint8_t> Vctr3ub;
typedef Vctr<4,uint8_t> Vctr4ub;

#ifdef VMMLIB_OLD_TYPEDEFS
typedef Vctr<2,float > vec2f;
typedef Vctr<2,double > vec2d;
typedef Vctr<3,float > dbl3f;
typedef Vctr<3,double > dbl3d;
typedef Vctr<4,float > vec4f;
typedef Vctr<4,double > vec4d;
_MT_ using vector = Vctr<M, T>;
#endif

#endif

//
//  some free functions for convenience
//

_MT_ bool equals(const Vctr<M,T>& a, const Vctr<M,T>& b) {	return a.equals(b); }


// allows float * vector, not only vector * float
_MT_ static Vctr<M,T> operator* (T factor, const Vctr<M,T>& vector_) {	return vector_ * factor; }


_MT_ inline T dot(const Vctr<M,T>& first, const Vctr<M,T>& second) {	return first.dot(second); }


_MT_ inline Vctr<M,T> cross(const Vctr< 3,T>& a, const Vctr< 3,T>& b) {	return a.cross(b); }


_MT_ inline Vctr<M,T> normalize(const Vctr<M,T>& vector_)
{
	Vctr<M,T> v(vector_);
	v.normalize();
	return v;
}

template< typename T >
inline Vctr< 4,T> compute_plane(const Vctr< 3,T>& a,  const Vctr< 3,T>& b,  const Vctr< 3,T>& c)
{
	const Vctr< 3,T> cb = b - c;
	const Vctr< 3,T> ba = a - b;

	Vctr< 4,T> plane = cb.cross(ba);
	plane.normalize();
	plane.w() = -plane.x() * a.x() - plane.y() * a.y() - plane.z() * a.z();
	return plane;
}

_MT_ Vctr<M,T>::Vctr(const T& _a)
{
	for(iterator it=begin(), it_end=end(); it!=it_end; ++it)
	{
		*it = _a;
	}
}

_MT_ Vctr<M,T>::Vctr(const T& _x, const T& _y)
{
	array[ 0 ] = _x;
	array[ 1 ] = _y;
}


_MT_ Vctr<M,T>::Vctr(const T& _x, const T& _y, const T& _z)
{
	array[ 0 ] = _x;
	array[ 1 ] = _y;
	array[ 2 ] = _z;
}



_MT_ Vctr<M,T>::Vctr(const T& _x, const T& _y, const T& _z, const T& _w)
{
	array[ 0 ] = _x;
	array[ 1 ] = _y;
	array[ 2 ] = _z;
	array[ 3 ] = _w;
}


_MT_ Vctr<M,T>::Vctr(const T* values)
{
	memcpy(array, values, M * sizeof(T));
}

#ifdef __OSG_MATH
_MT_ template< typename OSGVEC3 >
Vctr<M,T>::Vctr(const OSGVEC3& from, typename std::enable_if< M==3, OSGVEC3 >::type*)
{
	array[ 0 ] = from.x();
	array[ 1 ] = from.y();
	array[ 2 ] = from.z();
}
#endif

#ifndef SWIG // initializes the first M-1 values from vector_, the last from last_
_MT_ Vctr<M,T>::Vctr(const Vctr< M-1,T>& vector_, T last_)
{
	typename Vctr< M-1,T>::const_iterator	 it = vector_.begin(), it_end = vector_.end();

	iterator my_it = begin();

	for(; it!=it_end; ++it, ++my_it)		(*my_it) = *it;

	(*my_it) = last_;
}
#endif



// to-homogenous-coordinates ctor
_MT_  template< size_t N >
Vctr<M,T>::Vctr(const Vctr< N,T>& source_, typename std::enable_if< N==M - 1 >::type*)
{
	(*this) = source_;
}




// from-homogenous-coordinates ctor
_MT_   template< size_t N >
Vctr<M,T>::Vctr(const Vctr< N,T>& source_, typename std::enable_if< N==M + 1 >::type*)
{
	(*this) = source_;
}


_MT_  template< typename U >
Vctr<M,T>::Vctr(const Vctr<M, U >& source_)
{
	(*this) = source_;
}



_MT_ void Vctr<M,T>::set(T _a)
{
	for(iterator it=begin(), it_end=end(); it!=it_end; ++it)		*it = _a;
}


#ifndef SWIG
_MT_ void Vctr<M,T>::set(const Vctr< M-1,T>& v, T _a)
{
	memcpy(array, v.array, sizeof(T) * (M-1));
	at(M-1) = _a;
}
#endif

_MT_ template< size_t N >
void Vctr<M,T>::set(const Vctr< N,T>& v)
{
	size_t minimum = M;
	if (N < M) minimum = N;
	memcpy(array, v.array, sizeof(T) * minimum);
}

_MT_ void Vctr<M,T>::set(T _x, T _y)
{
	array[ 0 ] = _x;
	array[ 1 ] = _y;
}


_MT_ void Vctr<M,T>::set(T _x, T _y, T _z)
{
	array[ 0 ] = _x;
	array[ 1 ] = _y;
	array[ 2 ] = _z;
}



_MT_ void Vctr<M,T>::set(T _x, T _y, T _z, T _w)
{
	array[ 0 ] = _x;
	array[ 1 ] = _y;
	array[ 2 ] = _z;
	array[ 3 ] = _w;
}


_MT_ inline T&  Vctr<M,T>::operator()(size_t index) {	return at(index); }



_MT_ inline const T&  Vctr<M,T>::operator()(size_t index) const {	return at(index); }



_MT_ inline T&  Vctr<M,T>::at(size_t index)
{
	#ifdef VMMLIB_SAFE_ACCESSORS
	if (index >= M)
	{
		VMMLIB_ERROR("at() - index out of bounds", VMMLIB_HERE);
	}
	#endif
	return array[ index ];
}



_MT_ inline const T&  Vctr<M,T>::at(size_t index) const
{
	#ifdef VMMLIB_SAFE_ACCESSORS
	if (index >= M)
	{
		VMMLIB_ERROR("at() - index out of bounds", VMMLIB_HERE);
	}
	#endif
	return array[ index ];
}


#ifndef VMMLIB_NO_CONVERSION_OPERATORS

_MT_ Vctr<M,T>::operator T*() {	return array; }



_MT_ Vctr<M,T>::operator const T*() const {	return array; }
#else

_MT_ T&  Vctr<M,T>::operator[](size_t index) {	return at(index); }

_MT_ const T&  Vctr<M,T>::operator[](size_t index) const {	return at(index); }


#endif


#if 0
_MT_ inline T&  Vctr<M,T>::operator[](size_t index) {	return at(index); }
_MT_ inline const T&  Vctr<M,T>::operator[](size_t index) const {	return at(index); }
#endif


_MT_ Vctr<M,T> Vctr<M,T>::operator*(const Vctr<M,T>& other) const
{
	Vctr<M,T> res;
	for(size_t index = 0; index < M; ++index)		res.at(index) = at(index) * other.at(index);
	return res;
}



_MT_ Vctr<M,T> Vctr<M,T>::operator/(const Vctr<M,T>& other) const
{
	Vctr<M,T> res;
	for(size_t index = 0; index < M; ++index)		res.at(index) = at(index) / other.at(index);
	return res;
}



_MT_ Vctr<M,T> Vctr<M,T>::operator+(const Vctr<M,T>& other) const
{
	Vctr<M,T> res;
	for(size_t index = 0; index < M; ++index)		res.at(index) = at(index) + other.at(index);
	return res;
}



_MT_ Vctr<M,T> Vctr<M,T>::operator-(const Vctr<M,T>& other) const
{
	Vctr<M,T> res;
	for(size_t index = 0; index < M; ++index)		res.at(index) = at(index) - other.at(index);
	return res;
}




_MT_ void  Vctr<M,T>::operator*=(const Vctr<M,T>& other)
{
	for(size_t index = 0; index < M; ++index)		at(index) *= other.at(index);
}



_MT_ void  Vctr<M,T>::operator/=(const Vctr<M,T>& other)
{
	for(size_t index = 0; index < M; ++index)		at(index) /= other.at(index);
}



_MT_ void  Vctr<M,T>::operator+=(const Vctr<M,T>& other)
{
	for(size_t index = 0; index < M; ++index)		at(index) += other.at(index);
}



_MT_ void  Vctr<M,T>::operator-=(const Vctr<M,T>& other)
{
	for(size_t index = 0; index < M; ++index)		at(index) -= other.at(index);
}



_MT_ Vctr<M,T> Vctr<M,T>::operator*(const T other) const
{
	Vctr<M,T> res;
	for(size_t index = 0; index < M; ++index)		res.at(index) = at(index) * other;
	return res;
}



_MT_ Vctr<M,T> Vctr<M,T>::operator/(const T other) const
{
	Vctr<M,T> res;
	for(size_t index = 0; index < M; ++index)		res.at(index) = at(index) / other;
	return res;
}



_MT_ Vctr<M,T> Vctr<M,T>::operator+(const T other) const
{
	Vctr<M,T> res;
	for(size_t index = 0; index < M; ++index)		res.at(index) = at(index) + other;
	return res;
}



_MT_ Vctr<M,T> Vctr<M,T>::operator-(const T other) const
{
	Vctr<M,T> res;
	for(size_t index = 0; index < M; ++index)		res.at(index) = at(index) - other;
	return res;
}




_MT_ void  Vctr<M,T>::operator*=(const T other)
{
	for(size_t index = 0; index < M; ++index)		at(index) *= other;
}



_MT_ void  Vctr<M,T>::operator/=(const T other)
{
	for(size_t index = 0; index < M; ++index)		at(index) /= other;
}



_MT_ void  Vctr<M,T>::operator+=(const T other)
{
	for(size_t index = 0; index < M; ++index)		at(index) += other;
}



_MT_ void  Vctr<M,T>::operator-=(const T other)
{
	for(size_t index = 0; index < M; ++index)		at(index) -= other;
}



_MT_ Vctr<M,T> Vctr<M,T>::operator-() const
{
	Vctr<M,T> v(*this);
	return v.negate();
}



_MT_ const Vctr<M,T>&
Vctr<M,T>::negate()
{
	for(size_t index = 0; index < M; ++index)		array[ index ] = -array[ index ];
	return *this;
}



_MT_ inline T&  Vctr<M,T>::x() {	return array[ 0 ]; }



_MT_ inline T&  Vctr<M,T>::y() {	return array[ 1 ]; }



_MT_ inline T&  Vctr<M,T>::z() {	return array[ 2 ]; }



_MT_ inline T&  Vctr<M,T>::w() {	return array[ 3 ]; }



_MT_ inline const T&  Vctr<M,T>::x() const {	return array[ 0 ]; }



_MT_ inline const T&  Vctr<M,T>::y() const {	return array[ 1 ]; }



_MT_ inline const T&  Vctr<M,T>::z() const {	return array[ 2 ]; }



_MT_ inline const T&  Vctr<M,T>::w() const {	return array[ 3 ]; }


_MT_ inline T&  Vctr<M,T>::r() {	return array[ 0 ]; }



_MT_ inline T&  Vctr<M,T>::g() {	return array[ 1 ]; }



_MT_ inline T&  Vctr<M,T>::b() {	return array[ 2 ]; }



_MT_ inline T&  Vctr<M,T>::a() {	return array[ 3 ]; }



_MT_ inline const T&  Vctr<M,T>::r() const {	return array[ 0 ]; }



_MT_ inline const T&  Vctr<M,T>::g() const {	return array[ 1 ]; }



_MT_ inline const T&  Vctr<M,T>::b() const {	return array[ 2 ]; }



_MT_ inline const T&  Vctr<M,T>::a() const {	return array[ 3 ]; }

// res = vec1.cross(vec2) => res = vec1 x vec2
_MT_ template< typename TT >
inline Vctr<M,T> Vctr<M,T>::cross(const Vctr<M, TT >& rhs,  typename std::enable_if< M==3, TT >::type*) const
{
	Vctr<M,T> res;   	res.cross(*this, rhs);   	return res;
}



// res.cross(vec1, vec2) => (this) = vec1 x vec2
_MT_ template< typename TT >
void Vctr<M,T>::cross(const Vctr<M, TT >& aa, const Vctr<M, TT >& bb, typename std::enable_if< M==3, TT >::type*)
{
	array[ 0 ] = aa.y() * bb.z() - aa.z() * bb.y();
	array[ 1 ] = aa.z() * bb.x() - aa.x() * bb.z();
	array[ 2 ] = aa.x() * bb.y() - aa.y() * bb.x();
}



_MT_ inline T Vctr<M,T>::dot(const Vctr<M,T>& other) const
{
	T tmp = 0.;
	for(size_t index = 0; index < M; ++index)		tmp += at(index) * other.at(index);
	return tmp;
}


_MT_ inline T Vctr<M,T>::normalize()
{
	const T len = mag();

	if (len==0)		return 0;

	const T tmp = 1. / len;
	(*this) *= tmp;
	return len;
}

_MT_ inline T Vctr<M,T>::mag() const {	return std::sqrt(magSqr()); }

_MT_ inline T Vctr<M,T>::magSqr() const
{
	T _magSqr = 0.;
	for(const_iterator it=begin(), it_end=end(); it!=it_end; ++it)
		_magSqr += (*it) * (*it);

	return _magSqr;
}



_MT_ inline T
Vctr<M,T>::distance(const Vctr<M,T>& other) const {	return std::sqrt(squared_distance(other)); }



_MT_ inline T Vctr<M,T>::squared_distance(const Vctr<M,T>& other) const
{
	Vctr<M,T> tmp(*this);
	tmp -= other;
	return tmp.magSqr();
}

_MT_ inline T Vctr<M,T>::product() const
{
	T res = at(0);
	for(size_t i=1; i<M; ++i)		res *= at(i);
	return res;
}

_MT_ void Vctr<M,T>::compute_normal(const Vctr<M,T>& aa,  const Vctr<M,T>& bb,  const Vctr<M,T>& cc)
{
	Vctr<M,T> u,v;
	// right hand system, CCW triangle
	u = bb - aa;
	v = cc - aa;
	cross(u, v);
	normalize();
}



_MT_ Vctr<M,T> Vctr<M,T>::compute_normal(const Vctr<M,T>& bb, const Vctr<M,T>& cc) const
{
	Vctr<M,T> tmp;
	tmp.compute_normal(*this, bb, cc);
	return tmp;
}

_MT_ template< typename TT >
Vctr< 3,T> Vctr<M,T>::rotate(const T theta, Vctr<M, TT > axis, typename std::enable_if< M==3, TT >::type*) const
{
	axis.normalize();

	const T costheta = std::cos(theta);
	const T sintheta = std::sin(theta);

	return Vctr< 3,T>(
		(costheta + (1.0f - costheta) * axis.x() * axis.x()) * x()	+
		((1 - costheta) * axis.x() * axis.y() - axis.z() * sintheta) * y() +
		((1 - costheta) * axis.x() * axis.z() + axis.y() * sintheta) * z(),

		((1 - costheta) * axis.x() * axis.y() + axis.z() * sintheta) * x() +
		(costheta + (1 - costheta) * axis.y() * axis.y()) * y() +
		((1 - costheta) * axis.y() * axis.z() - axis.x() * sintheta) * z(),

		((1 - costheta) * axis.x() * axis.z() - axis.y() * sintheta) * x() +
		((1 - costheta) * axis.y() * axis.z() + axis.x() * sintheta) * y() +
		(costheta + (1 - costheta) * axis.z() * axis.z()) * z());
}


// sphere layout: center xyz, radius w
_MT_ template< typename TT >
inline Vctr< 3,T>
Vctr<M,T>::
project_point_onto_sphere(const Vctr< 3, TT >& point, typename std::enable_if< M==4, TT >::type*) const
{
	const Vctr< 3,T>& _center = get_sub_vector< 3 >(0);

	Vctr< 3,T> projected_point(point);
	projected_point -= _center;
	projected_point.normalize();
	projected_point *= w();
	return _center + projected_point;
}



// sphere layout: center xyz, radius w
_MT_ template< typename TT >
inline T  Vctr<M,T>::distance_to_sphere(const Vctr< 3, TT >& point, typename std::enable_if< M==4, TT >::type*) const
{
	const Vctr< 3,T>& center_ = get_sub_vector< 3 >(0);
	return (point - center_).mag() - w();
}

_MT_ template< size_t N > inline
Vctr< N,T>& Vctr<M,T>::get_sub_vector(size_t offset,    typename std::enable_if< M >= N >::type*)
{
	assert(offset <= M - N);
	return reinterpret_cast< Vctr< N,T>& >(*(begin() + offset));
}

_MT_ template< size_t N > inline
const Vctr< N,T>& Vctr<M,T>::get_sub_vector(size_t offset,  typename std::enable_if< M >= N >::type*) const
{
	assert(offset <= M - N);
	return reinterpret_cast< const Vctr< N,T>& >(*(begin() + offset));
}


// plane: normal xyz, distance w
_MT_ template< typename TT >
inline T Vctr<M,T>::distance_to_plane(const Vctr< 3, TT >& point, typename std::enable_if< M==4, TT >::type*) const
{
	const Vctr< 3,T>& normal = get_sub_vector< 3 >(0);
	return normal.dot(point) + w();
}



// plane: normal xyz, distance w
_MT_ template< typename TT >
Vctr< 3,T> Vctr<M,T>::project_point_onto_plane(const Vctr< 3, TT >& point, typename std::enable_if< M==4, TT >::type*) const
{
	const Vctr< 3,T>& normal = get_sub_vector< 3 >(0);
	return point - (normal * distance_to_plane(point));
}



_MT_ bool Vctr<M,T>::operator==(const Vctr<M,T>& other) const {	return memcmp(array, other.array, sizeof(array))==0; }


_MT_ bool Vctr<M,T>::operator!=(const Vctr<M,T>& other) const {	return ! this->operator==(other); }


_MT_ bool  Vctr<M,T>::equals(const Vctr<M,T>& other, T tolerance) const
{
	for(size_t index = 0; index < M; ++index)
		if(fabs(at(index) - other(index)) >= tolerance)
			return false;
	return true;

}


_MT_ bool  Vctr<M,T>::operator<(const Vctr<M,T>& other) const
{
	for(size_t index = 0; index < M; ++index)
	{
		if (at(index) < other.at(index)) return true;
		if (other.at(index) < at(index)) return false;
	}
	return false;
}


// to-homogenous-coordinates assignment operator non-chainable because of sfinae
_MT_ template< size_t N >  typename std::enable_if< N==M - 1 >::type*
Vctr<M,T>::operator=(const Vctr< N,T>& source_)
{
	std::copy(source_.begin(), source_.end(), begin());
	at(M - 1) = static_cast<T>(1.);
	return 0;
}


// from-homogenous-coordinates assignment operator non-chainable because of sfinae
_MT_ template< size_t N >  typename std::enable_if< N==M + 1 >::type*
Vctr<M,T>::operator=(const Vctr< N,T>& source_)
{
	const T w_reci = static_cast<T>(1.) / source_(M);
	iterator it=begin(), it_end=end();
	for(size_t index = 0; it!=it_end; ++it, ++index)		*it = source_(index) * w_reci;
	return 0;
}


_MT_ Vctr<M,T>& Vctr<M,T>::operator=(const T* c_array)
{
	iter_set(c_array, c_array + M);
	return *this;
}



_MT_ T Vctr<M,T>::operator=(T filler_value)
{
	for(size_t index = 0; index < M; ++index)
		at(index) = filler_value;
	return filler_value;
}




_MT_ Vctr<M,T>& Vctr<M,T>::operator=(const Vctr<M,T>& other)
{
	if(this!=&other)
		memcpy(array, other.array, M * sizeof(T));
	return *this;
}



// returns void to avoid 'silent' loss of precision when chaining
_MT_ template< typename U >
void Vctr<M,T>::operator=(const Vctr<M, U >& source_)
{
	typedef typename Vctr<M, U >::const_iterator u_c_iter;
	u_c_iter it = source_.begin(), it_end = source_.end();
	for(iterator my_it = begin(); it!=it_end; ++it, ++my_it)
		*my_it = static_cast<T>(*it);
}



_MT_ template< typename input_iterator_t >
void  Vctr<M,T>::iter_set(input_iterator_t begin_, input_iterator_t end_)
{
	input_iterator_t in_it = begin_;
	iterator it=begin(), it_end=end();
	for(; it!=it_end && in_it!=end_; ++it, ++in_it)
		(*it) = static_cast<T>(*in_it);
}

_MT_ void Vctr<M,T>::clamp(const T& min, const T& max)
{
	for(size_t i=0; i<M; ++i)
	{
		if(array[i] < min)			array[i] = min;
		if(array[i] > max) 			array[i] = max;
	}
}



_MT_ template< typename TT >
void  Vctr<M,T>::scale_to(Vctr<M, TT >& result_, T min_value, T max_value) const
{
	T range	   = max_value-min_value;
	T half_range  = range * 0.5;
	T scale	   = (1. / range) * static_cast<T>(std::numeric_limits< TT >::max());

	for(size_t index = 0; index < M; ++index)
	{
		result_.at(index)
			= static_cast< TT >((at(index) + half_range) * scale);
	}

}



_MT_ inline size_t   Vctr<M,T>::size() {	return M; }



_MT_ size_t   Vctr<M,T>::find_min_index() const {	return std::min_element(begin(), end()) - begin(); }



_MT_ size_t   Vctr<M,T>::find_max_index() const {	return std::max_element(begin(), end()) - begin(); }



//_MT_ //size_t
//Vctr<M,T>::find_abs_min_index() const
//{
	//return std::min_element(begin(), end(), vmml::math::abs_less<T>()) - begin();
//}



//_MT_ //size_t
//Vctr<M,T>::find_abs_max_index() const
//{
	//return std::max_element(begin(), end(), vmml::math::abs_greater<T>()) - begin();
//}



_MT_ T&        Vctr<M,T>::find_min()   {	return *std::min_element(begin(), end()); }

_MT_ const T&  Vctr<M,T>::find_min() const   {	return *std::min_element(begin(), end()); }

_MT_ T&        Vctr<M,T>::find_max()   {	return *std::max_element(begin(), end()); }

_MT_ const T&  Vctr<M,T>::find_max() const   {	return *std::max_element(begin(), end()); }


_MT_ inline typename Vctr<M,T>::iterator          Vctr<M,T>::begin() {	return array; }

_MT_ inline typename Vctr<M,T>::iterator          Vctr<M,T>::end() {	return array + M; ; }

_MT_ inline typename Vctr<M,T>::const_iterator    Vctr<M,T>::begin() const {	return array; }

_MT_ inline typename Vctr<M,T>::const_iterator    Vctr<M,T>::end() const {	return array + M; }

_MT_ inline typename Vctr<M,T>::reverse_iterator  Vctr<M,T>::rbegin() {	return array + M - 1; }

_MT_ inline typename Vctr<M,T>::reverse_iterator  Vctr<M,T>::rend() {	return array - 1; }

_MT_ inline typename Vctr<M,T>::const_reverse_iterator  Vctr<M,T>::rbegin() const {	return array + M - 1; }

_MT_ inline typename Vctr<M,T>::const_reverse_iterator  Vctr<M,T>::rend() const {	return array - 1; }



_MT_ bool  Vctr<M,T>::is_unit_vector() const
{
	const_iterator it=begin(), it_end=end();
	bool one = false;
	for(; it!=it_end; ++it)
	{
		if (*it==1.)
		{
			if (one)		return false;
			one = true;
		}
		else if (*it!=0.)
		{
			return false;
		}
	}
	return one;
}




_MT_ void  Vctr<M,T>::perturb(T perturbation)
{
	for(iterator it=begin(), it_end=end(); it!=it_end; ++it)
	{
		(*it) += (rand() & 1u) ? perturbation : -perturbation;
	}

}

_MT_ void  Vctr<M,T>::sqrt_elementwise()
{
	for(iterator it=begin(), it_end=end(); it!=it_end; ++it)
	{
		(*it) = std::sqrt(*it);
	}
}



_MT_ void  Vctr<M,T>::reciprocal()
{
	for(iterator it=begin(), it_end=end(); it!=it_end; ++it)
	{
		(*it) = static_cast<T>(1.) / (*it);
	}
}



_MT_ void  Vctr<M,T>::reciprocal_safe()
{
	for(iterator it=begin(), it_end=end(); it!=it_end; ++it)
	{
		T& v = *it;

		if (v==static_cast<T>(0))			v = std::numeric_limits<T>::max();
		else                     			v = static_cast<T>(1.) / v;
	}
}



_MT_  template< typename TT >
void  Vctr<M,T>::cast_from(const Vctr<M, TT >& other)
{
	typedef Vctr<M, TT > vector_tt_type ;
	typedef typename vector_tt_type::const_iterator tt_const_iterator;
	iterator it=begin(), it_end=end();
	tt_const_iterator other_it = other.begin();
	for(; it!=it_end; ++it, ++other_it)		*it = static_cast<T>(*other_it);

}

_MT_ size_t  Vctr<M,T>::nnz() const
{
	size_t counter = 0;

	const_iterator  it=begin(), it_end=end();
	for(; it!=it_end; ++it)		if (*it!=0)  ++counter;

	return counter;
}

_MT_ void  Vctr<M,T>::read_from_stream(std::istream& is)
{
	for(size_t i=0; i<M; ++i)		is >> at(i);
}

_MT_ void  Vctr<M,T>::write_to_stream(std::ostream& os) const
{
	for(size_t i=0; i<M; ++i)		 os << at(i) << " ";
}

_MT_ double  Vctr<M,T>::norm() const
{
	double norm_v = 0.;
	const_iterator it=begin(), it_end=end();
	for(; it!=it_end; ++it)		norm_v += *it * *it;
	return std::sqrt(norm_v);
}

_MT_ void  Vctr<M,T>::set_random(int seed)
{
	if (seed >= 0) 	srand(seed);
	for(size_t i=0; i<M; ++i)	at(i) = -1. + 2. *  double(rand())/double(RAND_MAX) ;
}


_MT_ inline std::istream& operator>> (std::istream& in, Vctr<M,T>& vec) {  for(T &v : vec)  in >> v;  return in;  }

_MT_ inline std::ostream& operator<< (std::ostream& ou, Vctr<M,T>& vec) {  for(T &v : vec)  ou << v<<' ';  return ou;  }




// compatibility with typses
#ifdef TYPSES_H
 template< size_t M, typename T >
 Vctr<M,T>::Vctr(const var3<T>& v3) 	{ array[ 0 ] = v3.x;	array[ 1 ] = v3.y;	array[ 2 ] = v3.z; 	}
#endif

#endif
