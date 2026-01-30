/*!
  \file auto_axes.hpp
  \brief Scalable Vector Graphic (SVG) autoscaling of axes.
  \details Inspect container or data values to find minimum and maximum,
    avoiding values that are NaN and/or 'at limit'.
    Scale axis using max and min values (calculated or user provided),
    optionally to include the orgin, and to set the ticks.
    Provide fine control over any overlap at the edges of the axes to avoid a tiny
    amount over the limit resulting in an ugly extra major tick.
    Also allow optional forcing of the ticks to be multiples of 1, 2, 5, 10.
  \author
*/

// Copyright Paul A. Bristow 2006 - 2013.
// Copyright Ali Q. Raeini 2019

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_SVPLOT_AXES_HPP
#define BOOST_SVPLOT_AXES_HPP


#include <algorithm>
 // minmax_element finds both min and max elements more efficiently than separately.

#include <cmath> // using std::fabs, std::pow, std::ceil, std::log10
#include <limits> // using std::numeric_limits;
#include <stdexcept> // using std::domain_error; sdtd::runtime_error;
#include <iterator> // using std::iterator_traits;
#include <utility> // using std::pair; using std::make_pair;
#include <typses.h> // using std::pair; using std::make_pair;

							namespace svg _begins_

class XAxis
{ public:
	const double& operator()(const dbl2& pr) const {return pr.a;}
	bool operator()(const dbl2& p1, const dbl2& p2) const {return p1.a<p2.a;}
	static const char name='X';
};
class YAxis
{ public:
	const double& operator()(const dbl2& pr) const {return pr.b;}
	bool operator()(const dbl2& p1, const dbl2& p2) const {return p1.b<p2.b;}
	static const char name='Y';
};


 /*! \brief Inspect values to find min and max.
     \details Inspect all values between begin and (one before) end to work out and update min and max.
       Similar to boost::minmax_element, but ignoring at 'limit': non-finite, +-infinity, max & min, & NaN).
    \tparam Iter InputIterator into STL container.
    \param frs Iterator to chosen first item in container.
    \param lst Iterator to chosen last item in container.
    \param minv Updated with Minimum value found (not 'at limit').
    \param maxv Updated with Maximum value found (not 'at limit').
    \return number of normal values (not 'at limit' neither too big, NaN nor infinite).
  */
template <typename Iter, class Axis>
int mnmx(const Axis& Ax,
  Iter frs, // iterator to chosen first item in container.
  Iter lst,  // iterator to chosen last item in container.
  double* minv, // Updated with Minimum value found (not 'at limit').
  double* maxv) // Updated with Maximum value found (not 'at limit').
{
  *maxv = std::numeric_limits<double>::quiet_NaN();
  *minv = std::numeric_limits<double>::quiet_NaN();
  int nGood = 0; // Count of values within limits.
  Iter pos = frs;
  while(pos != lst && !std::isfinite(Ax(*pos)))  ++pos;   // skip those outside limits
  if (pos != lst)  {
    double x = Ax(*pos);
    *maxv = x;    *minv = x;
    ++nGood;
    while( ++pos != lst )  {
      x = Ax(*pos);
      if (std::isfinite(x))  {
        if      (x > *maxv)  *maxv = x;
        else if (x < *minv)  *minv = x;
        ++nGood;
      }
    }
  }
  dAsrt(nGood==(lst-frs), _s(lst-frs-nGood)+" non finites out of "+_s(lst-frs)+" "+Ax.name+" values");

  return nGood; // for check
} //  mnmx


//! Scale axis from data series (usually to plot), internal
//! \tparam Iter Type of interator into STL container type: @c array, @c vector ...
//! \param frs First item in container to use to calculate autoscale mimimum or maximum.
//! \param lst Last item in container to use to calculate autoscale mimimum or maximum.
//! \param axis_min_value Computed minimum value for the axis, updated by scale_axis.
//! \param axis_max_value Computed maximum value for the axis, updated by scale_axis.
//! \param major_ticks_fraction_  Computed relative size of tick increments for the axis, updated by scale_axis.
//! \param check_limits Whether to check all values for infinity, NaN etc.
//! \param autoscale_plusminus Multiplier of uncertainty or standard deviations to allow for confidence ellipses.
//! \param origin If false, do not include the origin unless the range min_value <= 0 <= max_value.
//! \param tight fraction of 'overrun' allowed before another tick used. For visual effect up to about 0.001 might suit a 1000 pixel wide image, allowing values just 1 pixel over the tick to be shown.
//! \param precision Minimum number of major ticks. TODO clarify the link
//! \param scaleMinMax (?) forgots
template <typename Iter,class Axis>
void scale_axis(const Axis& Ax,
  Iter frs, // Iterator into begin in STL container.
  Iter lst, // Iterator into end in STL container.  // (not necessarily ordered by size, so will find min and max).
  double* axis_min_value, // Computed minimum value for the axis, updated here.
  double* axis_max_value, //  Computed maximum value for the axis, updated here.
  double* major_ticks_fraction_, //  Computed tick increment for the axis, updated here.
  bool check_limits, // Whether to check all values for infinity, NaN etc.
  double autoscale_plusminus, // NOT_USED Mutiplier of uncertainty or standard deviations to allow for confidence ellipses.
  bool origin = false, // Do not include the origin unless the range min_value <= 0 <= max_value.
  double tight = 0., // NOT_USED tightest - fraction of 'overrun' allowed before another tick used.
  int precision = 3, // Minimum number of major ticks.
  bool scaleMinMax = false
)
{
	bool scaleMTicks = *major_ticks_fraction_<=0.0;
	double minv;
	double maxv;
	if (!check_limits)
	{ // minmax_element is efficient for maps because can use knowledge of being sorted,
		// BUT only if it can be assumed that no values are 'at limits',
		// infinity, NaN, max_value, min_value, denorm_min.
		// Otherwise it is necessary to inspect all values individually.
		std::pair<Iter, Iter> result = std::minmax_element(frs, lst, Ax); // min & max
		// scale_axis (not check_limits version) forward declaration to ensure compiler finds right version.
		minv = result.first->b;
		maxv = result.second->b;
	}
	else
	{ // Must check limits.
		int nGood = mnmx(Ax,frs, lst, &minv, &maxv);
		dAsrt(nGood >=2, "Could not find min & max for "+Ax.name+" axis!"); //throw
		minv = std::min(*axis_min_value,minv); //
		maxv = std::max(*axis_max_value,maxv);
		if(nGood<2) {minv-=1.0e-64+1.0e-15*minv;  maxv+=1.0e-64+1.0e-15*maxv;}
	}

	using std::isfinite;
	dAsrt(isfinite(minv), Ax.name+"_min not finite!");
	dAsrt(isfinite(maxv), Ax.name+"_max not finite!");
	if(!isfinite(minv) || !isfinite(maxv)) return;

	if (origin == true)	{ // include zero  per user request
		if (minv > 0.)    minv = 0.;
		else if(maxv < 0.)  maxv = 0.;
	}

	double del = maxv-minv;
	int order = std::max(precision+5,-int(log10(maxv - minv)));
	minv=std::min(minv,(round((maxv-1e-16*del)*pow(10.0,order))-2.0)*pow(10.0,-order)-1.0e-72);
	maxv=std::max(maxv,(round((minv+1e-16*del)*pow(10.0,order))+4.0)*pow(10.0,-order)+2.0e-72);
	dAsrt(minv < maxv, "min > max!",-1);// max and min are transposed!


	double Delta = maxv - minv;
	order = int(floor(log10(Delta))); // 0 to 9.999, gives 0, 10 to 99.9 gives 1 ...
	double scaled_min = minv * pow(10., -order); // 0 to 9.99 is unchanged, 10 to 9.99 scaled down to 1. to 9.99
	double scaled_max = maxv * pow(10., -order); // 0 to 9.99 is unchanged, 10 to 9.99 scaled down to 1. to 9.99

	del = scaled_max - scaled_min;
	double offscal=del/std::max(-scaled_min, scaled_max); // low values require more decimal points in axis labels, ideally should be a function of font size
	if(!scaleMTicks) del *= *major_ticks_fraction_;
	else if(del<1.5) del = 0.25;// del converts from Delta(1-10) to delta(0-1)
	else if(del<1.999 || (del<2.999 && offscal>0.5)) del = 0.5;
	else if(del<4.999 || (del<5.999 && offscal>0.5))   del = 1;
	else del = 2;

	if(scaleMinMax)
	{
		maxv = ceil(scaled_max/del)*del;// * 10^... missing
		if((maxv-scaled_max)/del>1) maxv -= int((maxv-scaled_max)/del)*del; // between 1 and 1.5
		*axis_max_value = maxv *    pow(10., order); // Ensure both axis and ticks have the *same* range.

		*axis_min_value = (floor(scaled_min/del)*del)*    pow(10., order); // (To use the separation, made to give the potential for different ranges,
	}
	if(scaleMTicks)
		*major_ticks_fraction_ = -del*pow(10., order) / (*axis_max_value - *axis_min_value); //absoloute value will be taken later, - indicates autoscale

	dAsrt(*major_ticks_fraction_<0.5, "major_ticks_fraction_ < 2.01");

} //   scal axis



								_end_of_(namespace svg)

#endif // BOOST_SVPLOT_AXES_HPP
