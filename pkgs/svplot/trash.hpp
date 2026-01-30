
// from svplot.hpp
#ifdef _MSC_VER
#pragma warning(push)
#  pragma warning (disable : 4800) // Forcing value to bool 'true' or 'false' (performance warning).
#  pragma warning (disable : 4512) // Assignment operator could not be generated.
#endif



/* chaining seters
     //inline svplot&   y_autoscale(bool b)
     //{ //! Set @c true if to autoscale minimum and maximum for Y-axis.
     //  y_autoscale_ = b;
     //  return *this; .
     //}



//
//		 Derived& x_major_tick(double d)
//		 { //! Set interval (Cartesian units) between major ticks.
//			x_ticks_.major_ticks_fraction_ = d;
//		 }


		 //
		 //Derived& x_minor_interval(double interval)
		 //{ //! Set interval between X-axis minor ticks.			// aka x_minor_tick
			//x_ticks_.minor_interval_ = interval;
			//return *this;
		 //}


//
//		 Derived& x_min(double min_x)
//		 { //! Set the minimum value on the X-axis.
//			// Not useful to check here that x_max_ > x_min_ because may not have set x_min_ yet.
//			x_ax_.min_ = min_x;
//			return *this;
//		 }
//
//
//		 Derived& x_max(double x)
//		 { //! Set the maximum value on the X-axis.
//			// Not useful to check here that x_max_ > x_min_ because may not have set x_min_ yet.
//			x_ax_.max_ = x;
//			return *this;
//		 }

//      inline svplot&   y_axis_color(const svg_color& col)
//      { //! Set Y-axis linecolor. (set only stroke color).
//        Self_.g(g_Y_AXIS).style().stroke_color(col);
//        return *this; .
//      }


      //inline svplot&   y_num_major_ticks(double inter)
      //{ //! Set major interval between ticks on Y-axis.
        //y_ticks_.major_ticks_fraction_ = inter;
        //return *this; //! \return reference to svplot to make chainable.
      //}


		// Derived& x_min_ticks(int min_ticks)
		// { //! Set X-axis autoscale to include at least minimum number of ticks (default = 6).
		//	//! Must preceed x_autoscale(data) call.
		//	x_min_ticks_ = min_ticks;
		//	return *this;
		// }


		// Derived& x_steps(int steps)
		// {	//! Set autoscale to set ticks in steps multiples of:\n
		//	 //! 2,4,6,8,10, if 2\n
		//	 //! or 1,5,10 if 5\n
		//	 //! or 2,5,10 if 10.\n
		//	 //! default = 0 (none).
		//	 //! \note: Must \b preceed x_autoscale(data) call).
		//	x_steps_ = steps;
		//	return *this;
		// }




	  //Derived& limit_color(const svg_color& col)
	  //{ //! Set the color for 'at limit' point stroke color.		  Need to set the series
		 //g(g_LIMIT_POINTS).style().stroke_color(col);
		  //////serieses_[0].limit_point_color(col); // Would require to add some data first!
		 //return *this;
	  //}



	  //Derived& limit_fill_color(const svg_color& col)
	  //{ //! Set the color for 'at limit' point fill color.
		 //g(g_LIMIT_POINTS).style().fill_on(true);
		 //g(g_LIMIT_POINTS).style().fill_color(col);
		 //////serieses_[0].limit_point_style_.fill_color(col);
		 //return *this;
	  //}



//#endif // _PLOT_CHAINING_;


      //inline svplot&   y_ticks_right_on(bool cmd)
      //{ //! Set true if ticks on the Y-axis are to be on right of axis line.
        //y_ticks_.right_ticks_on_ = cmd;	 return *this;
      //}



		 // Get the minimum and maximum (cartesian data units) for the plot window axes.


//
//		 Derived& x_label_on(bool cmd)
//		 { //! Set true if want to show X-axis label text.
//			//! \details Also switched on by setting label text.
//			//! (on the assumption that if label text is set, display is also wanted,
//			//! but can be switched off if \b not required).
//			x_ax_.label_on_ = cmd;
//			return *this;
//		 }
//
//		 Derived& x_data_font_size(unsigned int i)
//		 { //! Set X tick value label font size (svg units, default pixels).
//			x_data_value.textstyle().font_size(i);
//			return *this;
//		 }


		// Derived& x_tick_values_font_size(unsigned int i)
		// { //! Set X ticks value label font size (svg units, default pixels).
		//	x_tick_.textstyle().font_size(i);
		//	return *this;
		// }



		// Derived& x_tick_values_font_family(const std::string& family)
		// { //! Set X ticks value label font family.
		//	x_tick_.textstyle().font_family(family);
		//	return *this;
		// }

//
//		 Derived& x_ticks_location(int cmd)
//		 { //!  Set  X ticks on window or axis
//			//! \arg cmd -1 bottom of plot window,
//			//! \arg cmd 0 on X axis.
//			//! \arg cmd +1 top of plot window.
//		
//			x_ticks_.location_ = cmd;
//			return *this;
//		 }



		 //Derived& axes_on(bool is)
		 //{ //! If set true, draw \b both x and y axes (note plural axes).
			//x_ax_.axis_line_on() = is;
			//y_ax_.axis_line_on() = is;
			//return *this;
		 //}

		 ////////////////////////////////


		 //Derived& x_axis_on(bool is)
		 //{ //! If set true, draw a horizontal X-axis line.
			//x_ax_.axis_line_on() = is;
			//return *this;
		 //}



		 //Derived& y_axis_on(bool is)
		 //{ //! If set true, draw a vertical Y-axis line.
			//y_ax_.axis_line_on() = is;
			//return *this;
		 //}


		 // enums like g_TITLE provide a std:string like "title"
		 // colors .stroke_color, .stroke_width and font are set in the appropriate g_element.


		// Derived& x_autoscale(bool b)
		// { //! Set @c true if to use autoscaled values for X-axis.
		//	 //if (b && x_ticks_.major_ticks_fraction_ < 0)
		//	 //{ // No autoscale values have been calculated, so not safe to make x_autoscale true.
		//	//	 throw std::runtime_error("X autoscale has not been calculated yet!" );
		//	 //}
		//	x_autoscale_ = b;
		//	return *this; //! \return Reference to caller to make chainable.
		// }


		// Derived& autoscale(bool b)
		// { //! Set whether to use X autoscaled values.
		//	//! Same as x_autoscale, and used by boxplot too.
		//	 //if (x_ticks_.major_ticks_fraction_ < 0)
		//	 //{ // No autoscale values have been calculated, so not safe to make x_autoscale true.
		//	//	 throw std::runtime_error("X-axis autoscale has not been calculated yet!" );
		//	 //}
		//	x_autoscale_ = b;
		//	return *this;
		// }

		//  Derived& confidence(double alpha)
		//  { //! Set alpha for displaying confidence intervals.
		// 	//! Default is 0.05 for 95% confidence.
		// 	if (alpha <= 0.)
		// 	{ // Warn and leave alpha_ unchanged.
		// 		std::cout << "alpha must be > 0." << std::endl;
		// 	}
		// 	else if (alpha > 0.5)
		// 	{ // Warn and leave alpha_ unchanged.
		// 		std::cout << "alpha must be fraction < 0.5 (for example, 0.05 for 95% confidence)" << std::endl;
		// 	}
		// 	else
		// 	{
		// 	  alpha_ = alpha;
		// 	}
		// 	return *this;
		//  }

//
//		 Derived& x_labels_strip_e0s(bool cmd)
//		 { //! Set if to strip redundant zeros, signs and exponents, for example, reducing "1.2e+000" to "1.2"
//			//! This markedly reduces visual clutter, and is the default.
//			x_ticks_.strip_e0s_ = cmd;
//			return *this;
//		 }

		//
		//Derived& legend_lines(bool is)
		//{ //! Set true if legend should include samples of the lines joining data points.
		//	//! \details This allows different series of data points to be distinguished by different color and/or width.
		//	//! This is especially useful to show plots of different functions and/or different parameters in different colors.
		//	legend_lines_ = is;
		//	return *this;
		//}

//
//		 Derived& x_num_major_ticks(double inter)
//		 { //! Set the interval between X-axis major ticks.
//			x_ticks_.major_ticks_fraction_ = inter;
//			return *this;
//		 }

 //     inline svplot&   x_tick_values_side(int side)
 //     { //! Set which side (up, down or none) for major ticks label values:
 //       //! \param side -1 labels downward, 0 no labels, +1 labels upward.
 //      
 //       x_ticks_.tick_values_side_ = side;
 //       return *this; .
 //     }



		 //
		 //Derived& x_tick_values_ioflags(std::ios_base::fmtflags f)
		 //{ //! Set iostream format flags of data point X values near data points markers.//! Useful to set hexadecimal, fixed and scientific, (std::ios::scientific).
			//x_ticks_.value_ioflags_ = f;
			//return *this;
		 //}


		 //
		 //Derived& x_size(unsigned int i)
		 //{ //! Set SVG image X-axis size (SVG units, default pixels).
		//	svchart::x_size(i);
		//	return *this;
		 //}


//
//		 Derived& image_x_size(unsigned int i) //!< Obselete - deprecated.
//		 { //! Set SVG image X-axis size (SVG units, default pixels).
//			// Can't this be x_size(unsigned int i)
//			imageContainer_.x_size(i);
//			return *this;
//		 }
//
//
//		 unsigned int image_y_size() //!< Obselete - deprecated.
//		 { //! \return SVG image Y-axis size as vertical height (SVG units, default pixels).
//			return imageContainer_.y_size();
//		 }
//
//
//		 Derived& image_y_size(unsigned int i) //!< Obselete - deprecated.
//		 {//! Set SVG image Y-axis size (SVG units, default pixels).
//			imageContainer_.y_size(i);
//			return *this;
//		 }


//
//		 Derived& description(const std::string d)
//		 { //! Writes description to the document for header, for example:
//			//! \verbatim
//			//!  <desc> My Data </desc>
//			//! \endverbatim
//		
//			image_.description(d);
//			return *this;
//		 }
//
//
//
//
//		 Derived& document_title(const std::string d)
//		 { //! Write document title to the SVG document for header as \verbatim <title> My Title </title>  \endverbatim
//			image_.document_title(d);
//			return *this;
//		 }

//
//		 Derived& license(std::string repro,
//			std::string distrib,
//			std::string attrib,
//			std::string commercial,
//			std::string derivative)
//		 { //! Set license conditions for reproduction, attribution, commercial use, and derivative works,
//			//! usually "permits", "requires", or "prohibits",
//			//! and set license_on == true.
//			// Might check these are "permits", "requires", or "prohibits"?
//			image_.license(repro, distrib, attrib, commercial, derivative);
//			return *this;
//		 }
//
//
//		 Derived&  license_on(bool l)
//		 { //! Set if license conditions should be included in the SVG document.
//			//! \see axis_pl_frame::license
//			image_.license_on(l);
//			return *this;
//		 }
//
//		 Derived& boost_license_on(bool l)
//		 { //! Set if the Boost license conditions should be included in the SVG document.
//			image_.boost_license_on(l);
//			return *this;
//		 }
//
//
//
//		 Derived& coord_precision(int digits)
//		 { //! Precision of SVG coordinates in decimal digits (default 3).
//			//!  3 decimal digits precision is sufficient for small images.
//			//!  4 or 5 decimal digits precision will give higher quality plots,
//			//!  especially for larger images, at the expense of larger .svg files,
//			//!  particularly if there are very many data points.
//		
//			image_.coord_precision(digits);
//			return *this;
//		 }
* 
*/

//class pair_double_2d_convert
//{ //! \class pair_double_2d_convert
//  //! \brief This functor allows any 2 D data convertible to type std::pair<double, double> to be plotted.
// public:
//    typedef std::pair<double, double> result_type; //!< result type is a pair (X and Y) of doubles.
//
//    double i; //!< Current value, 1st set by start(double i0).
//
//    void start(double i0)    {  i = i0;  } //! Set a start value.
//
//
//     //! Convert a pair of X and Y (whose types can be converted to double values) to a pair of doubles.
//     //! \tparam T type whose value can be converted to double.
//     //! \tparam U type whose value can be converted to double.
//    template <typename T, typename U>
//    std::pair<double, double> operator()(const std::pair<T, U>& a) const
//    { //! Assumes that a conversion from double yields just the value component of the uncertain value.
//        return std::pair<double, double>((double)(a.first), (double)(a.second));
//    }
//
//    template <typename T>
//    std::pair<double, double> operator()(T a)
//    {  //! Convert a pair of X and Y values to a pair of doubles.
//        return std::pair<double, double>(i++, (double)a); //! \return pair of doubles.
//    }
//}; // class pair_double_2d_convert
//

	//template <typename T, typename U, typename std::enable_if<has_begin<T>::value,bool>::type = 0>
	//svplot_series& plot(const T& container, const std::string& title = "", U functor = pair_double_2d_convert() );

	//template <typename T, typename U, typename std::enable_if<svplot::is_ptr<T>::value,bool>::type = 0>
	//svplot_series& plot(const T& begin, const T& end, const std::string& title = "", U functor = pair_double_2d_convert() );


  ////!Add a container of a data series to the plot.\n
  ////!  This version permits a custom functor (rather than default conversion to @c double).\n
  ////!  \note that this version assumes that @b ALL the data values in the container is used.
  ////!  template <typename T, typename U, typename std::enable_if<svplot::has_begin<T>::value,bool>::type>
  //svplot_series& svplot::plot(const T& container, const std::string& title /* = "" */, U functor /* = pair_double_2d_convert*/)
  //{
    //serieses_.push_back(
      //svplot_series(
      //boost::make_transform_iterator(container.begin(), functor),
      //boost::make_transform_iterator(container.end(),   functor),
      //title, auto_point_shape(), auto_color(), auto_stroke_dash())
    //);
    //return serieses_[serieses_.size()-1]; //! \return Reference to data series just added to make chainable.
  //}
// //! Add (part of) a container of a data series to the plot, using a functor.
//      //!This version permits part of the container to be used, a partial range, using iterators begin to end.\n
//      //!Version with custom functor, rather than automatically converting to double).
  //template <typename T, typename U, typename std::enable_if<svplot::is_ptr<T>::value,bool>::type>
  //svplot_series& svplot::plot(const T& begin, const T& end, const std::string& title, U functor)
  //{
    //serieses_.push_back(
      //svplot_series(
      //boost::make_transform_iterator(begin, functor),
      //boost::make_transform_iterator(end,   functor),
      //title, auto_point_shape(), auto_color(), auto_stroke_dash())
    //);
    //return serieses_[serieses_.size() - 1]; //! \return Reference to data series just added to make chainable.
  //}







  // http://www.croczilla.com/~alex/conformance_suite/svg/text-align-02-b.svg
  // tests for baseline shifted text.  This is needed for subscript and superscript,
  // vital for nice display of units like m^2 and chemical formulae like H2O
  // IE (Adobe SVG viewer) and Opera conforms but not Firefox (yet).

  //// operators needed for testing at least.
  //bool operator==(const text_style& ts)
  //{ //! Compare text_style for equality (needed for testing).
  // return (ts.font_size_ == font_size_)
  //   && (ts.font_family_ == font_family_)
  //   && (ts.stretch_ == stretch_)
  //   && (ts.style_ == style_)
  //   && (ts.weight_ == weight_)
  //   && (ts.decoration_ == decoration_);
  //} // operator==
  //
  //bool operator!=(const text_style& ts)
  //{ //! Compare text_style for inequality (needed for testing).
  // return (ts.font_size_ != font_size_)
  //   || (ts.font_family_ != font_family_)
  //   || (ts.stretch_ != stretch_)
  //   || (ts.style_ != style_)
  //   || (ts.weight_ != weight_)
  //   || (ts.decoration_ != decoration_);
  //} //  operator!=
  //
  //bool operator==(const text_style& rhs) const
  //{ //! Compare two text_style for equality
  //  //! Note operator== and operator << both needed to use Boost.text.
  //  //! (But can be avoided with a macro define).
  //   return (font_size_ == rhs.font_size_)
  //     && (font_family_ == rhs.font_family())
  //     && (stretch_ ==  rhs.stretch_)
  //     && (style_ ==  rhs.style_)
  //     && (weight_ ==  rhs.weight_)
  //     && (decoration_ ==  rhs.decoration_);
  //} //   bool operator==(const text_style& lhs, const text_style& rhs)
  //
  //bool operator!= (const text_style& rhs) const
  //{ //! Compare two text_style for equality.
  //  //! Note operator== and operator << both needed to use Boost.Test.
  //  //! (But can be avoided with a macro define).
  //    return (font_size_ != rhs.font_size_)
  //     && (font_family_ != rhs.font_family())
  //     && (stretch_ !=  rhs.stretch_)
  //     && (style_ !=  rhs.style_)
  //     && (weight_ !=  rhs.weight_)
  //     && (decoration_ !=  rhs.decoration_);
  //} //   bool operator!= (const text_style& lhs, const text_style& rhs)
