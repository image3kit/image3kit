
// Copyright Paul A. Bristow 2006 - 2013.
// Copyright Ali Q. Raeini 2019

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)


svplot& title(const std::string title)  {
	//! Set a title for plot.
	//! The string may include Unicode for greek letter and symbols.
	//! @b example:
	//! A title that includes a greek omega and degree symbols:
	//! \code
	//!	 my_plot.title("Plot of &#x3A9; function (&#x00B0;C)");
	//! \endcode
	//! Unicode symbols are at http://unicode.org/charts/symbols.html.
	// Plot title.  TODO
	// new text parent code push_back
	// effectively concatenates with any existing title.
	// So clear the existing string first but doesn't work,
	// so need to clear the whole g_element.
	//g(g_TITLE).clear();
	title_info_.text(title);
	return *this;
 }


svplot& title_font_size(unsigned int i)  { //! Sets the font size for the title (svg units, default pixels).
	title_info_.textstyle().font_size(i);
	return *this;
 }


svplot& legend_on(bool cmd)  { //! Set true if a legend is wanted.
	legend_on_ = cmd;
	return *this;
 }


svplot& legend_place(legend_places l)  { //! Set the position of the legend, \see  svg::legend_places
	legend_place_ = l;
	return *this;
 }


svplot& legend_title_font_size(unsigned int size)  { //! Set the font size for the legend title (svg units, default pixels).
	legend_header_.textstyle().font_size(size);
	return *this;
 }






svplot& x_tick_values_precision(int p)  { //! Set iostream decimal digits precision of data point X values near data points markers.
	x_ticks_.value_precision_ = p;
	return *this;
 }






//svplot& size(unsigned int x, unsigned int y)  { //! Set SVG image size (SVG units, default pixels).
	//setsize(x, y);			// Check on sanity of these values?
	//return *this;
 //}



svplot& background_color(const svg_color& col)  { //! Set plot background color.
	g(g_BACKGROUND).style().fill_color(col);
	return *this;
 }


svplot& background_border_color(const svg_color& col)  { //! Set plot background border color.
	g(g_BACKGROUND).style().stroke_color(col);
	/*!
	  background_border_color, for example:
	  \verbatim svplot my_plot(my_data, "My Data").background_border_color(red).background_color(azure);
	  \endverbatim
	*/
	return *this;
 }



svplot& background_border_width(double w)  { //! Set plot background border width.
	g(g_BACKGROUND).style().stroke_width(w);
	return *this;
 }




svplot& x_data_values_precision(int digits)  { /*! Precision of X tick label values in decimal digits (default 3).
	  3 decimal digits precision is sufficient for small images.
	  4 or 5 decimal digits precision will give more cluttered plots.
	  If the range of labels is very small, then more digits will be essential.
	*/

	x_ticks_.value_precision_ = digits;
	return *this;
 }



svplot& x_data_ioflags(std::ios_base::fmtflags flags)  { //! Set iostream std::ios::fmtflags for X value label (default decimal == 0X201).
	//! Mainly useful for changing to scientific, fixed or hexadecimal format.
	//! For example: .x_data_ioflags(std::ios::dec | std::ios::scientific)

	x_ticks_.value_ioflags_ = flags;
	return *this;
 }


svplot& title_font_family(const std::string& family)  { //! Set the font family for the title (for example: .title_font_family("Lucida Sans Unicode");
	title_info_.textstyle().font_family(family);
	return *this;
 }


svplot& title_font_style(const std::string& style)  { //! Set the font style for the title (default normal).
	title_info_.textstyle().font_style(style);
	return *this;
 }



svplot& title_font_weight(const std::string& weight)  { //! Set the font weight for the title (default normal).
	title_info_.textstyle().font_weight(weight);
	return *this;
 }



svplot& title_font_stretch(const std::string& stretch)  { //! Set the font stretch for the title (default normal), wider or narrow.
	title_info_.textstyle().font_stretch(stretch);
	return *this;
 }



svplot& title_font_decoration(const std::string& decoration)  { //! Set the font decoration for the title (default normal, or underline, overline or strike-thru).
	title_info_.textstyle().font_decoration(decoration);
	return *this;
 }


svplot& title_font_rotation(rotate_style rotate)  { //! Set the rotation for the title font (degrees, 0 to 360).
	title_info_.rotation(rotate);
	return *this;
 }



svplot& title_font_alignment(align_style alignment)  { //! Set the alignment for the title.
	title_info_.alignment(alignment);
	return *this;
 }



svplot& legend_width(double width)  { //! Set the width for the legend.
	legend_width_ = width;
	return *this;
 }



svplot& legend_title(const std::string title)  { //! Set the title for the legend.
	legend_header_.text(title);
	return *this;
 }



svplot& legend_font_weight(const std::string& weight)  { //! Set the font weight for the legend title.
	legend_header_.textstyle().font_weight(weight);
	return *this;
 }


svplot& legend_font_family(const std::string& family)  { //! Set the font family for the legend title.
	legend_header_.textstyle().font_family(family);
	return *this;
 }


svplot& legend_top_left(double x, double y)  { //! Set position of top left of legend box (svg coordinates, default pixels).
	//! Bottom right is controlled by contents, so the user cannot set it.
	if((x < 0) || (x > x_size()) || (y < 0) || (y > y_size()))
	{
	  throw std::runtime_error("Legend box position outside image!");
	}
	legend_left_ = x;
	legend_top_ = y;
	return *this;
 }




svplot& legend_header_font_size(int size)  { //! Set legend header font size (svg units, default pixels).
	legend_header_.textstyle().font_size(size);
	return *this;
 }


svplot& plot_window_on(bool cmd)  { //! Set true if a plot window is wanted (or false if the whole image is to be used).
	pl_window_on_ = cmd;
	if(cmd)
	{ // Set plot window.
	  g(g_AREA).style()
		 .fill_color(area_.fill_) // background color and
		 .stroke_color(area_.stroke_); // border color.
	}
	//legend_place_ = outside_right;
	return *this;
 }


svplot& plot_border_color(const svg_color& col)  { //! Set the color for the plot window background.
	area_.stroke_ = col;
	g(g_AREA).style().stroke_color(col);
	return *this;
 }


svplot& plot_border_width(double w)  { //! Set the width for the plot window border (svg units, default pixels).
	area_.width_ = w;
	g(g_AREA).style().stroke_width(w);
	return *this;
 }


svplot& image_border_margin(double w)  { //! Set the margin around the plot window border (svg units, default pixels).
	//! \details This prevents the plot window getting too close to other elements of the plot.
	pl_box_.margin_ = w;
	return *this;
 }


svplot& image_border_width(double w)  { //! Set the svg image border width (svg units, default pixels).
	pl_box_.width_ = w;
	return *this;
 }


svplot& plot_window_x(double min_x, double max_x)  { //! Set the minimum and maximum (cartesian data units) for the plot window X axis.
	//! This is normally calculated from other plot values.
	if(max_x <= min_x)
	{
	  throw std::runtime_error("plot_window X: x_max_ <= x_min_");
	}
	if((max_x - min_x) < std::numeric_limits<double>::epsilon() * 1000)
	{ // Range too small to display.
	  throw std::runtime_error("plot_window X range too small!" );
	}
	pl_left_ = min_x;
	pl_right_ = max_x;
	return *this;
 }


svplot& plot_window_y(double min_y, double max_y)  { //! Set the minimum and maximum (cartesian data units) for the plot window Y axis.
	//! This is normally calculated from other plot values.

	if(max_y <= min_y)
	{
	  throw std::runtime_error("plot_window Y: y_max_ <= x_min_");
	}
	if(max_y <= min_y)
	{
	  throw std::runtime_error("plot_window Y range too small!");
	}
	pl_top_ = min_y;
	pl_bottom_ = max_y;
	return *this;
 }



svplot& x_label_font_size(unsigned int i)  { //! Set X axis label font size (svg units, default pixels).
	x_label_info_.textstyle().font_size(i);
	// Also duplicated at
	// x_label_style_.font_size(i);
	return *this;
 }



svplot& x_label_font_family(const std::string& family)  { //! Set X tick value label font family.
	x_label_info_.textstyle().font_family(family);
	return *this;
 }



svplot& x_axis_label_color(const svg_color& col)  { //! Set X axis label color.
	g(g_X_LABEL).style().fill_color(col);
	//g(g_X_LABEL).style().stroke_color(col);
	// Setting the stroke color produces fuzzy characters :-(
	// Set BOTH stroke and fill to the same color?
  return *this;
 }




svplot& x_tick_values_color(const svg_color& col)  { //! Set X axis tick value label color.
	// Set BOTH stroke and fill to the same color.
	g(g_X_TICKS_VALUES).style().fill_color(col);
	//g(g_X_TICK_VALUE_LABELS).style().stroke_color(col);
	// Setting the stroke color produces fuzzy characters :-(
	//x_ticks_.ticks_values_color = col;
	return *this;
 }



svplot& x_tick_values_rotation(rotate_style rot)  { //! Set rotation for X ticks major value labels. (Default horizontal).		\see rotate_style
	x_ticks_.label_rotation_ = rot;
	return *this;
 }



svplot& x_major_grid_on(bool is)  { //! If set true, will include a major X-axis grid.
	x_ticks_.major_grid_on_ = is;
	return *this;
 }


svplot& x_minor_grid_on(bool is)  { //! If set true, will include a minor X-axis grid.
	x_ticks_.minor_grid_on_ = is;
	return *this;
 }



svplot& title_color(const svg_color& col)  { //! Set the color of any title of the plot.
	// Function title_color could set both fill (middle) and stroke (outside),
	// but just setting fill if simplest,
	// but does not allow separate inside & outside colors.
	g(g_TITLE).style().fill_color(col);
	//g(g_TITLE).style().stroke_color(col);
	return *this;
 }




svplot& legend_color(const svg_color& col)  { //! Set the color of the title of the legend.
	// g(g_LEGEND_TEXT).style().stroke_color(col);
	g(g_LEGEND_TEXT).style().fill_color(col);
	return *this;
 }


svplot& legend_background_color(const svg_color& col)  { //! Set the background fill color of the legend box.
	legend_.fill(col);
	g(g_LEGEND_BOX).style().fill_color(col);
	return *this;
 }


svplot& legend_border_color(const svg_color& col)  { //! Set the border stroke color of the legend box.
	legend_.stroke(col);
	g(g_LEGEND_BOX).style().stroke_color(col);
	return *this;
 }


svplot& plot_background_color(const svg_color& col)  { //! Set the fill color of the plot window background.
	g(g_AREA).style().fill_color(col);
	return *this;
 }



svplot& x_axis_color(const svg_color& col)  { //! Set the color of the X-axis line.
	// Note only stroke color is set.
	g(g_X_AXIS).style().stroke_color(col);
	return *this;
 }



svplot& y_axis_color(const svg_color& col)  { //! Set the color of the Y-axis line.
	g(g_Y_AXIS).style().stroke_color(col);
	return *this;
 }


svplot& x_label_color(const svg_color& col)  { //! Set the color of X-axis label (including any units).
	// add fill as well PAB Oct 07
	g(g_X_LABEL).style().fill_color(col);
	g(g_X_LABEL).style().stroke_color(col);
	return *this;
 }


svplot& x_label_width(double width)  { //! Set the width (boldness) of X-axis label (including any units).
	//! (not recommended until browsers implement better).
	// width of text is effectively the boldness.
	g(g_X_LABEL).style().stroke_width(width);
	return *this;
 }


svplot& y_label_color(const svg_color& col)  { //! Set the color of Y-axis label (including any units).
	g(g_Y_LABEL).style().fill_color(col);
	g(g_Y_LABEL).style().stroke_color(col);
	return *this;
 }



svplot& x_major_tick_color(const svg_color& col)  { //! Set the color of X-axis major ticks.
	g(g_X_MAJOR_TICKS).style().stroke_color(col);
	return *this;
 }



svplot& x_minor_tick_color(const svg_color& col)  { //! Set the color of X-axis minor ticks.
	g(g_X_MINOR_TICKS).style().stroke_color(col);
	return *this;
 }


svplot& x_major_grid_color(const svg_color& col)  { //! Set the color of X-axis major grid lines.
	g(g_X_MAJOR_GRID).style().stroke_color(col);
	return *this;
 }



svplot& x_major_grid_width(double w)  { //! Set the width of X-axis major grid lines.
	g(g_X_MAJOR_GRID).style().stroke_width(w);
	return *this;
 }



svplot& x_minor_grid_color(const svg_color& col)  { //! Set the color of X-axis minor grid lines.
	g(g_X_MINOR_GRID).style().stroke_color(col);
	return *this;
 }


svplot& x_minor_grid_width(double w)  { //! Set the width of X-axis minor grid lines.
	g(g_X_MINOR_GRID).style().stroke_width(w);
	return *this;
 }



svplot& x_axis_width(double width)  { //! Set the width of X-axis lines.
	g(g_X_AXIS).style().stroke_width(width);
	return *this;
 }


svplot& data_lines_width(double width)  { //! Set the width of lines joining data points.
	g(g_DATA_LINES).style().stroke_width(width);
	return *this;
 }


svplot& x_values_on(bool b)  { //! \return true if values of X data points are shown (for example: 1.23).
	// (Want override xy_values_on that would otherwise cause overwriting).
	// So the last values_on setting will prevail.
	// But this is only defined in 2D
	//if(xy_values_on())
	//{ // Would be overwritten by XY pair.
	//  xy_values_on(false);
	//}
	x_values_on_ = b;
	return *this;
 }


svplot& x_values_font_size(unsigned int i)  { //! Set font size of data point X values near data points markers.
	x_values_style_.val_style_.font_size(i);
	return *this;
 }



svplot& x_values_font_family(const std::string& family)  { //! Set font family of data point X values near data points markers.
	x_values_style_.val_style_.font_family(family);
	return *this;
 }




svplot& x_num_minor_ticks(unsigned int num)  { //! Set number of X-axis minor ticks between major ticks.
	x_ticks_.num_minor_ticks_ = num;
	return *this;
 }


svplot& x_label(const std::string& str)  { //! Set the text to label the X-axis (and set x_label_on(true)).
	x_label_info_.text(str);
	// Might switch label_on false if null string?
	return *this;
 }


svplot& y_label(const std::string& str)  { //! Set the text for the Y-axis label (and set y_label_on(true)).
	y_label_info_.text(str);
	return *this;
 }

svplot& x_range(double min_x, double max_x)  { //! Set the range of values on the X-axis.
	//! The minimum and maximum values must be finite and not too near
	//! to the minima or maxima that can be represented by floating point doubles,
	//! and the range must not be too small.
	if (!std::isfinite(min_x))
	{
	  throw std::runtime_error("X range: min not finite!");
	}
	if (!std::isfinite(max_x))
	{
	  throw std::runtime_error("X range: max not finite!");
	}
	if(max_x < min_x)
	{ // max_x <= min_x.
	  std::stringstream message("X range: max <= min! ");
	  message << max_x << " <= " << min_x << std::ends;
	  throw std::runtime_error(message.str());
	  //throw std::runtime_error("X range: max <= min!");
	}
	if(  (std::abs(max_x - min_x) < std::numeric_limits<double>::epsilon() * 1000 * std::abs(max_x))
	  )
	{ // Range too small to display.
	  std::cerr<<"Warning: X range too small, readjusting!"<<std::endl;
	  max_x = 0.5*(max_x+min_x)+std::numeric_limits<double>::epsilon() * 500;
	  min_x = max_x - std::numeric_limits<double>::epsilon() * 1000;
	}
	x_ticks_.min_ = min_x;
	x_ticks_.max_ = max_x; 
	return *this;
 }



svplot& x_values_color(const svg_color& col)  { //! Set the color of data point X values near data points markers.
	// Function could set both fill (middle) and stroke (outside),
	// but just setting fill is simplest,
	// but does not allow separate inside & outside colors.
	// Might be better to set in x_values_style
	g(g_X_POINT_VALUES).style().fill_color(col);
	//g(g_X_POINT_VALUES).style().stroke_color(col);

	return *this;
 }


svplot& x_values_rotation(rotate_style rotate)  { //! \return  the rotation (rotate_style) of data point X values near data points markers.
	//! (Degrees: 0 to 360 in 45 steps).
	x_values_style_.value_label_rotation_ = rotate;
	return *this;
 }


svplot& x_values_precision(int p)  { //! Set iostream decimal digits precision of data point X values near data points markers.
	x_values_style_.value_precision_ = p;
	return *this;
 }


svplot& x_values_ioflags(std::ios_base::fmtflags f)  { //! Set iostream format flags of data point X values near data points markers.
	//! Useful to set hexadecimal, fixed and scientific, (std::ios::scientific).
	x_values_style_.value_ioflags_ = f;
	return *this;
 }




svplot& x_plusminus_on(bool b)  { //! Set if to append std_dev estimate to data point X values near data points markers.
	//! (May not be implemented yet).
	x_values_style_.plusminus_on_ = b;
	return *this;
 }



svplot& x_plusminus_color(const svg_color& col)  { //! Set the color of X std_dev of value, for example, the color of 0.02 in "1.23 +-0.02 (9)".
	x_values_style_.plusminus_color_ = col;
	return *this;
 }


svplot& x_addlimits_on(bool b)  { //! Set if to append confidence limits to data point X values near data points markers.
	//! (May not be implemented yet).
	x_values_style_.addlimits_on_ = b;
	return *this;
 }



svplot& x_addlimits_color(const svg_color& col)  { //! Set the color of X confidence limits value, for example, the color of "<1.23 , 1.34>".
	x_values_style_.addlimits_color_ = col;
	return *this;
 }



svplot& x_df_on(bool b)  { //! Set true if to append a degrees of freedom estimate to data point X values near data points markers.
	x_values_style_.df_on_ = b;
	return *this;
 }


 svg_color x_addlimits_color()  { //! Get the color of X confidence limits of value, for example, the color of "<1.23 , 1.34>"".
	return x_values_style_.addlimits_color_;
 }

 bool x_df_on()  { //! \return true if to append a degrees of freedom estimate to data point X values near data points markers.
	//! (May not be implemented yet).
	return x_values_style_.df_on_;
 }


 svg_color x_df_color()  { //! Get the color of X degrees of freedom, for example, the color of 9 in "1.23 +-0.02 (9)".
	return x_values_style_.df_color_;
 }


 bool x_id_on()  { //! \return true if to append an ID or name to data point X values near data points markers.
	//! (May not be implemented yet).
	return x_values_style_.id_on_;
 }


svplot& x_df_color(const svg_color& col)  { //! Set the color of X degrees of freedom, for example, the color of 9 in "1.23 +-0.02 (9)".
	x_values_style_.df_color_ = col;
	return *this;
 }


svplot& x_id_on(bool b)  { //! Set true if to append a id or name to data point X values near data points markers.
	//! (May not be implemented yet).
	x_values_style_.id_on_ = b;
	return *this;
 }



svplot& x_id_color(const svg_color& col)  { //! Set the color of X ID or name, for example, the color of text in "My_id".
	x_values_style_.id_color_ = col;
	return *this;
 }


svplot& x_datetime_on(bool b)  { //! Set true if to append a datetime to data point X values near data points markers.
	//! (May not be implemented yet).
	x_values_style_.datetime_on_ = b;
	return *this;
 }


svplot& x_datetime_color(const svg_color& col)  { //! Set the color of X point datetime, for example, the color of text in "2004-Jan-1 05:21:33.20".
	x_values_style_.datetime_color_ = col;
	return *this;
 }


svplot& x_order_on(bool b)  { //! Set true if to append an order # to data point X values near data points markers.
	x_values_style_.order_on_ = b;
	return *this;
 }



svplot& x_order_color(const svg_color& col)  { //! Set the color of X point order in sequence, for example, #3.
	x_values_style_.order_color_ = col;
	return *this;
 }



svplot& x_decor(const std::string& pre, const std::string& sep, const std::string& suf)  { /*! Set prefix, separator and suffix for x_style
	\note if you want a space, you must use a Unicode space "\&#x00A0;",
	for example, ",\&#x00A0;" rather than just ", ".
	*/
	x_values_style_.prefix_ = pre;
	x_values_style_.separator_ = sep;
	x_values_style_.suffix_ = suf;
	return *this;
 }



svplot& x_major_tick_length(double length)  { //! Set length of X major ticks.
	x_ticks_.major_tick_length_ = length;
	return *this;
 }



svplot& x_major_tick_width(double width)  { //! Set width of X major ticks.
	x_ticks_.major_tick_width_ = width; // Redundant?
	g(g_X_MAJOR_TICKS).style().stroke_width(width);
	return *this;
 }



svplot& x_minor_tick_length(double length)  { //! Set length of X minor ticks.
	x_ticks_.minor_tick_length_ = length;
	return *this;
 }



svplot& x_minor_tick_width(double width)  { //! Set width of X minor ticks.
	x_ticks_.minor_tick_width_ = width;
	g(g_X_MINOR_TICKS).style().stroke_width(width);
	return *this;
 }




svplot& autoscale_check_limits(bool b)  { //! Set to check that values used for autoscale are within limits.
	//! Default is true, but can switch off checks for speed.
	autoscale_check_limits_ = b;
	return *this;
 }


svplot& autoscale_plusminus(double pm)  { //! Set how many std_dev or standard deviation to allow for ellipse when autoscaling.
	//! Default is 3 for 99% confidence.
	autoscale_plusminus_ = pm;
	return *this;
 }



svplot& x_with_zero(bool b)  { //! Set X-axis autoscale to include zero (default = false).
	//! Must preceed x_autoscale(data) call.
	x_include_zero_ = b;
	return *this;
 }

svplot& x_tight(double tight)  { //! Set tolerance to autoscale to permit data points slightly outside both end ticks. default 0. Must preceed x_autoscale(data) call.
	x_tight_ = tight;
	return *this;
 }




svplot& y_values_on(bool b)  { //! Set @c true if values of Y data points are shown (for example: 1.23, 2.34 or (x,y) if x_values_on is set too).
  y_values_on_ = b;	 return *this;
}


svplot& y_plusminus_on(bool b)  { //! Set true if values of Y data points are to include uncertainty estimates.
  y_values_style_.plusminus_on_ = b;
  return *this;
}

svplot& y_plusminus_color(const svg_color& col)  {  //! Set color of Y uncertainty of value.
  y_values_style_.plusminus_color_ = col;	 return *this;
}

svplot& y_addlimits_on(bool b)  { //! Set true if values of Y data points are to include confidence interval.
  y_values_style_.addlimits_on_ = b;	 return *this;
}

svplot& y_addlimits_color(const svg_color& col)  {  //! Set color of Y confidence interval.
  y_values_style_.addlimits_color_ = col;	 return *this;
}

svplot& y_df_on(bool b)  { //! Set @c true if values of Y data points are to include degrees of freedom estimates.
  y_values_style_.df_on_ = b;	 return *this;
}

svplot& y_df_color(const svg_color& col)  {  //! Set color of Y degrees of freedom.
  y_values_style_.df_color_ = col;	 return *this;
}

/*! \brief Set prefix, separator and suffix for Y-axis.\n @b Example:
\code        
	my_1d_plot.x_decor("[ x = ", "", "&#x00A0;sec]");
\endcode
 \note If you want a space, you must use a @b Unicode space.
*/    
svplot& y_decor(const std::string& pre, const std::string& sep, const std::string& suf)  {
  y_values_style_.prefix_ = pre;
  y_values_style_.separator_ = sep;
  y_values_style_.suffix_ = suf;	 return *this;
}



svplot& y_major_tick_length(double length)  { //! Set major tick length on Y-axis.
  y_ticks_.major_tick_length_ = length;	 return *this;
}


svplot& y_minor_tick_length(double length)  { //! Set minor tick length on Y-axis.
  y_ticks_.minor_tick_length_ = length;	 return *this;
}

svplot& y_num_minor_ticks(unsigned int num)  { //! Set number of minor ticks on Y-axis.
  y_ticks_.num_minor_ticks_ = num;	 return *this;
}

svplot& y_label_axis(const std::string& str)  { //! Set text to label Y-axis.
  y_label_info_.text(str);	 return *this;
}

svplot& y_major_tick_width(double width)  { //! Set width of major ticks on Y-axis.
  y_ticks_.major_tick_width_ = width;
  Self_.g(g_Y_MAJOR_TICKS).style().stroke_width(width);
  return *this; //! \return reference to svplot to make chainable.
}

svplot& y_minor_tick_width(double width)  { //! Set width of minor ticks on Y-axis.
  y_ticks_.minor_tick_width_ = width;
  Self_.g(g_Y_MINOR_TICKS).style().stroke_width(width);	 return *this;
}


svplot& x_axis_location(axis_intersect pos)  {  /*! Set position of the horizontal X-axis line (on the border).\n
	But controlled by the intersection with Y-axis,
	so this only changes the default position from bottom to top,
	but will be changed if X-axis intersects the Y-axis
	(that is, if Y-axis includes zero).     */
 x_ax_.location_ = pos; // top or bottom
 return *this;  //! \return Reference to svg_boxplot to make chainable.
}

svplot& y_axis_location(axis_intersect pos)  { /*! Set position of the vertical Y-axis line (on the border).
	But controlled by the intersection with X-axis,
	so this only changes the default position from bottom to top,
	but will be changed if X-axis intersects the X-axis
	(that is if X-axis includes zero).    */
 y_ax_.location_ = pos; // left or right
 return *this; //! \return Reference to svg_boxplot to make chainable.
}

svplot& x_ticks_location(tick_place side)  { // Set if ticks on the plot window or on the X-axis.
  //! \param side -1 ticks downward, 0 no ticks, +1 ticks upward.
  x_ticks_.ticks_loc_ = side;
  return *this; //! \return reference to svplot to make chainable.
}

svplot& y_ticks_location(tick_place cmd)  { //!  Set Y ticks on window or axis
	//! \param cmd -1 left of plot window, 0 on Y-axis, +1 right of plot window.
  y_ticks_.ticks_loc_ = cmd;
  return *this; //! \return reference to svplot to make chainable.
}

svplot& x_tick_values_side(tick_side cmd)  { //! Set true if ticks on the Y-axis are to be on left of axis line.
  x_ticks_.values_side_ = cmd;	 return *this;
}

svplot& y_tick_values_side(tick_side side)  { //! Position of labels for major ticks on vertical Y-axis line. 
  y_ticks_.values_side_ = side;	 return *this;
}



svplot& y_ticks_sides(tick_side cmd)  { //! Set true if X major ticks should mark upwards.
	y_ticks_.ticks_sides_ = cmd;
	return *this;
 }


svplot& x_ticks_sides(tick_side cmd)  { //! Set true if X major ticks should mark downwards.
	x_ticks_.ticks_sides_ = cmd;
	return *this;
 }



svplot& y_major_grid_on(bool is)  { //! Set true to include major grid lines.
  y_ticks_.major_grid_on_ = is;	 return *this;
}

svplot& y_minor_grid_on(bool is)  {//! Set true to include minor grid lines.
  y_ticks_.minor_grid_on_ = is;	 return *this;
}

svplot& y_minor_grid_width(double width)  { //! Set width of minor grid lines.
  y_ticks_.minor_grid_width_ = width;
  Self_.g(g_Y_MINOR_GRID).style().stroke_width(width);
  return *this; //! \return reference to svplot to make chainable.
}

svplot& y_major_grid_width(double width)  { //! Set width of major grid lines.
  y_ticks_.major_grid_width_ = width;
  Self_.g(g_Y_MAJOR_GRID).style().stroke_width(width);	 return *this;
}

svplot& y_label_font_size(unsigned int i)  { //! Set Y-axis label text font size.
  // May be best to tie label & unit font sizes together?
  // y_label_style_.font_size(i);
  y_label_info_.textstyle().font_size(i);
  return *this;
}

svplot& y_label_weight(std::string s)  { //! Set Y-axis label text font weight (for example: "bold").
  //! ("bold" is only one that works so far, and quality may be poor for some browsers).
  //y_label_style_.font_weight(s);
  y_label_info_.textstyle().font_weight(s);	 return *this;
}

svplot& y_label_font_family(const std::string& family)  { /*! Set Y-axis label text font family (for example: "Lucida Sans Unicode").
	 Available fonts depend on the program rendering the SVG XML, usually a browser.
	 The default font (usually "Lucida Sans Unicode") is used
	 if a renderer (in a browser or a converter to PDF like RenderX)
	 does not provide the font specified.
	 A Unicode font has a better chance of providing Unicode symbols, for example, specified as @c \&\#x221E;.
	 These fonts are probably usable:
	 \code
		"arial", "impact", "courier", "lucida console",  "Lucida Sans Unicode", "Verdana", "calibri", "century",
		"lucida calligraphy", "tahoma", "vivaldi", "informal roman", "lucida handwriting", "lucida bright", "helvetica"
	 \endcode
  */
  //y_label_style_.font_family(family);
  y_label_info_.textstyle().font_family(family);	 return *this;
}



svplot& y_tick_values_color(const svg_color& col)  { //! Set color for Y_axis tick values.        // Function could set both fill (middle) and stroke (outside), but just setting fill if simplest,  but does not allow separate inside & outside colors.
  y_ticks_.values_color_ = col;
  Self_.g(g_Y_TICKS_VALUES).style().fill_color(col);	 return *this;
}


svplot& y_data_values_precision(int p)  { //! Set @c std::iostream decimal digits precision of ticks Y values.
	y_ticks_.value_precision_ = p;	 return *this;
}

svplot& y_tick_values_ioflags(std::ios_base::fmtflags f)  { //! Set @c std::iostream format flags of ticks Y values.
  //! Useful to set hexadecimal, fixed and scientific, (std::ios::scientific).
  y_ticks_.value_ioflags_ = f;
  return *this;
}

svplot& y_tick_values_font_size(unsigned int i)  { //! Set font size for Y-axis ticks values (svg units, default pixels).
  y_ticks_.value_label_style_.font_size(i);	 return *this;
}


svplot& y_tick_values_font_family(const std::string& family)  { /*! Set font family for Y-axis ticks values.
	 Available fonts depend on the program rendering the SVG XML, usually a browser.
	The default font (usually "verdana") is used if a render program does not provide the font specified.

	 These are probably usable:
  */

  /*!        
	 \code
"arial", "impact", "courier", "lucida console",  "Lucida sans unicode", "verdana", "calibri", "century",
"lucida calligraphy", "tahoma", "vivaldi", "informal roman", "lucida handwriting", "lucida bright", "helvetica"
	 \endcode
  */

  y_ticks_.value_label_style_.font_family(family);	 return *this;
}

svplot& y_values_font_size(unsigned int i)  { //! Set font size for Y-axis values.
  y_values_style_.val_style_.font_size(i);	 return *this;
}

svplot& y_values_font_family(const std::string& family)  { //! Set font family for Y-axis values.
  //! Available fonts depend on the program rendering the SVG XML, usually a browser.
  //! The default font (usually "verdana") is used if a render program does not provide the font specified.
  //! These are probably usable:
  //! "arial", "impact", "courier", "lucida console",  "Lucida sans unicode", "verdana", "calibri", "century",
  //! "lucida calligraphy", "tahoma", "vivaldi", "informal roman", "lucida handwriting", "lucida bright", "helvetica"
  y_values_style_.val_style_.font_family(family);	 return *this;
}


svplot& y_values_color(const svg_color& col)  { //! Set color for Y-axis values.
  // Function could set both fill (middle) and stroke (outside),
  // but just setting fill if simplest,
  // but does not allow separate inside & outside colors.
  Self_.g(g_Y_POINT_VALUES).style().fill_color(col);
  //svplot().Self_.g(g_Y_POINT_VALUES).style().stroke_color(col);
  return *this; //! \return reference to svplot to make chainable.
}

svplot& y_values_rotation(rotate_style rotate)  { //! Set rotation for value labels on Y-axis ticks.  \see rotate_style
  y_values_style_.value_label_rotation_ = rotate;	 return *this;
}

svplot& y_values_precision(int p)  { //! Set @c iostream precision for data points Y values.
  y_values_style_.value_precision_ = p;	 return *this;
}


svplot& y_values_ioflags(std::ios_base::fmtflags f)  { //! Set @c iostream format flags for data point values.
  y_values_style_.value_ioflags_ = f;
  return *this; //! \return reference to svplot to make chainable.
}




svplot& y_tick_values_rotation(rotate_style rot)  { /*! Rotation or orientation of labels for major ticks on vertical Y-axis line.
  \param rot Default orientation is horizontal.
  \see @c rotate_style for possible values: horizontal, uphill...
  */
  y_ticks_.label_rotation_ = rot;	 return *this;
}

svplot& y_axis_width(double width)  { //! Set width of Y-axis line.
  Self_.g(g_Y_AXIS).style().stroke_width(width);	 return *this;
}

svplot& y_tick_values_precision(int digits)  { //! Set precision of Y tick label values in decimal digits (default 3).
  //! @b Example:
	//! \code
	//!   my_plot.x_data_ioflags(ios::dec | ios::scientific).x_tick_values_precision(2);
	//! \endcode

  y_ticks_.value_precision_ = digits;	 return *this;
}


svplot& y_data_ioflags(std::ios_base::fmtflags flags)  { //! Set std::ioflags of Y tick label values (default 0x201 == dec).
  //! @b Example:
  //!
	//! \code
	//!	my_plot.x_data_ioflags(ios::dec | ios::scientific).x_tick_values_precision(2);
	//! \endcode
  y_values_style_.value_ioflags_ = flags;	 return *this;
}

svplot& y_labels_strip_e0s(bool cmd)  { //! If @c true then strip unnecessary zeros, signs from labels.
  y_ticks_.strip_0es_ = cmd;	 return *this;
}


svplot& y_axis_label_color(const svg_color& col)  { //! Set y_axis stroke color.
  //! \note Setting the stroke color may produce fuzzy characters :-(
  Self_.g(g_Y_LABEL).style().fill_color(col);	 return *this;
}

svplot& y_axis_value_color(const svg_color& col)  { //! Set color of Y-axis @b value labels.
  Self_.g(g_Y_TICKS_VALUES).style().stroke_color(col);	 return *this;
}

svplot& y_label_width(double width)  { //! Set width of Y-axis value labels.
  Self_.g(g_Y_LABEL).style().stroke_width(width);	 return *this;
}

svplot& y_major_grid_color(const svg_color& col)  { //! Set color of Y major grid lines.
  Self_.g(g_Y_MAJOR_GRID).style().stroke_color(col);	 return *this;
}

svplot& y_minor_grid_color(const svg_color& col)  { //! Set color of Y minor grid lines.
  Self_.g(g_Y_MINOR_GRID).style().stroke_color(col);	 return *this;
}

svplot& y_major_tick_color(const svg_color& col)  { //! Set color of Y major tick lines.
  Self_.g(g_Y_MAJOR_TICKS).style().stroke_color(col);	 return *this;
}

svplot& y_minor_tick_color(const svg_color& col)  {  //! Set color of Y minor tick lines.
  Self_.g(g_Y_MINOR_TICKS).style().stroke_color(col);	 return *this;
}

svplot& y_range(double min_y, double max_y)  { //! Set the range (max and min) for Y-axis from the parameters provided.
  if (!std::isfinite(min_y))
  {
	 throw std::runtime_error("Y range: min not finite!");
  }
  if (!std::isfinite(max_y))
  {
	 throw std::runtime_error("Y range: max not finite!");
  }

  if(max_y <= min_y)
  { // max <= min.
	 throw std::runtime_error("Y range: y max <= y min!");
  }
  if((max_y - min_y) < std::numeric_limits<double>::epsilon() * 1000)
  { // Range too small to display.
	 throw std::runtime_error("Y range too small!" );
  }
  y_ticks_.min_ = min_y;
  y_ticks_.max_ = max_y;
  //y_autoscale_ = false;
 return *this;
}



