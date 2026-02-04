
// Copyright Paul A. Bristow 2006 - 2013.
// Copyright Ali Q. Raeini 2019

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#define thisAxes   y_ax_


void draw_y_axis()
{ //! Draw the Y-axis line, grids and ticks with labels.	// Perhaps to left or right of plot window if X values do not include zero.

	const auto& tiks = y_ticks_;
	if (thisAxes.axis_line_on()) // Want a Y-axis line.
	{
		double y_1 = pl_bottom_;
		double y_2 = pl_top_;
		if (y_ax_.location_ & at_zero) // Draw the vertical Y-axis line at cartesian xx = 0).
		{
			if (x_ticks_.ticks_side()&down_side)
			{
				if (x_ticks_.ticks_loc_ < 0) // at bottom  // Extend the vertical line down in lieu of longest tick.
				  y_1 += std::max(x_ticks_.minor_tick_length_, x_ticks_.major_tick_length_);// Avoid macro max trap!
			}
			else if (x_ticks_.ticks_side()&up_side)
			{
				if (x_ticks_.ticks_loc_ > 0) // at top  // Extend the vertical line up in lieu of longest tick.
					y_2 += std::max(x_ticks_.minor_tick_length_, x_ticks_.major_tick_length_);// Avoid macro max trap!
			}
			double xx = thisAxes.position_; // Y-axis (xx = 0) transformed into X SVG coordinates.
			g(g_Y_AXIS).line(xx, y_1, xx, y_2);
			if (tiks.ticks_loc_ < 0) //(y_ax_.location_ == left)  // Draw vertical line holding the ticks on the left of plot window.
				g(g_Y_AXIS).line(pl_left_, y_1, pl_left_, y_2);
			else if (tiks.ticks_loc_ > 0)  // Draw vertical line holding the ticks on the right of plot window.
				g(g_Y_AXIS).line(pl_right_, y_1, pl_right_, y_2);
		}
		else if (y_ax_.location_ == at_min)  // Draw on the left of plot window.
			g(g_Y_AXIS).line(pl_left_, pl_top_, pl_left_, pl_bottom_);
		else if (y_ax_.location_ == at_max) // Draw on the right of plot window.
			g(g_Y_AXIS).line(pl_right_, pl_top_, pl_right_, pl_bottom_);
		else {   }// ??? Warn that things have gone wrong?
	} // axis_line_on()



	// Access the paths for the ticks & grids, ready for additions.
	path_element& minor_tick_path = g(g_Y_MINOR_TICKS).path();
	path_element& major_tick_path = g(g_Y_MAJOR_TICKS).path();
	path_element& minor_grid_path = g(g_Y_MINOR_GRID).path();
	path_element& major_grid_path = g(g_Y_MAJOR_GRID).path();

	// dMinor is the interval between minor ticks.
	double dMajor = (tiks.max_-tiks.min_)*tiks.major_ticks_frac_;
	double dMinor = dMajor / (double(tiks.num_minor_ticks_ + 1.) );
	// TODO Problem here with using floating point?
	// Was i < y_max_; but didn't show the tick and value at y_max_ so now i <= y_max_;
	// But may still fail if a ls or few bits out? Seems to fail for yy = 100, for example.

	// Draw the ticks on the positive side.
	for (double yy = y_include_zero_ ? 0. : std::max(0.,tiks.min_); yy <= tiks.max_; yy += dMajor)
	{
	  for(double j = yy + dMinor; j < (yy + dMajor) * (1. - 2*std::numeric_limits<double>::epsilon()); j += dMinor)
	  {
	    if (j != 0. || ! x_ax_.axis_line_on()) // Avoid a major tick at yy == 0 where there *is* a horizontal X-axis line.
	      draw_y_minor_tick(j, minor_tick_path, minor_grid_path);
	  }
	  // Draw major tick.
	  if ((yy != 0. || ! x_ax_.axis_line_on()) // axis line requested.
	    || (tiks.ticks_loc_ != 0)) // ticks & labels on plot window.
	  { // Avoid a major tick at yy == 0 where there *is* a horizontal X-axis line.
	    draw_y_major_tick(yy, major_tick_path, major_grid_path);
	  } //else  // std::cout << "Missed yy " << yy << std::endl; // Only miss 0s

	}

	// Draw the ticks on the negative side.
	for(double yy = y_include_zero_ ? 0.0 : std::min(0.0,tiks.max_); yy >= tiks.min_; yy -= dMajor)
	{
	  for(double j = yy; j > yy - dMajor; j-= dMajor / (tiks.num_minor_ticks_ + 1))
	  { // Draw minor ticks.
	    if ( (j != 0. || ! thisAxes.axis_line_on())
	      || (tiks.ticks_loc_ != 0) // ticks & labels on plot window.
	      )
	    { // Avoid a major tick at yy == 0 where there *is* a horizontal X-axis line.
	      draw_y_minor_tick(j, minor_tick_path, minor_grid_path);
	    }
	  }
	  if ((yy != 0. || ! x_ax_.axis_line_on())
	    || (tiks.ticks_loc_ != 0) ) // ticks & labels on plot window.
	  { // Avoid a major tick at yy == 0 where there *is* a horizontal X-axis line.
	    draw_y_major_tick(yy, major_tick_path, major_grid_path);
	  }
	}
} // void draw_axis






void draw_y_ax_label()
{ //! Draw a vertical Y-axis label, and optional yy units.
	const auto& tiks = y_ticks_;
	const  float  valFntSiz = y_ticks_.value_label_style_.font_size();
	const  float  lblFntSiz = y_label_info_.textstyle().font_size();
	rotate_style rot = tiks.label_rotation_;
	 // Glyphs for western characters are aligned with the left bottom of capital letter, so need to allow for any descenders.
	double xx = pl_left_; // left edge of plot window.
	if (tiks.ticks_loc_ < 0) // -1 means left
	{ // Ticks value labels left of plot window.
	  if (tiks.values_side_ < 0) // -1 means value label to left of Y-axis.
	  { // tick values labels are to left of Y-axis.,  Shift right to allow for any tick value labels.
	    if ((rot == downward) || (rot == upward))
	    { // downward tick value label direction 90 vertical up or down, or 60 steep degrees (might handle 60 separately).
	      if (tiks.values_side_ < 0) // tick value labels are to left of plot window.
	        xx -= valFntSiz * 1.3; // Allow space for tick value labels font size to left of Y-axis or plot window.
	      if (tiks.ticks_side()&left_side)  // Allow for any righttoleft ticks.
	        xx -= 1.1 * std::max(tiks.major_tick_length_, tiks.minor_tick_length_);// And avoid macro max trap!
	      xx -= 0.3 * (lblFntSiz + valFntSiz); // best compromise
	    }
	    else if ((rot == steepdown) || (rot == steepup))
	    { // downward tick value label direction 90 vertical up or down, or 60 steep degrees (might handle 60 separately).
	      if (tiks.values_side_ < 0) // tick value labels are to left of plot window.
	        xx -= valFntSiz * 1.3; // Allow space for tick value labels font size to left of Y-axis or plot window.
	      if (tiks.ticks_side()&left_side)  // Allow for any righttoleft ticks.
	        xx -= 1.1 * std::max(tiks.major_tick_length_, tiks.minor_tick_length_);// And avoid macro max trap!
	      xx -= 0.3 * (lblFntSiz + valFntSiz); // best compromise?
	    }
	    else if ((rot == uphill)  || (rot == downhill))
	    { // sloping 45 degrees .
	       xx -= tiks.label_max_space_ * sin45;
	       if (tiks.ticks_side()&left_side) // Move left for any righttoleft ticks, and a small space.
	         xx -= 0.3 * (lblFntSiz + valFntSiz); // best compromise?	        //xx -= 1.2 * (std::min)(lblFntSiz, valFntSiz ); // better

	    }
	    else if ((rot == slopeup) || (rot == slopedownhill))
	    { // sloping 30 degrees.
	       xx -= tiks.label_max_space_ * sin45;
	       if (tiks.ticks_side()&left_side)// Move left for any righttoleft ticks, and a small space.
	         xx -= 1.1 * std::max(tiks.major_tick_length_, tiks.minor_tick_length_) // And avoid macro max trap!
	               + 0.3 * (lblFntSiz + valFntSiz); // best compromise? //xx -= 1.2 * (std::min)(lblFntSiz, valFntSiz ); // better
	    }
	    else if  (rot == horizontal)
	    { // horizontal
	      if (tiks.ticks_side()&left_side)   // Move left for any righttoleft ticks, and a small space.
	        xx -= 1.1 * std::max(tiks.major_tick_length_, tiks.minor_tick_length_); // And avoid macro max trap!

	      xx -= tiks.label_max_space_; // Move left for the longest tick value label. (Might be zero?)
	      xx -= 0.3 * (lblFntSiz + valFntSiz); // best compromise?
	    }
	    else
	      std::cout << " Rotation of Y label" << rot << "not yet implemented!" << std::endl;
	  }
	  else if (tiks.values_side_ & right_side)// +1 means Y Tick labels to right of Y-axis.
	     xx -= lblFntSiz * 1.5;  // Move left from Y-axis.
	  else // tiks.y_tick_values_side_ == 0 means no tick value labels.
	    xx -= lblFntSiz * 1.5;  // Move left from Y-axis.


	  if (tiks.ticks_side()&left_side)  // Shift right for biggest of any righttoleft ticks.
	    xx += std::max(tiks.minor_tick_length_, tiks.major_tick_length_);
	}
	else if (tiks.ticks_loc_ > 0)
	{ // tick values labels are to right of Y-axis.
	   xx = 0. +  pl_box_.width_ // Start Y Label just right of the image left side.
	        + pl_box_.margin_ + lblFntSiz * 1.; // Shift right to suit Y labels.
	}
	else if (tiks.ticks_loc_ == 0)
	{  // Ticks are ON the Y-axis line, so Y label is just right the plot left.
	   // Character starts at bottom of capital letter, so allow for descenders.
	   xx = 0. +  pl_box_.width_; // Start Y Label just right of the image left side.
	   xx += pl_box_.margin_;
	   xx += lblFntSiz * 1.; // Shift right to suit Y labels.
	}
	// Glyph is at bottom left of western characters.

	g(g_Y_LABEL).text(xx, // distance from left side of image.
	   (pl_bottom_ + pl_top_) / 2., // center on the plot window.
	  y_label_info_.text(), // "Y-Axis" for example.
	  y_label_info_.textstyle(), // font and size
	  center_align, // One might want it to left or right_align?
	  upward); // Y label must be drawn vertically.

   } // draw_thisAxeslabel


void draw_y_major_tick(double value, path_element& tick_path, path_element& grid_path)
{ //! Draw a Y-axis major tick, tick value labels & grids.
	const  float  valFntSiz = y_ticks_.value_label_style_.font_size();
	const auto& tiks = y_ticks_;
	double yy(value); // for tick and/or grid.
	transform_y(yy); // Cartesian to SVG coordinates.
	if((yy < pl_top_ - 0.1) || (yy > pl_bottom_ + 0.1))  // Allow a bit extra to allow for round-off errors. tick value is way outside plot window, so nothing to do.
	{    (std::cout << " Yob ").flush();  return;     }


// TODO: merge these with draw_axis location computation
	double o_1(0.); // Left end of tick.
	double o_2(x_size()); // Right end of tick.
	if(tiks.major_grid_on_)
	{ // Draw horizontal major Y grid line.
	  if(!pl_window_on_)
	  {
	    if(tiks.values_side_ < 0) // left  // Start further right to give space for Y-axis value label.
	      yy -= valFntSiz * text_margin_;
	    if(tiks.ticks_side()&left_side)  // And similarly space for left ticks.
	      yy -= tiks.major_tick_length_;
	  }
	  else
	  { // pl_window_on_ to use full width of plot window.
	    o_1 = pl_left_ + area_.width_; // Don't write over either border.
	    o_2 = pl_right_ - area_.width_;
	  }
	  grid_path.M(o_1, yy).L(o_2, yy); // Horizontal grid line.
	 } // major_grid_on

	// Draw major ticks & value label, if necessary.
	double y_tick_length = tiks.major_tick_length_;
	if (tiks.ticks_loc_ < 0)
	{ // Start ticks on the plot window border left.
			o_1 = pl_left_; // o_1 = left,
			o_2 = pl_left_; //  o_2 = right.
	}
	else if (tiks.ticks_loc_ > 0)
	{ // Start ticks on the plot window border right.
	  o_1 = pl_right_;
	  o_2 = pl_right_;
	}
	else // tiks.ticks_loc_== 0
	{ // Internal style ticks on vertical Y-axis.
	  o_1 = thisAxes.position_; // Y-axis line.
	  o_2 = thisAxes.position_;
	//TODO THIS IS WRONG PROBABLY
	}
	if (tiks.ticks_side()&left_side)     o_1 -= y_tick_length; // left
	if (tiks.ticks_side()&right_side)  o_2 += y_tick_length; // right.
	tick_path.M(o_1, yy).L(o_2, yy); // Draw the major tick.
	// leaving o_1 at the left most end of any tick,  and o_2 at the rightmost end of any tick.
	// These may be on the axis line.  yy is the vertical tick position.


	if(tiks.values_side_ != 0)
	{ // Label the tick with a value, like "1.2"
	  std::stringstream labels;
	  labels.precision(tiks.value_precision_);
	  labels.flags(tiks.value_ioflags_); // set ALL IOflags.
	  labels << value; // Example: labels.str() == "20" or "0.25" or "1.2e+015"
	  if (tiks.strip_0es_)
	  { // Remove unecessary e, +, leadings 0s.
	    std::string v = strip_e0s(labels.str());
	    labels.str(v);
	  }
		if (tiks.e_to_x10_)
		{ // Remove unecessary e, +, leadings 0s.
		 std::string v = e_to_x10(labels.str());
		 labels.str(v);
		}
	  double xx = 0; // Where to start writing from, at end of left or right tick, if any.
	  // = 0 is only to avoid unitialised warning.
	  align_style alignment = center_align;
	  rotate_style rot = tiks.label_rotation_;
	  // Adjustments to provide space from end of tick before or after writing labels.
	  if (rot == horizontal)
	  {  // Shift to center value digits on tick.
	    if (tiks.values_side_ < 0)
	    { // labels to left, so start a little to left of o_1.
	      yy += valFntSiz * 0.2;
	      xx = o_1 - valFntSiz * 0.5;
	      alignment = right_align;
	    }
	    else if(tiks.values_side_ > 0)
	    { // labels to right, so start a little to right of o_2.
	     yy += valFntSiz * 0.2;
	     xx = o_2 + valFntSiz * 0.5;
	      alignment = left_align;
	    }
	  }
	  else if (rot == upsidedown)
	   {  // Just shift up to center value digits on tick.
	    if (tiks.values_side_ < 0)
	    { // labels to left, so start a little to left of o_1.
	      yy -= valFntSiz * 0.1;
	      xx = o_1 - valFntSiz * 0.5;
	      alignment = left_align;
	    }
	    else if(tiks.values_side_ > 0)
	    { // labels to right, so start a little to right of o_2.
	      yy -= valFntSiz * 0.1;
	      xx = o_2 + valFntSiz * 0.5;
	      alignment = right_align;
	    }
	  }
	  else if (rot == uphill)
	  { // Assume some 45 slope, so need about sqrt(2) less space.
	    if (tiks.values_side_ < 0)
	    { // labels to left, so start a little to left of o_1.
	      yy -= valFntSiz * 0.2;
	      xx = o_1 - valFntSiz * 0.2;
	      alignment = right_align;
	    }
	    else if(tiks.values_side_ > 0)
	    { // labels to right, so start a little to right of o_2.
	      yy += valFntSiz * 0.2;
	      xx = o_2 + valFntSiz * 0.7;
	      alignment = left_align ;
	    }
	  }
	  else if (rot == slopeup)
	  { // Assume some 30 slope, so need about sqrt(2) less space.
	    if (tiks.values_side_ < 0)
	    { // labels to left, so start a little to left of o_1.
	      yy -= valFntSiz * 0.2;
	      xx = o_1 - valFntSiz * 0.2;

	      alignment = right_align;
	    }
	    else if(tiks.values_side_ > 0)
	    { // labels to right, so start a little to right of o_2.
	      yy += valFntSiz * 0.2;
	      xx = o_2 + valFntSiz * 0.7;
	      alignment = left_align;
	    }
	  }
	  else if (rot == downhill)
	  { // Assume some 45 slope, so need about sqrt(2) less space.
	    if (tiks.values_side_ < 0)
	    { // labels to left, so start a little to left of o_1.
	      yy += valFntSiz * 0.3;
	      xx = o_1 - valFntSiz * 0.7;

	      alignment = right_align;
	    }
	    else if(tiks.values_side_ > 0)
	    { // labels to right, so start a little to right of o_2.
	      yy -= valFntSiz * 0.3;
	      xx = o_2 + valFntSiz * 0.1;
	      alignment = left_align;
	    }
	  }
	  else if (rot == slopedownhill)
	  { // Assume some 30 slope, so need about sqrt(2) less space.
	    if (tiks.values_side_ < 0)
	    { // labels to left, so start a little to left of o_1.
	      yy += valFntSiz * 0.3;
	      xx = o_1 - valFntSiz * 0.7;
	      alignment = right_align;
	    }
	    else if(tiks.values_side_ > 0)
	    { // labels to right, so start a little to right of o_2.
	      yy -= valFntSiz * 0.3;
	      xx = o_2 + valFntSiz * 0.1;
	      alignment = left_align;
	    }
	  }
	  else if (rot == steepdown)
	  { // Assume some 45 slope, so need about sqrt(2) less space.
	    if (tiks.values_side_ < 0)
	    { // labels to left, so start a little to left of o_1.
	      yy += valFntSiz * 0.3;
	      xx = o_1 - valFntSiz * 0.5;
	      alignment = right_align;
	    }
	    else if(tiks.values_side_ > 0)
	    { // labels to right, so start a little to right of o_2.
	      yy -= valFntSiz * 0.3;
	      xx = o_2 + valFntSiz * 0.1;
	      alignment = left_align;
	    }
	  }
	  else if (rot == upward)
	  { // Tick value labels straight up vertically on Y-axis.
	    yy -= valFntSiz * 0.1;
	    if (tiks.values_side_ < 0)
	    { // labels to left, so start a little to left of o_1.
	      xx = o_1 - valFntSiz * 0.7;
	      alignment = center_align;
	    }
	    else if(tiks.values_side_ > 0)
	    { // labels to right, so start a little to right of o_2.
	      xx = o_2 + valFntSiz * 1.5;
	      alignment = center_align;
	    }
	  }
	  else if (rot == steepup)
	  { // Tick value labels straight up vertically on Y-axis.
	    yy -= valFntSiz * 0.1;
	    if (tiks.values_side_ < 0)
	    { // labels to left, so start a little to left of o_1.
	      xx = o_1 - valFntSiz * 0.5;
	      alignment = center_align;
	    }
	    else if(tiks.values_side_ > 0)
	    { // labels to right, so start a little to right of o_2.
	      xx = o_2 + valFntSiz * 1.5;
	      alignment = center_align;
	    }
	  }
	  else if (rot == downward)
	  { // Tick value labels straight down vertically on Y-axis.
	    yy -= valFntSiz * 0.1;
	    if (tiks.values_side_ < 0)
	    { // labels to left, so start a little to left of o_1.
	      xx = o_1 - valFntSiz * 1.2;
	      alignment = center_align;
	    }
	    else if(tiks.values_side_ > 0)
	    { // labels to right, so start a little to right of o_2.
	      xx = o_2 + valFntSiz * 0.7;
	      alignment = center_align;
	    }
	  }
	  else // Other rotations not yet implemented.
	    return; // Without any value labels.

	  // Sanity checks on svg coordinates.
	  if (xx <= 0)  throw std::runtime_error("Y-tick X location negative!");
	  if (yy <= 0)  throw std::runtime_error("Y-tick Y location negative!");


	  if(tiks.ticks_loc_ != 0) // either on plot window or 'on axis'.
	  { // External to plot window style left or right. Always want all values including "0", if labeling external to plot window.
	      g(g_Y_TICKS_VALUES).text( xx, yy, labels.str(), tiks.value_label_style_,  alignment, rot);
	  }
	  else
	  { // ! tiks.thisTikson_pl_window_ == 0 'Internal' - value labels either side of vertical Y-axis.
	    if ((value != 0) && thisAxes.axis_line_on())
	    { // Avoid a zero ON the Y-axis if it would be cut through by any horizontal X-axis line.
	      g(g_Y_TICKS_VALUES).text( xx, yy, labels.str(), tiks.value_label_style_, alignment, rot );
	    }
	  }
	} // want value labels on tick
} // draw__major_tick


void draw_y_minor_tick(double value, path_element& tick_path, path_element& grid_path)
{ //! Draw a Y-axis minor tick and optional grid. (minor ticks do not have value labels).
	double yy(value); // Tick position and value label,
	transform_y(yy); // Convert yy to svg.
	double o_1(0.); // Start on vertical Y-axis line.
	double o_2(y_size()); // right edge of image.
	const auto& tiks = y_ticks_;

	if(tiks.minor_grid_on_)
    {
	  if(!pl_window_on_)
	  {
	    if(x_label_info_.text() != "")
	    {
	      o_1 += tiks.value_label_style_.font_size() * text_margin_;
	      o_2 -= tiks.value_label_style_.font_size() * text_margin_;
	    }
	  }
	  else // pl_window_on_ == true.
	  {
	    o_1 = pl_left_ + area_.width_;
	    o_2 = pl_right_ - area_.width_; // Ensure just *inside* window?
	  }

	    // Note comparisons are 'upside-down' - yy is increasing downwards!
	    grid_path.M(o_1, yy).L(o_2, yy); // Draw horizontal grid line.
	} // y_minor_grid

	// Draw yy minor ticks.
	if(tiks.ticks_loc_ < 0)
	{ // Put yy minor ticks on the plot window border left.
	  o_1 = pl_left_;
	  o_2 = pl_left_;
	}
	else if (tiks.ticks_loc_ > 0)
	{ // Put yy minor ticks on the plot window border left.
	  o_1 = pl_right_;
	  o_2 = pl_right_;
	}
	else // Internal style, tiks.ticks_loc_ == 0
	{
	  o_1 = thisAxes.position_; // On the Y-axis line itself.
	  o_2 = thisAxes.position_;
	}
	if(tiks.ticks_side()& left_side)		o_1 -= tiks.minor_tick_length_;
	if(tiks.ticks_side()& right_side)		o_2 += tiks.minor_tick_length_;

	if((yy <= pl_bottom_) && (yy >= pl_top_)) // Make sure that we are drawing inside of the allowed plot window.
	  tick_path.M(o_1, yy).L(o_2, yy); // Draw the horizontal tick.

} //  draw__minor_tick






void y_axis_auto_locate()
{
	auto& tiks = y_ticks_;
        // calculates the min & max, increments & ticks.

	bool scalMinMax = (tiks.max_<tiks.min_);
	if(scalMinMax)
	{
		tiks.min_ = std::numeric_limits<double>::max();
		tiks.max_ = std::numeric_limits<double>::lowest();
	 }
	for(const auto& ser:serieses_) 	if (ser.series_.size())	//autoscale_y(ser.series_);
	{
		scale_axis(YAxis(),ser.series_.begin(), ser.series_.end(), // All the container.
			&tiks.min_, &tiks.max_, &tiks.major_ticks_frac_,
			autoscale_check_limits_, // autoscale_plusminus_,
			y_include_zero_, // y_tight_,
			tiks.value_precision_, scalMinMax);
	}
	tiks.major_ticks_frac_ = std::abs(tiks.major_ticks_frac_); // in case


	//! X-axis position is determined by the range of Y min and max label values.
	if (tiks.min_ > std::numeric_limits<double>::min()) // All Y values definitely > zero.
	{ // y_min_ > 0, so X-axis will not intersect Y-axis, so put X-axis line on bottom plot window.
		if(x_ax_.location_ & at_auto)           x_ax_.location_ = at_min; // X-axis to bottom.
		if(x_ticks_.ticks_loc_ & auto_place)    x_ticks_.ticks_loc_ = left_of_plot;
		if(x_ticks_.values_side_& auto_side)    x_ticks_.values_side_ = auto_left;
		if(x_ticks_.ticks_sides_& auto_side)    x_ticks_.ticks_sides_ = auto_right;
	}
	else if(tiks.max_ < -(std::numeric_limits<double>::min)())  // all Y values definitely < zero.
	{ // // y_max_ < 0, so X-axis will not intersect Y-axis, so put X-axis line on top plot window.
		if(x_ax_.location_ & at_auto)            x_ax_.location_ = at_max;
		if(x_ticks_.ticks_loc_ & auto_place)  	x_ticks_.ticks_loc_ = right_of_plot;
		if(x_ticks_.values_side_& auto_side)    x_ticks_.values_side_ = auto_right;
		if(x_ticks_.ticks_sides_& auto_side)    x_ticks_.ticks_sides_ = auto_left;
	}
	else
	{//WARNING: on_axis does not reserve space for thick values
		if(x_ax_.location_ & at_auto)
			x_ax_.location_ = (tiks.min_ < epsT(double) && tiks.max_ > -epsT(double) )	?	at_zero  :  at_min; // Assume Y-axis will intersect X-axis (X range includes zero).
		if(x_ticks_.ticks_loc_   & auto_place)    x_ticks_.ticks_loc_ = left_of_plot;
		if(x_ticks_.values_side_& auto_side)    x_ticks_.values_side_ = auto_left;
		if(x_ticks_.ticks_sides_& auto_side)    x_ticks_.ticks_sides_ = auto_right;
	}


	double longstLabl = tiks.longest_label(y_include_zero_);// Calculate the number of chars of the longest value labels.


	// Check that labels won't collide and advise if they will - seems very difficult.
	// Change rotation to avoid collision - not practical.

	tiks.label_max_space_ = 0.; // Work out space for yy labels, depending on orientation.
	if (tiks.label_rotation_ == horizontal) // Move plot left edge right to give space for y_tick_values_precision_ digits.
	  tiks.label_max_space_ += longstLabl; // SVG units (default pixels).
	else if((tiks.label_rotation_ == upward) || (tiks.label_rotation_ == downward)) // Only need one char & 1 space width from Y-axis value label.
	 tiks.label_max_space_ += 2 * tiks.value_label_style_.font_size();
	else // Assume some slope 45, so diagonally down from tick, and takes a bit less room.
	 tiks.label_max_space_ += longstLabl * sin45;

}


#undef thisAxes
