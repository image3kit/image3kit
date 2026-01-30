
// Copyright Paul A. Bristow 2006 - 2013.
// Copyright Ali Q. Raeini 2019

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#define thisAxes   x_ax_


void draw_x_axis()
{ //! Draw horizontal X-axis line & plot window line to hold, and ticks and grids.
	const auto& tiks = x_ticks_;
	  //const double eps = std::numeric_limits<double>::epsilon();
	if(thisAxes.axis_line_on()) // Want a horizontal X-axis line drawn.
	{
		double x_1 = pl_left_;
		double x_2 = pl_right_;
		if (x_ax_.location_ == at_zero) // Draw the horizontal X-axis line the full width of the plot window, // perhaps including an addition in lieu of a major tick.
		{
		  if (y_ticks_.ticks_side()&left_side)
		  {
			 if (y_ticks_.ticks_loc_ < 0) // at left  // Extend the horizontal line left in lieu of longest tick.
				x_1 -= std::max(y_ticks_.minor_tick_length_, y_ticks_.major_tick_length_);
		  }
		  else if (y_ticks_.ticks_side()&right_side)//side of axis, not plot
		  {
			 if (y_ticks_.ticks_loc_ > 0) // at right // Extend the horizontal X-axis line right in lieu of longest tick.
				x_2 += std::max(y_ticks_.minor_tick_length_, y_ticks_.major_tick_length_);
		  }
		  double yy = thisAxes.position_; // yy = 0, (provided yy range includes zero).
		  g(g_X_AXIS).line(x_1, yy, x_2, yy);
		  if (tiks.ticks_loc_ < 0) // bottom // Draw a vertical line holding the ticks on the top of plot window.
			 g(g_X_AXIS).line(x_1, pl_bottom_, x_2, pl_bottom_);
		  else if (tiks.ticks_loc_ > 0)  // top,  Draw a vertical line holding the ticks on the bottom of plot window.
			 g(g_X_AXIS).line(x_1, pl_top_, x_2, pl_top_);
		}
		else if (x_ax_.location_ & at_max)
			g(g_X_AXIS).line(x_1, pl_top_, x_2, pl_top_);
		else if (x_ax_.location_ & at_min)
			g(g_X_AXIS).line(x_1, pl_bottom_, x_2, pl_bottom_);
		else {	} // warn that things have gone wrong?
	} // axis_line_on()



	// Access the paths for the ticks & grids, ready for additions.
	path_element& minor_tick_path = g(g_X_MINOR_TICKS).path();
	path_element& major_tick_path = g(g_X_MAJOR_TICKS).path();
	path_element& minor_grid_path = g(g_X_MINOR_GRID).path();
	path_element& major_grid_path = g(g_X_MAJOR_GRID).path();

	// dMinor is the interval between minor ticks.
	double dMajor = (tiks.max_-tiks.min_)*tiks.major_ticks_frac_;
	double dMinor = dMajor / (double(tiks.num_minor_ticks_ + 1.) );
	// TODO Problem here with using floating point?
	// Was i < x_max_; but didn't show the tick and value at x_max_ so now i <= x_max_;
	// But may still fail if a ls or few bits out? Seems to fail for xx = 100, for example.

	// Draw the ticks on the positive side.
	//int nMTicks=0;
	double xx = x_include_zero_ ? 0.0 : std::max(0.0,std::floor(tiks.min_/dMajor)*dMajor);
	double j  = x_include_zero_ ? 0.0 : std::max(0.0,std::floor(tiks.min_/dMinor)*dMinor) + dMinor;
	for(; xx <= tiks.max_ + 2*std::numeric_limits<double>::epsilon(); xx += dMajor)
	{
	  for(; j < (xx + dMajor) * (1. - 2*std::numeric_limits<double>::epsilon()); j+=dMinor)
	  {
	    if (j != 0. || ! y_ax_.axis_line_on()) // Avoid a major tick at xx == 0 where there *is* a horizontal X-axis line.
	      draw_x_minor_tick(j, minor_tick_path, minor_grid_path);
	  }
	  // Draw major tick.
	  if ((xx != 0. || ! y_ax_.axis_line_on()) // axis line requested.
	    || (tiks.ticks_loc_ != 0)) // ticks & labels on plot window.
	  { // Avoid a major tick at xx == 0 where there *is* a horizontal X-axis line.
	    draw_x_major_tick(xx, major_tick_path, major_grid_path);
	      //++nMTicks;
	  } //else  // std::cout << "Missed xx " << xx << std::endl; // Only miss 0s
		j = xx + dMinor;
	}

	// Draw the ticks on the negative side.
	for(xx = x_include_zero_ ? 0.0 : std::min(0.0,tiks.max_);
	     xx >= tiks.min_; xx -= dMajor)
	{
	  for(j = xx; j > xx - dMajor; j-= dMajor / (tiks.num_minor_ticks_ + 1))
	  { // Draw minor ticks.
	    if ( (j != 0. || ! thisAxes.axis_line_on())
	      || (tiks.ticks_loc_ != 0) // ticks & labels on plot window.
	      )
	    { // Avoid a major tick at xx == 0 where there *is* a horizontal X-axis line.
	      draw_x_minor_tick(j, minor_tick_path, minor_grid_path);
	    }
	  }
	  if ((xx != 0. || ! y_ax_.axis_line_on())
	    || (tiks.ticks_loc_ != 0) ) // ticks & labels on plot window.
	  { // Avoid a major tick at xx == 0 where there *is* a horizontal X-axis line.
	    draw_x_major_tick(xx, major_tick_path, major_grid_path);
	      //++nMTicks;
	  }
	}
	//if(nMTicks<2)
		//std::cout<<"Error only "<<nMTicks<<" major ticks"<<std::endl;
} // void draw_axis






void draw_x_ax_label()
{ //! Draw the X-axis label text (for example,  length (km) ),

	const auto& tiks = x_ticks_;
	 const  float  valFontSize = tiks.value_label_style_.font_size();
	 const  float  lblFontSize = x_label_info_.textstyle().font_size();
	 // Glyphs for western characters are aligned with the left bottom of capital letter, so need to allow  for any descenders.
	 double yy = pl_bottom_;
	 if (tiks.ticks_loc_ & bottom_of_plot) // -1 means ticks on bottom of plot window.
	 { // Ticks value labels below plot window.
		if (tiks.values_side_ < 0) // bottom
		{ // Shift down to allow for any tick value labels.
		  if ((tiks.label_rotation_ == downward) || (tiks.label_rotation_ == upward))
		  { // downward tick value label direction 90 up or down.
			 yy += tiks.label_max_space_;
			 if (tiks.ticks_side()&down_side)
			 {  // Move down for any downward ticks.	// and a small space.
				yy += 1.1 * std::max(tiks.major_tick_length_, tiks.minor_tick_length_);
				yy += 0.7 * (lblFontSize + valFontSize); // best compromise?
			 }
		  }
		  else if ((tiks.label_rotation_ == steepdown) || (tiks.label_rotation_ == steepup))
		  { // downward tick value label direction 60 up or down.
			 yy += tiks.label_max_space_;
			 if (tiks.ticks_side()&down_side)
			 {  // Move down for any downward ticks, and a small space..
				yy += 1.1 * std::max(tiks.major_tick_length_, tiks.minor_tick_length_);
				yy += 0.5 * (lblFontSize + valFontSize); // best compromise?
			 }
		  }
		  else if ((tiks.label_rotation_ == uphill)  || (tiks.label_rotation_ == downhill))
		  { // sloping 45 degrees up or down.
			 yy += tiks.label_max_space_ * sin45; // Move down from end of tick.
			 if (tiks.ticks_side()&down_side)
			 {  // Move down for any downward ticks, and a small space.
				yy += 1.1 * std::max(tiks.major_tick_length_, tiks.minor_tick_length_);
				yy += 0.7 * (lblFontSize + valFontSize); // best compromise?
			 }
		  }
		  else if ((tiks.label_rotation_ == slopeup)  || (tiks.label_rotation_ == slopedownhill))
		  { // sloping 30 degrees.
			 yy += tiks.label_max_space_ * sin45; // Move down from end of tick.
			 if (tiks.ticks_side()&down_side)
			 {  // Move down for any downward ticks.
				yy += 1.1 * std::max(tiks.major_tick_length_, tiks.minor_tick_length_);
				// and a small space.
				yy += 0.5 * (lblFontSize + valFontSize); // best compromise?
			 }
		  }
		  else if (tiks.label_rotation_ == horizontal)
		  { // horizontal X ticks value labels (default).
			 if (tiks.values_side_ < 0)  //  Move down to allow space for font size of tick value labels below X-axis.
				yy += valFontSize ;

			 yy += lblFontSize * 1.3; // Allow for the X-axis label font and space.
			 // See also 1.3 factor drawing ticks.
		  }
		  else
			 std::cout << " Rotation of X label rotation" << tiks.label_rotation_ << "not yet implemented!" << std::endl;
		}
		else if (tiks.values_side_ > 0)
		{ // Tick labels above, only ticks below, so just move down for height of label font.
		  yy += lblFontSize * 1.3; // Allow for the X-axis label font and space.
		}
		else
		{ // tiks.values_side_ == 0  // So no change for labels.
		  //yy += lblFontSize * 1.3; // Allow for the X-axis label font and space.
		}

		if (tiks.ticks_side()&down_side)// Shift down for biggest of any ticks, and bit of space.
		  yy += 1.1 * std::max(tiks.minor_tick_length_, tiks.major_tick_length_);
	 }
	 else if (tiks.ticks_loc_ > 0)
	 {  // = +1 means ticks are on top of plot window.
		 // Shift down from plot window bottom to suit X-axis label.
		 yy += lblFontSize * 1.7;
	 }
	 else if (tiks.ticks_loc_ == 0)
	 { // Ticks are ON the X-axis line, so X label is just below the plot bottom.
		// No space needed for ticks.
		 // Character starts at bottom of capital letter, so allow for descenders.
		 //yy -= pl_box_.margin_;
		 yy += lblFontSize * 1.3;
	 }

	 g(g_X_LABEL).push_back(new text_element(
		( // xx position relative to the xx-axis which is middle of plot window.
		pl_right_ + pl_left_) / 2,  // xx coordinate - middle.
		yy, // Down from plot window.
		x_label_info_.text(),
		x_label_info_.textstyle(),
		center_align, horizontal)
		);
  } // void draw_thisAxeslabel()



void draw_x_major_tick(double value, path_element& tick_path, path_element& grid_path)
{ //! Draw major ticks - and grid too if wanted.
	const auto& tiks = x_ticks_;
	const  float  valFontSize = tiks.value_label_style_.font_size();
	double xxx(value); //
	transform_x(xxx); // xx value in svg.
	if((xxx < pl_left_ - 0.01) || (xxx > pl_right_ + 0.01)) return;

	double o_1(0.); // upper end of tick.
	double o_2(svchart::x_size()); // o_2 = lower end of tick.
	if(tiks.major_grid_on_)
	{ // Draw major grid vertical line.
		if(!pl_window_on_)
		{ // Allow a modest margin around text of title and X-axis labels, if in use.
		  if(title_info_.text() != "")
			 o_1 += title_info_.textstyle().font_size() * text_margin_;
		  if(tiks.values_side_ != 0) // Value may be shown either side the major tick.
			o_2 -= x_label_info_.textstyle().font_size() * text_margin_;
		}
		else
		{ // pl_window_on_ == true
		  o_1 = pl_top_; // Bottom of plot window.
		  o_2 = pl_bottom_; // Top of plot window.
		}
		//if((o_2 <= pl_bottom_) && (o_1 >= pl_top_) && (xx >= pl_left_) && (xx <= pl_right_)) // Make sure that we are drawing inside the allowed window.
		  grid_path.M(xxx, o_1).L(xxx, o_2); // Vertical grid line.
      } // major_grid_on

	 // Draw major tick (perhaps as well as grid - ticks might be wider than grid).
		double x_tick_length = tiks.major_tick_length_;
		if(tiks.ticks_loc_ < 0)
		{ // Put the ticks on the plot window border (was external).
			o_1 = pl_bottom_; // on the window line.
			o_2 = pl_bottom_;      // o_1 = upper, o_2 = lower.
		}
		else if(tiks.ticks_loc_ > 0)
		{ // Put the ticks on the plot window border (was external).
		  o_1 = pl_top_; // on the window line.
		  o_2 = pl_top_; // o_1 = upper, o_2 = lower.
		}
		else
		{ // Draw tick from the central X axis line (Internal_style).
		  o_1 = thisAxes.position_; // X-axis line.
		  o_2 = thisAxes.position_;
		}
		//TODO THIS IS WRONG PROBABLY
		if(tiks.ticks_side()&up_side)     o_1 -=  x_tick_length; // up
		if (tiks.ticks_side()&down_side)  o_2 += x_tick_length; // down.

		tick_path.M(xxx, o_1).L(xxx, o_2);
		// Leaving current position at the bottom end of the tick.
		// o_1 and yy-down are the ends of the tick.
		// These may be on the axis line, or the plot window.

		if(tiks.values_side_ != 0)
		{ // Label the tick with a value, like "1.2"
		  std::stringstream labels;
		  labels.precision(tiks.value_precision_);
		  labels.flags(tiks.value_ioflags_);
		  labels << value; // for tick "4", "1.2" or "3.4e+000"...
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
		  double yy = 0; // Where to start writing from, at end of bottom or top tick, if any.
		  // = 0 is only to avoid unitialised warning.
		  align_style alignment = center_align;
		  // Adjustments to provide space from end of tick before or after writing label.
		  if (tiks.label_rotation_ == upward) // vertical writing up.
		  {  // Shift to center value digits and minus sign on tick.
			 xxx += valFontSize * 0.2;
			 if (tiks.values_side_ < 0)
			 { // labels to bottom, so start a little below o_2.
				yy = o_2 + valFontSize * 0.6;
				alignment = right_align;
			 }
			 else if(tiks.values_side_ > 0)
			 { // labels to top, so start a little above o_1.
				yy = o_1 - valFontSize * 0.5;
				alignment = left_align;
			 }
		  }
		  else if (tiks.label_rotation_ == downward)
		  {  // Should handle other directions too.
			 xxx -= valFontSize * 0.3;
			 if (tiks.values_side_ < 0)
			 { // labels to bottom, so start a little below o_2.
				yy = o_2 + valFontSize * 0.5;
				alignment = left_align;
			 }
			 else if(tiks.values_side_ > 0)
			 { // labels to top, so start a little above o_1.
				yy = o_1 - valFontSize * 0.5;
				alignment = right_align;
			 }
		  }
		  else if (tiks.label_rotation_ == steepup)
		  {  // Should handle other directions too.
			 xxx -= valFontSize * 0.3;
			 if (tiks.values_side_ < 0)
			 { // labels upward, so start a little below o_2.
				yy = o_2 + valFontSize * 0.5;
				alignment = left_align;
			 }
			 else if(tiks.values_side_ > 0)
			 { // labels to top, so start a little above o_1.
				yy = o_1 - valFontSize * 0.5;
				alignment = right_align ;
			 }
		 }
//
		  else if (tiks.label_rotation_ == uphill)
		  { // Assume some 45 slope, so need about sqrt(2) less space.
			 xxx += valFontSize * 0.5;
			 if (tiks.values_side_ < 0)
			 { // labels to bottom, so start a little to bottom of o_1.
				yy = o_2 + valFontSize * sin45;
				// Seems to need a bit more space for top than bottom if rotated.
				alignment = right_align;
			 }
			 else if(tiks.values_side_ > 0)
			 { // labels to top, so start a little to top of o_2.
				yy = o_1 - valFontSize * 0.3;
				alignment = left_align;
			 }
		  }
		  else if (tiks.label_rotation_ == slopeup)
		  { // Assume for 30 degree slope, need about sqrt(2) less space.
			 xxx += valFontSize * 0.5;
			 if (tiks.values_side_ < 0)
			 { // labels to bottom, so start a little to bottom of o_1.
				yy = o_2 + valFontSize * sin45;
				// Seems to need a bit more space for top than bottom if rotated.
				alignment = right_align;
			 }
			 else if(tiks.values_side_ > 0)
			 { // labels to top, so start a little to top of o_2.
				yy = o_1 - valFontSize * 0.2;
				alignment = left_align;
			 }
		  }
		  else if (tiks.label_rotation_ == downhill)
		  { // Assume some 45 slope, so need about sqrt(2) less space.
			 xxx -= valFontSize * 0.3;
			 if (tiks.values_side_ < 0)
			 { // labels to bottom, so start a little to bottom of o_2.
				yy = o_2 + valFontSize * 0.7;  // Seems to need a bit more space for top than bottom if rotated.
				alignment = left_align;
			 }
			 else if(tiks.values_side_ > 0)
			 { // labels to top, so start a little to top of o_1.
			  yy = o_1 - valFontSize * 0.3;
				alignment = right_align;
			 }
		  }
		  else if (tiks.label_rotation_ == slopedownhill)
		  { // Assume some 30 slope, so need about sqrt(2) less space.
			 xxx -= valFontSize * 0.3;
			 if (tiks.values_side_ < 0)
			 { // labels to bottom, so start a little to bottom of o_2.
				yy = o_2 + valFontSize * 0.7;  // Seems to need a bit more space for top than bottom if rotated.
				alignment = left_align;
			 }
			 else if(tiks.values_side_ > 0)
			 { // labels to top, so start a little to top of o_1.
			  yy = o_1 - valFontSize * 0.3;
				alignment = right_align;
			 }
		  }
		  else if (tiks.label_rotation_ == steepdown)
		  {  // Should handle other directions too.
			 xxx -= valFontSize * 0.3;
			 if (tiks.values_side_ < 0)
			 { // labels to bottom, so start a little below o_2.
				yy = o_2 + valFontSize * 0.5;
				alignment = left_align;
			 }
			 else if(tiks.values_side_ > 0)
			 { // labels to top, so start a little above o_1.
				yy = o_1 - valFontSize * 0.5;
				alignment = right_align;
			 }
		  }
		  else if (tiks.label_rotation_ == horizontal)
		  { // Tick value label on X-axis is normal default horizontal.
			 if (tiks.values_side_ < 0)
			 { // labels to bottom of tick, so start a little below bottom of o_2.
				yy = o_2 + valFontSize * 1.3; // 1.3 allows 1/3 font space.
				alignment = center_align; // center on the tick.
			 }
			 else if(tiks.values_side_ > 0)
			 { // labels to top, so start a little to top of o_1.
			  yy = o_1 - valFontSize * 0.7;
				alignment = center_align;
			 }
		  }
		  else // Other rotations not yet implemented.
            return; // Without any value label.

		  if (xxx <= 0) // Sanity checks on svg coordinates.
			 throw std::runtime_error("X-tick X location negative!");
		  if (yy <= 0)
			 throw std::runtime_error("X-tick Y location negative!");

		  // Draw the X ticks value labels, "1", "2" "3" ...
		  if(tiks.ticks_loc_ != 0)
		  { // External to plot window style bottom or top.
			 // Always want all values including "0", if labeling external to plot window.
			 // tiks.ticks_loc_ == true != 0
			 g(g_X_TICKS_VALUES).text(
				xxx, yy, labels.str(),
				tiks.value_label_style_, // font, size etc
				alignment, tiks.label_rotation_);
		  }
		  else
		  {
			 //if(((xxx < pl_left_ + 0.01) || (xxx > pl_right_ - 0.01)) || (int(xxx*100.0+0.5)==0  && thisAxes.axis_line_on()))
			 {

				g(g_X_TICKS_VALUES).add_g_element(svg_style(blank,svg_color(255,255,255,200),1)).rect(
					xxx-valFontSize/2., yy-valFontSize*1.0, valFontSize, valFontSize*1.2);
			 }
			 //else
			 { // Avoid a "0" below the X-axis if it would be cut through by any internal vertical Y-axis line.
				g(g_X_TICKS_VALUES).text(
				  xxx, yy, labels.str(), // "1.23"
				  tiks.value_label_style_, // font, size etc
				  alignment, tiks.label_rotation_);
			 }
		  } // on plot window or 'on axis'.
	 }
} // draw_x_major_tick


void draw_x_minor_tick(double value, path_element& tick_path, path_element& grid_path)
{ //!< Draw X-axis minor ticks, and optional grid. (Value is NOT (yet) shown beside the minor tick).
	double xx(value); // Tick position and tick value label,
	transform_x(xx); // Convert to svg.
	if((xx < pl_left_ - 0.01) || (xx > pl_right_ + 0.01))  return;// allow = too? add epsilon?
	double o_1(0.); // Start on the horizontal X-axis line.
	double o_2(svchart::y_size()); // Image top.
	const auto& tiks = x_ticks_;

	if(tiks.minor_grid_on_) // Draw the minor grid, if wanted.
	{
		if(!pl_window_on_)
		{ // Use whole image.  Make space for title and X-axis labels.
		  if(title_info_.text() != "") // Allow text_margin_ * font_size around text (pixels).
			 o_1 += title_info_.textstyle().font_size() * text_margin_;
		  if(x_label_info_.text() != "")
			 o_2 -= x_label_info_.textstyle().font_size() * text_margin_;
		}
		else // pl_window_on_ == true.
		{
		  o_1 = pl_top_ + area_.width_; // Top.
		  o_2 = pl_bottom_ - area_.width_; // Bottom. Ensure *inside* window.
		}
		// Make sure that we are drawing inside the allowed window.
        grid_path.M(xx, o_1).L(xx, o_2); // Draw grid line.
	} // x_minor_grid

	 // Draw xx minor ticks.
	 if (tiks.ticks_loc_ < 0)
	 { // Put minor ticks on the plot window border bottom.
		o_1 = pl_bottom_; // on the window line.
		o_2 = pl_bottom_; // o_1 = upper, o_2 = lower end of tick.
	 }
	 else if (tiks.ticks_loc_ > 0)
	 { // Put minor ticks on the plot window border top.
		o_1 = pl_top_; // on the window line.
		o_2 = pl_top_; // o_1 = upper, o_2 = lower end of tick.
	 }
	 else // tiks.ticks_loc_ == 0
	 { // Internal style, draw tick up and/or down from the X-axis line.
		o_1 = thisAxes.position_; // ON X-axis horizontal line.
		o_2 = thisAxes.position_;
	 }
	if(tiks.ticks_side()&up_side)		o_1 -= tiks.minor_tick_length_; // up
	if(tiks.ticks_side()&down_side)		o_2 += tiks.minor_tick_length_; // down.

	 tick_path.M(xx, o_1).L(xx, o_2);	// No value labels on minor ticks, at present.

} // void draw_x_minor_tick



void x_axis_auto_locate()
{
	auto& tiks = x_ticks_;
        // calculates the min & max, increments & ticks.
	//if (tiks.autoscale() && serieses_.size())
	//{ // Use calculated autoscale values.


	bool scalMinMax = (tiks.max_<tiks.min_);
	if(scalMinMax)
	{
		tiks.min_ = std::numeric_limits<double>::max();
		tiks.max_ = std::numeric_limits<double>::lowest();
	}
	for(const auto& ser:serieses_)	 if (ser.series_.size())//		autoscale_x(ser.series_);
	{ //! Data series (all values) to use to calculate autoscaled X-axis values.
		scale_axis(XAxis(),ser.series_.begin(), ser.series_.end(), // All the container.
			&tiks.min_, &tiks.max_, &tiks.major_ticks_frac_,
			autoscale_check_limits_, autoscale_plusminus_,
			x_include_zero_, x_tight_, tiks.value_precision_, scalMinMax);
	}

	tiks.major_ticks_frac_ = std::abs(tiks.major_ticks_frac_); // in case

	// Check if the axes will intersect.
	// Y-axis position is determined by the range of X values.
	if (tiks.min_ > std::numeric_limits<double>::min()) // All Y values definitely > zero.
	{ // y_min_ > 0, so X-axis will not intersect Y-axis, so put X-axis line on bottom plot window.
		if(y_ax_.location_ & at_auto)          y_ax_.location_ = at_min; // X-axis to bottom.
		if(y_ticks_.ticks_loc_ & auto_place)     y_ticks_.ticks_loc_ = left_of_plot;
		if(y_ticks_.values_side_& auto_side)    y_ticks_.values_side_ = auto_left;
		if(y_ticks_.ticks_sides_& auto_side)    y_ticks_.ticks_sides_ = auto_right;
	}
	else if(tiks.max_ < -(std::numeric_limits<double>::min)())  // all Y values definitely < zero.
	{ // // y_max_ < 0, so X-axis will not intersect Y-axis, so put X-axis line on top plot window.
		if(y_ax_.location_ & at_auto)            y_ax_.location_ = at_max;
		if(y_ticks_.ticks_loc_ & auto_place)  	y_ticks_.ticks_loc_ = right_of_plot;
		if(y_ticks_.values_side_& auto_side)    y_ticks_.values_side_ = auto_right;
		if(y_ticks_.ticks_sides_& auto_side)    y_ticks_.ticks_sides_ = auto_left;
	}
	else
	{
		if(y_ax_.location_ & at_auto)
			y_ax_.location_ = (tiks.min_ < epsT(double) && tiks.max_ > -epsT(double) )
			?	at_zero  :  at_min; // Assume Y-axis will intersect X-axis (X range includes zero).
		if(y_ticks_.ticks_loc_   & auto_place)    y_ticks_.ticks_loc_ = left_of_plot;
		if(y_ticks_.values_side_& auto_side)    y_ticks_.values_side_ = auto_left;
		if(y_ticks_.ticks_sides_& auto_side)    y_ticks_.ticks_sides_ = auto_right;
	}

	double longstLabl = tiks.longest_label(x_include_zero_);

	// Check that labels won't collide and advise if they will - seems very difficult.
	// Change rotation to avoid collision - not practical.


	tiks.label_max_space_ = 0.; // Work out the longest ticks values label for X-Axis.
	if (tiks.label_rotation_ == horizontal)
	{ // Only 1 char height & small space needed if labels are horizontal.
		tiks.label_max_space_ += 1.5 * tiks.value_label_style_.font_size(); // 2 SVG chars
	}
	else if ((tiks.label_rotation_ == upward) || (tiks.label_rotation_ == downward)) // ! X_axis ticks labels vertical so will need enough for all the characters in the label.
		tiks.label_max_space_ += longstLabl; // in SVG units pixels.

	else // Assume label is sloping, say 45, so * sin(45) = 0.707.
		tiks.label_max_space_ += longstLabl * sin45; // SVG 'chars'.


 }


#undef thisAxes
