
// Copyright Paul A. Bristow 2006 - 2013.
// Copyright Ali Q. Raeini 2019

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)









void transform_point(double& x, double& y)
{ //!< Scale & shift both X & Y to graph Cartesian coordinates.
	x = x_scale_ * x + x_shift_;
	y = y_scale_ * y + y_shift_;
	//void adjust_limits(double& x, double& y)	adjust_limits(x, y); // In case either hits max, min, infinity or NaN.
	{ //! If value of a data point reaches limit of max, min, infinity, use the appropriate plot min or max value instead.
		// If value is NaN, use zero instead. TODO Do we want/get a different color or shape for NaNs?
		if(std::isnan(x))
		{
		  x = 0.;		  transform_x(x);
		}
		if(std::isnan(y))
		{
		  y = 0.;		  transform_y(y);
		}

		if(!std::isfinite(4.0*x))
		{
		  if(x>0.0)			x = pl_right_;
		  else			    x = pl_left_;
		}

		if(!std::isfinite(4.0*y))
		{
		  if(x<0.0)  y = pl_top_;
		  else       y = pl_bottom_;
		}
	}

}



  void transform_x(double& x)  {  x = x_scale_ * x + x_shift_;  }//!< Scale and shift X value only.
  void transform_y(double& y)  {  y = y_scale_ * y + y_shift_;  } //!< Scale and shift Y value only.





  void draw_title()
  { //! Draw title (for the whole plot). Update title_info_ with position.  Assumes align = center_align
	 if(!title_info_.data_size()) return;

	 title_info_.x(svchart::x_size() / 2.); // Center of image.
	 title_info_.y(title_info_.textstyle().font_size() * text_margin_);// Leave a linespace above.
	 svchart::g(g_TITLE).push_back(new text_element(title_info_));
  }



//Annotations:
	//#define svplot_EXTRA_SHAPES
	//#include <svg_elements_extra.hpp>



 //private:
 public: // Temporary for gil experiments.

 // svplot Implementation Member Functions.

 //! \cond DETAIL


// Legends   ************************************************************

void place_legend_update_window()
{ //! Calculate how big the legend box needs to be.
	if(legend_on_ == false)
	{ // No legend wanted, so set values to show legend positions invalid?
		legend_height_ = 0.; // At least set the size to zero.
		legend_width_ = 0.;
		return;
	}
	else // legend_on_ == true
	{	// Work out the size the legend box needs to be to hold the header, markers, lines & text.

		float font_size = legend_header_.textstyle().font_size();
		float point_size =  serieses_[0].marker_.size();
		if(serieses_[0].pointScale_.size()>1) point_size=std::max(point_size, serieses_[0].marker2_.size());
		// Use height of whichever is the biggest of point marker and font.
		double spys = std::max(font_size, point_size);
		bool draw_series_lines = false;
		double longest = string_svg_length(legend_header_.text(), legend_style_);

		// legend_height must be *at least* enough for any legend header and text_margin(s) around it (if any) plus a text_margin_ top and bottom.
		legend_height_ = spys; // At top

		legend_width_ = 2 * (legend_.margin() * legend_.width());
		// Don't plan to write on either side border, or on the 'forbidden' margins of the box.
		for(size_t i = 0; i < serieses_.size(); ++i)
		 if(serieses_[i].title_.size())
		 { // Find the longest text (if longer than header) in all the data sries.
			double siz = string_svg_length(serieses_[i].title_, legend_style_);
			if (siz > longest)  longest = siz;
			draw_series_lines |= serieses_[i].line_style_.is_on();

			if(serieses_[i].pointScale_.size()>1 && serieses_[i].seriesClr_.size()>1)  legend_height_  += 2*spys; //colour/size-bar
			legend_height_ += spys * 2; // Space for the series legend text.

		 }

		legend_width_ += longest * 0.8; // Space for longest text. Kludge factor to allow for not knowing the real length.

		legend_width_ += spys * 0.5;		// Allow width for a leading space, and trailing space before box margin.


		if (draw_series_lines) // Add for colored line marker in legend.
		  legend_width_ += spys * 3.5;
		else if(serieses_[0].marker_.shape() != no_point)   // Add for any colored data marker, cross, round... & space.
		  legend_width_ += 1. * point_size;


		// Add height depending on the number of lines of text.
		if ( (legend_header_.text() != "") )  legend_height_ += 2 * font_size; // text & space after header.


 //! Place legend box (if required).
	//outside_legend_on_ = true; // Unless proves to be inside.
	//double spys = legend_.margin(); // Might be better to use this, but needs redoing.
	spys = y_label_info_.textstyle().font_size() * 1.;  //TODO check // Around any legend box - beyond any border.
	switch (legend_place_)
	{
	 case nowhere:
	  return; // Actually places it at (0, 0), probably overwriting the plot.
	 case somewhere:  // Assume legend_top_left will place it somewhere where there is nothing else.
		break; //TODO
	 case inside_left:
		legend_left_ = pl_left_ + 2*pl_box_.margin_+0.25*legend_width_; // FIXME: what the hell
		legend_top_ = pl_top_+ 2*pl_box_.margin_; // Level with top of plot window.
		break;
	 case inside_right:
		//outside_legend_on_ = false;
		if (legend_left_ < 0.0)// Legend box position NOT been set by legend_top_left.
		{	// Default inside position is top left level with plot window.
			legend_left_ = pl_right_ - 2*pl_box_.margin_ - legend_width_; // left image edge + space.
			legend_top_ = pl_top_+ 2*pl_box_.margin_; // Level with top of plot window.
		}  // else Legend position has been specified by legend_top_left.
		break;
		 // If outside then reserve space for legend by reducing plot window.
	 case outside_right: // Default legend position is outside_right,  so that it isn't too close to the image edge or the plot window.
	  pl_right_ -= legend_width_ + spys; // Narrow plot window from right.
	  legend_left_ = pl_right_  + spys; // plot + border.
	  legend_top_ = pl_top_; // Level with top of plot window.
	  break;
	case outside_left:
	  pl_left_ += legend_width_ + spys /2 ; // Push plot window right same amount to make room,
	  legend_left_ = pl_box_.width_ + pl_box_.margin_; // left image edge + space.
	  legend_top_ = pl_top_; // Level with top of plot window.
	  break;
	case outside_top:  // Centered.
		legend_left_ = svchart::x_size()/2. - legend_width_/2; // Center.
		pl_top_ += legend_height_ + spys;
		legend_top_ = title_info_.y() + title_info_.textstyle().font_size() * text_margin_;
		legend_top_ += spys;
	  break;
	case outside_bottom: // Centered.
		legend_top_ = ( svchart::y_size() - (pl_box_.width_ + pl_box_.margin_) ) - legend_height_;
		legend_left_ = svchart::x_size()/2. - legend_width_/2; // Center.
		pl_bottom_ = legend_top_;
		pl_bottom_ -= 2*spys;
	  break;
	} // switch

	  // Check if the location requested will fit, 	  // now that we know the size of box needed.
	  if ( (legend_left_ < 0) || (legend_left_+legend_width_ > svchart::x_size())) // left outside image?
		 std::cout << "Legend top left " << legend_left_ << " is outside image size = " << svchart::x_size() << std::endl;
	  if ((legend_top_ < 0) || (legend_top_+legend_height_ > svchart::y_size())) // top outside image?
		 std::cout << "Legend top " << legend_top_ << " outside image!" << svchart::y_size() << std::endl;

	} // if legend_on_
} //  place_legend_update_window()


	//! AQR: update_internals needs improvements
  void update_internals()
  { //! The plot window is used to set a clip path:
	//! this ensures that data points and lines (and anything else)
	//! outside this window are NOT drawn.
	//! All calculation use svg units, pixels by default.

	//std::reverse( serieses_.begin(), serieses_.end() );


	// Start by assuming we can use all the svg image,  but reduce by the width of any image border.
	pl_left_ = 0 + pl_box_.width_; // Top left of image.
	pl_top_ = 0 + pl_box_.width_;
	pl_right_ = Self_.Dx_ - pl_box_.width_; // Bottom right of image.
	pl_bottom_ = Self_.Dy_ - pl_box_.width_;

	if(title_info_.text() != "")
	{ // Leave space at top for title.  TODO what if want to put title at bottom?
	  pl_top_ += title_info_.textstyle().font_size() * (text_margin_ + 0.7);
	}


	// Assume that X-axis labels are always at bottom.
	if(x_label_info_.text() != "")
	{ // Leave space at bottom for X-axis label.
	  pl_bottom_ -= x_label_info_.textstyle().font_size() * text_margin_;
	   pl_right_ -= x_label_info_.textstyle().font_size() * text_margin_*1.2; //  space for tick values AQR: TODO activate if x_ticks_values are on
	}
	// Assume that Y-axis labels are always at left.
	if(y_label_info_.text() != "")
	{ // Leave space at left for Y-axis label.
	  pl_left_ += y_label_info_.textstyle().font_size() * text_margin_;
	  pl_top_ += y_label_info_.textstyle().font_size() * text_margin_*0.6;//  space for tick values AQR: TODO activate if y_ticks_values are on
	}

	if(pl_window_on_)
	{
	  // A margin is needed to allow any plot window border rectangle to show OK.
	  // A small margin is to prevent it overlapping the image border.
	  // Also allows for axis value labels that mark the min and max
	  // that must extend beyond the plot window border,
	  // if writing is vertical need only half a font, but half x_ticks_.label_max_space_ if horizontal.
	  // x_ticks_.label_max_space_ calculated this - but is not calculated yet!
	  // So just allow a few chars.


	  //double border_margin = 0.;
	  //border_margin = std::max(pl_box_.margin_, value_space);
	  pl_left_ += plot_margin;
	  pl_right_ -= plot_margin;

	  // Might need to do similar for Y-axis if anyone complains.
	  //border_margin = std::max(pl_box_.margin_, static_cast<double>(y_tick_values_style_.font_size()/2) );
	  pl_top_ += plot_margin;
	  pl_bottom_ -= plot_margin;
	}

	place_legend_update_window(); // set legend box size and location,  according to options chosen.

	x_axis_auto_locate();
	y_axis_auto_locate();



	if (x_ticks_.values_side_ != 0)
	{ // Some tick value labels.
	  if ((x_ticks_.ticks_loc_ < 0) // ticks on bottom of plot window.
		 && (x_ticks_.values_side_ & down_side) ) // & labels on bottom too.
	  {  // Contract plot window bottom edge up to make space for X value labels on bottom.
		pl_bottom_ -= x_ticks_.label_max_space_; // Move up.
	  }
	  else if ((x_ticks_.ticks_loc_ > 0) && (x_ticks_.values_side_ & up_side) ) // & x labels to top.
	  { // Move top of plot window down to give space for x value labels.
		pl_top_ += x_ticks_.label_max_space_; // Move down.
	  }
	  // else no labels on plot window (may be on mid-plot X-axis).

	} // x_ticks_. major_value_labels_side

	// Make space for any X ticks.
	if(x_ticks_.ticks_side()&down_side) // Start bottom of plot higher to give space for any external down ticks.
	  pl_bottom_ -= std::max(x_ticks_.major_tick_length_, x_ticks_.minor_tick_length_);// Avoid macro max trap!


	if (x_ax_.axis_line_on())
	{ // Want an horizontal X-axis line, so check if range includes zero, so axes intersect,
	  // and x_ax_ is svg coordinate of Y-axis (usually y = 0).
	  // If not fiX-axis to bottom (or top) of the plot window.
	  if (x_ax_.location_ <= at_min) // All Y values definitely > zero.
		//&& !(x_ticks_.ticks_loc_ < 0) ) // & not already at bottom.
	  { // y_min_ > 0 so X-axis will not intersect Y-axis, so use plot window.
		//pl_bottom_ -= x_ticks_.label_max_space_; // Move up for the ticks value labels.
		x_ax_.position_ = pl_bottom_; // Put X-axis on bottom.
	  }
	  else if (x_ax_.location_ >= at_max)  // All Y values definitely < zero.
		//&& !(x_ticks_.ticks_loc_ > 0) ) // & not already at top.
	  { // // y_max_ < 0 so X-axis will not intersect Y-axis, so use plot window.
		 pl_top_ += x_ticks_.label_max_space_; // Move down for labels.
		 x_ax_.position_ = pl_top_; // Put X-axis on top.
	  }
	  //else // y_ax_.location_ == y_intersects_x,  Calculate below after transform is calculated.

	} // if (use_x_ax_line_)






	if (y_ticks_.values_side_ != 0)
	{ // Some major tick value labels wanted.
	  if ((y_ticks_.ticks_loc_ < 0) && (y_ticks_.values_side_ < 0) ) // On left of plot window. & labels on left.
	  {  // Contract plot window left edge to right to make space for value labels on left.
		pl_left_ += y_ticks_.label_max_space_;
	  }
	  else if ((y_ticks_.ticks_loc_ > 0) // On right of plot window.
		&& (y_ticks_.values_side_ > 0) ) // & labels to right.
	  {  // Contract plot window right to left to make space for value labels on right.
	   pl_right_ -= y_ticks_.label_max_space_;
	  }
	  else
	  { // y_ticks_.ticks_loc_ == 0
		// no value labels on plot window (may be on mid-plot Y-axis line).
		// Ignore the unusual case of Y-axis line too close to the axis label.
		// In this case the value labels may overflow the plot window
		// and collide with the axis label!
		// User must change to put value label downward, or on other side of the axis line.
		// using major_value_labels_side(int d)
		// to set tick value labels to left (<0), none (==0) or right (>0).
	  }
	} // y_ticks_. major_value_labels_side
	// Make space for any Y ticks.
	if(y_ticks_.ticks_side()&left_side)
	{ // Start left of plot to right to give space for biggest of any external left ticks.
	  pl_left_ += std::max(y_ticks_.major_tick_length_, y_ticks_.minor_tick_length_); // Avoid macro max trap!
	}
	if (y_ax_.axis_line_on())
	{ // Want a vertical Y-axis line, so check if range includes zero, so axes intersect,
	  // and y_ax_ is svg coordinate of X-axis (usually x = 0).
	  // If not fiX-axis to left (or right) of the plot window.
	  if (y_ax_.location_ == at_min) {// All X values definitely > zero.
		 //&& !(y_ticks_.ticks_loc_ < 0) // & not already at left.
	   // Y-axis will not intersect X -axis, so put Y-axis line on plot window.

//int error_replace_withMaxlblWidth;
		//pl_left_ += y_ticks_.label_max_space_; // with a space.
		y_ax_.position_ = pl_left_; // Y-axis to left,
	  }
	  else if (y_ax_.location_ == at_max) {// All X values definitely < zero.
		//&& !(y_ticks_.ticks_loc_ > 0) // & not already at right.
		//pl_right_ -= y_ticks_.label_max_space_;
		 y_ax_.position_ = pl_right_;
	  }// Y-axis to right of plot window,

	  else // x_ax_.location_ == x_intersects_y,
	  { }  //Calculate below after transform is calculated.

	} // if (use_y_ax_line_)







	if (pl_top_ >= pl_bottom_) 	{  throw std::runtime_error("Plot window top >= bottom!"); 	}
	if (pl_right_ <= pl_left_) 	{     throw std::runtime_error("Plot window right <= left!"); 	}

	// Calculate scale and shift factors for transform from Cartesian (assigned to ticks) to plot.
	// SVG image is 0, 0 at top left,y increases *downwards*
	x_scale_ = (pl_right_ - pl_left_) / (x_ticks_.max_ - x_ticks_.min_);
	x_shift_ = pl_left_ - x_ticks_.min_ * x_scale_;

	y_scale_ = (pl_top_ - pl_bottom_) / (y_ticks_.max_ - y_ticks_.min_);
	y_shift_ = pl_top_ - y_ticks_.max_ * y_scale_;



	if (x_ax_.axis_line_on() && x_ax_.location_ == at_zero)
	{  // Y Range *does* include zero, so x_ax_ not yet calculated.
		double y(0.); // Use y = 0
		transform_y(y);
		x_ax_.position_ = y; // svg Y coordinate of horizontal X-axis line.
	}
	if (y_ax_.axis_line_on() && y_ax_.location_ == at_zero)
	{ // May need to calculate axes, if will intersect.
		// X Range *does* include zero, so y_ax_ not yet calculated.
		double x(0.);
		transform_x(x);
		y_ax_.position_ = x; // SVG x coordinate of vertical Y-axis.
	}










	if (pl_window_on_)  // Draw plot window rectangle with border and/or background.
	  Self_.g(g_AREA).push_back( new rect_element(pl_left_, pl_top_, (pl_right_ - pl_left_), pl_bottom_ - pl_top_));



	updated_=true;
  } // update_internals

  //! Check if a point is within the plot window (or not too far outside).
  bool is_in_window(double x, double y)
  {
	if ((x < pl_left_ - plot_margin) || (x > pl_right_ + plot_margin)
	   ||(y < pl_top_ - plot_margin) || (y > pl_bottom_ + plot_margin))
	{  (std::cout<<" ErrorOutOfWindow ").flush();   return false;    }
	else      return true;
  } // bool is_in_window(double x, double y)




 void clear_all()
 { /*! \brief Clear all layers of the plot.
	  \details
		 When writing to multiple documents, the contents of the plot
		 may change significantly between. Rather than figuring out what
		 has and has not changed, just erase the contents of the
		 legend, title... in the document and start over.
	*/
	 { //!< Clear the legend layer of the SVG plot.
		g(g_LEGEND_BOX).clear();
		g(g_LEGEND_POINTS).clear();
		g(g_LEGEND_TEXT).clear();
	 }
		g(g_BACKGROUND).clear(); //!< Clear the whole image background layer of the SVG plot.
	 { //!< Clear the X axis layer of the SVG plot.
		g(g_X_AXIS).clear();
		g(g_X_MINOR_TICKS).clear();
		g(g_X_MAJOR_TICKS).clear();
		g(g_X_LABEL).clear();
		g(g_X_TICKS_VALUES).clear();
	 }
	 { //!< Clear the Y axis layer of the SVG plot.
		g(g_Y_AXIS).clear();
		g(g_Y_MINOR_TICKS).clear();
		g(g_Y_MAJOR_TICKS).clear();
		g(g_Y_LABEL).clear();
	 }
		g(g_TITLE).clear();//!< Clear the plot title layer of the SVG plot.
		g(g_DATA_POINTS).clear(); //!< Clear the data points layer of the SVG plot.
		g(g_AREA).clear(); //!< Clear the plot area background layer of the SVG plot.
	 { //!< Clear the  grids layer of the SVG plot.
		g(g_X_MAJOR_GRID).clear();
		g(g_X_MINOR_GRID).clear();
		g(g_Y_MAJOR_GRID).clear();
		g(g_Y_MINOR_GRID).clear();
	 }

 }



void update_image()
{ //! Draw the whole SVG image.
	clear_all();
	// svg paint rules are that later 'painting' writes over  previous painting, so the order of drawing is important.


	// Draw image background (perhaps with border and/or fill color).
	Self_.g(g_BACKGROUND).push_back( new rect_element(0, 0, Self_.x_size(),  Self_.y_size()));
	if(serieses_.empty()) {std::cout<<"Error no srieses in plot"<<std::endl; return;}


	update_internals();

	draw_title(); // Moved to ensure pl_X and Y are valid.

	// Define the clip path for the plot window.
	// We don't want to allow too much overlap of the plot window border lines, thus the minor adjustments.
 	// allow round point can lie on the axis line without being chopped in half or not show at all!!!
	Self_.clip_path(rect_element(pl_left_-1, pl_top_-1, pl_right_-pl_left_+2, pl_bottom_-pl_top_+2),  pl_window_clip_);

	Self_.g(g_DATA_POINTS).clip_id(pl_window_clip_);

	// Draw axes, labels & legend, as required.
	draw_x_axis(); // Must do X-axis first.
	draw_y_axis();
	if(legend_on_)                 	draw_legend();
	if(x_label_info_.text() != "") 	draw_x_ax_label();
	if(y_label_info_.text() != "") 	draw_y_ax_label();


	for(const auto& sri:serieses_)
	{
		if(sri.line_style_.bezier_curve_)	draw_bezier_lines(sri); // curved.
		else if(sri.line_style_.is_on()) 	draw_straight_lines(sri);// if order is important
	  //else  {   }// No line joining points.
	}
	draw_series_points();
	draw_series_bars();
	draw_series_histograms();
	draw_series_error_bars();
} // void update_image()
//! \endcond


//! Add line between sries of data points (straight rather than a Bezier curve).
//! Area fill with color if specified.
void draw_straight_lines(const svplot_series& sries)
{
	g_element& g_ptr = Self_.g(g_DATA_LINES).add_g_element(sries.line_style_);
	g_ptr.clip_id(pl_window_clip_);
	path_element& path = g_ptr.path();
	path.style().fill_color(sries.line_style_.fill_);
	bool is_fill = sries.line_style_.fill_.is_visible();
	path.style().fill_on(is_fill); // Ensure includes a fill="none" if no fill.
	//path.style().fill_color(sries.line_style_.fill_); // Duplicates so no longer needed?

	size_t inside_window = 0;  // OK data points that lie inside the plot window.

	// If required to fill the area under the plot,
	// we first have to move from the X-axis (y = 0) to the first point,
	// and again to the X-axis (y = 0) at the end after the last point.

	if (sries.series_.size() < 2)
	{ // Need at least two points for a line joining them.
	  std::cout << "Only " << sries.series_.size() << " points in sries " << sries.title_ << ", so no line drawn!" << std::endl;
	  // Also need two point  *inside window*, but that is checked later.
	}
	else
	{
	  auto j = sries.series_.cbegin();

	  double prev_x = 0; // Previous X and
	  double prev_y = 0; // Y of a data point.
	  auto prev_ux = (*j).a; // x
	  auto prev_uy = (*j).b; // y
	  // Assigned to keep compiler quiet at line 1817 "warning C4701: potentially uninitialized local variable 'prev_x' used."
	  double y0 = 0.; // y = 0, so is start point for fill area on horizontal X-axis.

	  // Try to find a first point inside the plot window. why do we care!
	  // It may not be the first point in the sries.
	  while (j != sries.series_.end())
	  {
		prev_ux = (*j).a; // x
		prev_uy = (*j).b; // y
		prev_x = prev_ux; // 1st point X-value.
		prev_y = prev_uy; // 1st point Y-value.

		transform_point(prev_x, prev_y);
		//if ((prev_x < pl_left_) || (prev_x > pl_right_) || (prev_y < pl_top_) || (prev_y > pl_bottom_))
		//if (is_in_window(prev_x, prev_y) == false)
		//{  ++outside_window;   ++j; 	} // Data point is OUTside plot window, so can't draw a line from y = 0 to this point. so try the next point to see if that is 'good' - inside the window.
		//else
		{ // Point is inside plot window, so is usable as a 1st point.
		  ++inside_window;
		  if (is_fill == true)
		  { // Move to 1st point.
			//path.style().fill_color(sries.line_style_.fill_); // Duplicates so no longer needed?
			transform_y(y0);
			path.M(prev_x, y0); // Start on X-axis
			path.L(prev_x, prev_y); // and draw line to 1st point.
			// This is to ensure fill.
		  }
		  else
		  {
			path.M(prev_x, prev_y);  // Just move to 1st good X point.
		  }

		  break;  // and continue drawing lines below in another while loop below.
		  // leaving j pointing to 2nd point.
		}
	  } // while j != series_.end()
	  if (inside_window == 0)
	  {
		std::cout << "No start point in sries " << sries.title_ << " is within plot window!" << std::endl;
		return;		// So no point trying to draw a line!
	  }//  else  {  }

	  double temp_x(0.);
	  double temp_y;
	  for(; j != sries.series_.end(); ++j)
	  {
		temp_x = (*j).a;
		temp_y = (*j).b;

		transform_point(temp_x, temp_y);
		//if (is_in_window(temp_x, temp_y) == true)
		//if ((temp_x >= pl_left_) && (temp_x <= pl_right_) && (temp_y >= pl_top_) && (temp_y <= pl_bottom_))
		{ // Data point is inside or on plot window, so draw a line to the point.
		  ++inside_window;
		  path.L(temp_x, temp_y); // Line to next point.
		  prev_x = temp_x;
		  prev_y = temp_y;
		}
		//else  ++outside_window; // Ignore any data point values outside the plot window. Not sure what this will do if area fill chosen.

	  } // for j'th point

	   //ensure(inside_window -1 + outside_window == sries.series_.size());

	  if(is_fill == true)
	  { // Area fill wanted.
		temp_y = 0; // y = 0, so on horizontal X-axis.
		transform_y(temp_y); // X-axis line SVG coordinate.
		path.L(temp_x, temp_y).z(); // Draw line to X-axis & closepath with Z.
		// to ensure area below is filled.
	  }
	}

  } // draw_straight_lines

  void draw_bezier_lines(const svplot_series& sri)
  { //! Add Bezier curve line between data points.
	//! Warning: not tested,  At present it is assumed that all data points lie within the plot window.  If this is not true, then strange and unpredictable curves will be produced!

	g_element& g_ptr = Self_.g(g_DATA_LINES).add_g_element(sri.line_style_);
	g_ptr.clip_id(pl_window_clip_);
	path_element& path = g_ptr.path();


	std::pair<double, double> n; // current point.
	std::pair<double, double> n_minus_1; // penultimate.
	std::pair<double, double> n_minus_2; // 'pen-penultimate'.
	//std::pair<double, double> fwd_vtr;
	std::pair<double, double> back_vtr;

	if( ! sri.line_style_.fill_.is_visible())
	{
	  path.style().fill_on(false); // Default path constructor is false.
	}
	else
	{ // !is_blank so DO want area fill.
	  path.style().fill_color(sri.line_style_.fill_);
	}

	if(sri.series_.size() > 2)
	{ // Need >= 3 points for a cubic curve (start point, 2 control points, and end point).
	  auto iter = sri.series_.cbegin();
	  auto un_minus_1 = *(iter++); // 1st unc X & Y data.
	  n_minus_1 = std::make_pair(un_minus_1.a, un_minus_1.b); // X and Y values.
	  //n_minus_1 = *(iter++);  // begin()
	  transform_point(n_minus_1.first, n_minus_1.second);
	  // Should check that point is inside plot window. TODO?
	  {
		auto un = *(iter++); // Middle point of trio for bezier.
		n = std::make_pair(un.a, un.b); // X and Y values.
	  }
	  transform_point(n.first, n.second);
	  // Should check that point is inside plot window. TODO?
	  path.M(n_minus_1.first, n_minus_1.second); // move m_minus_1, the 1st data point.

	  double control = 0.1;
	  // 0.2 is a scaling factor that Jake used to define the magnitude of the  vector of the current control point to be placed, basically
	  // taking advantage of the auto-drawing of Bezier curves that exists in  the SVG format, and this is his attempt to give the control point the proper length.

	  // Experiment suggests that 0.2 gives distorsions with exp curves.  0.05 is just visually OK with 50 points, but 100 are better.

	  for(; iter != sri.series_.end(); ++iter)
	  {
		n_minus_2 = n_minus_1;
		n_minus_1 = n;
		{
		  auto un = *iter; // middle
		  n = std::make_pair(un.a, un.b); // X and Y values.
		}
		transform_point(n.first, n.second);
		// Should check that point is inside plot window. TODO?

		back_vtr.first = ((n_minus_1.first - n.first) + // (x diff - x previous diff) * control
		  (n_minus_2.first - n_minus_1.first)) * control;
		back_vtr.second = ((n_minus_1.second - n.second) + // y
		  (n_minus_2.second - n_minus_1.second)) * control;

		// 8.3.6 The cubic Bezier curve commands path.S(x, y).
		// Start point, end point, & two control points.
		// Example: S378.5,519.3 381,519.3 ...
		// S end_control_point, end point
		// Start is reflection of last point's control point.
		path.S(n_minus_1.first + back_vtr.first, // x
		  n_minus_1.second + back_vtr.second, // y - end control point
		  n_minus_1.first, n_minus_1.second); // x, y - end point.
	  } // for
	  // Last point.
	  back_vtr.first = 0.;
	  back_vtr.second = (n.second - n_minus_1.second) * control;
	  path.S(n.first + back_vtr.first, // x
		n.second + back_vtr.second, // y
		n.first, n.second); // x, y
	}
	else
	{ // Only one or two points, so no curving possible!
	  draw_straight_lines(sri);
	}
  } // draw_bezier_lines


  void draw_series_bars()
  { //! Draw normal bar chart for 'good' non-limit points.
	double x(0.);
	double y(0.); // Cartesian coord y = 0.
	double x0(0.); // Y-axis line.
	double y0(.0); // X-axis line.
	transform_y(y0); // SVG coordinate of horizontal X-axis line.
	transform_x(x0); // SVG coordinate of vertical Y-axis line.
	for(unsigned int i = 0; i < serieses_.size(); ++i)
	{	const auto& sri = serieses_[i];

	  //if (sri.marker_.is_bar_style()) continue;// No bars wanted for this series.
	  g_element& g_ptr = Self_.g(g_DATA_POINTS).add_g_element(sri.marker_); // Moved up out of loop.

	  path_element& path = g_ptr.path();
	  //path.fill(sri.bar_style_.fill_ != blank);
	  path.fill_on(false);

	  double h_w = sri.marker_.size(); // For block bar chart.

	  for(auto j = sri.series_.cbegin(); j != sri.series_.end(); ++j)
	  {
		auto ux = j->a;
		x = ux;//.value();
		auto uy = j->b;
		y = uy;//.value();
		transform_point(x, y);
		if((x > pl_left_)  && (x < pl_right_) && (y > pl_top_)  && (y < pl_bottom_))
		{ // Is inside plot window, so some bar to draw.
		  switch(sri.marker_.shape_)
		  { // -2 block to Y-axis,-1 stick to Y-axis, none, +1 stick to X-axis, -2 block to X-axis.
		  case bar_y_block: // Draw a rectangle centered on the data point horizontally to Y-axis.
			 {
			   g_ptr.style().stroke_width(sri.line_style_.width_) // line_width used for rectangle line width.
				 .fill_color(sri.marker_.fill_);
			   double h_left = x;
			   double h_top = y - h_w / 2; // Start a half-width above the data point center.
			   path.M(h_left, h_top).L(h_left, h_top + h_w).L(x0, h_top + h_w).L(x0, h_top).z();
			 }
			 break;
		  case bar_y_stick:
			 path.style().stroke_width(sri.marker_.size()); // bar_width used for stick line width.
			 path.M(x, y).L(x0, y); // Draw a line from point horizontally to Y-axis.
			 break;
		  case no_point:
			 break; // Already handled above, so should not get here.
		  case bar_x_stick:
			 path.style().stroke_width(sri.marker_.size()); // bar_width used for stick line width.
			 path.M(x, y).L(x, y0); // Draw a line from point vertically to X-axis.
			 break;
		  case bar_x_block: // Draw a rectangle centered on the data point vertically to X-axis.
		   {
			 g_ptr.style().stroke_width(sri.line_style_.width_) // line_width used for rectangle line width.
			   .fill_color(sri.marker_.fill_);
			 double h_left = x - h_w / 2; // Start a half width left of the data point center.
			 double h_top = y;
			 path.M(h_left, h_top).L(h_left + h_w, h_top).L(h_left + h_w, y0).L(h_left, y0).z();
		   }
			break;
		  default:
			break;
		  } // switch
		} // for
	  } // for normal points
	}
	// Ignore all the 'bad' at_limit points.
  } //  void draw_series_bars()

  void draw_series_histograms() //! obselete: use normal plot with area_fill instead, it looks like this one takes a CDF!
  { /*!
	 Draw a histogram with variable width but contiguous bins.
	 Histograms differ from bar charts in the *area* denotes the value,
	 whereas the bar *height* denotes the value for a bar chart.
	 bin widths are provided from the X-axis data sries values.
	 The 1st data X-value provides the start of the 1st bin,
	 the 2nd data X-value provides the end of the 1st bin,
	 and the 1st Y-value the area of the 1st bin,
	 and the start of the second bin, and so on, until the
	 width of last bin is calculated from the last data point in sries,
	 that must have a zero area.  ? NaN
	 Bins can be the same (most common) or different widths.
	 Intervals must not overlap and bins must be adjacent.
	 http://en.wikipedia.org/wiki/Histogram

	 Attempts to allow a row or horizontal were abandoned because of complications
	 with the use of map which orders the x values providing the bins.
	 Using the y values for the bins implies changing the Y axes labeling and scaling too.
	*/

	g_element& g_ptr = Self_.g(g_DATA_POINTS).add_g_element(); // Moved up out of loop.
	for(unsigned int i = 0; i < serieses_.size(); ++i)
	{ const auto& sri = serieses_[i]; // for each data sries.
	  if (sri.marker_.shape_ != column_histogram)	continue;// No histogram wanted for this sries.

	  // Get the color scheme.
	  g_ptr.style().stroke_color(sri.line_style_.stroke_); // stroke color around bin blocks.
	  g_ptr.style().fill_color(sri.line_style_.fill_);
	  g_ptr.style().stroke_width(sri.line_style_.width_); // line_width used for stick line width.

	  path_element& path = g_ptr.path();
	  path.fill_on(sri.line_style_.fill_.is_visible());
	  if (path.fill_on() == true)
		path.style().fill_color(sri.line_style_.fill_);
	  else
		path.style().fill_color(blank);


	  auto last = sri.series_.cend();
	  last--; // Final pair with first the last bin end, and value zero or NaN.
	  auto u = last->b;
	  if (u != 0)//.value()
	  {
		std::cout << "Last bin end " << last->a << " should have zero value! but is "  << last->b << std::endl;
		// Or Throw? or skip this sries?
	  }
	  for(auto j = sri.series_.cbegin(); j != last; ++j)
	  { // All the 'good' 'real' data points.
		auto ux = j->a;
		double x = ux;//.value();
		auto uy =  j->b;
		double y = uy;//.value();
		auto j_next = j;
		j_next++;
		if (j != last)
		{ // Draw a column (perhaps filled) to show bin.
		  auto ux_next= j_next->a;
		  double x_next = ux_next;//.value();
		  double w = x_next - x;
		  double h = y / w;
		  double yy = h;
		  double y0(0.); // X-axis line.
		  transform_y(y0); // SVG y coordinate of horizontal X-axis line.
		  transform_x(x); // SVG x coordinate of start of bin,
		  transform_x(x_next);  // SVG x coordinate of end of bin,
		  transform_y(yy); // SVG y coordinate of height of bin.
		  //if((x > pl_left_)  && (x < pl_right_) && (y > pl_top_)  && (y < pl_bottom_))
		  //{ // Is inside plot window, so some columns to draw. TODO checks?
			path.M(x, y0).L(x, yy) // Draw a line from point vertically from X-axis.
			  .L(x_next, yy) // & horizonally to next bin end (next x value).
			  .L(x_next, y0) // back to X-axis.
			  .Z(); // So will fill.
		 } // if
	  } // for sries
	} // for normal points.
	// Ignore all the 'bad' at_limit points.
  } //  void draw_series_histograms()


  void draw_series_error_bars()
  {
	g_element& g_ptr = Self_.g(g_DATA_POINTS).add_g_element(); // Moved up out of loop.
	for(unsigned int i = 0; i < serieses_.size(); ++i)
	{ const auto& sri = serieses_[i]; // for each data sries.
	  if (sri.marker2_.shape_ != y_error_bar && sri.marker2_.shape_ != x_error_bar)	continue;// No histogram wanted for this sries.
	  int width=sri.marker2_.stroke_width();

	  // Get the color scheme.
	  g_ptr.style().stroke_color(sri.marker2_.stroke_); // stroke color around bin blocks.
	  //g_ptr.style().fill_color(sri.line_style_.fill_);
	  g_ptr.style().stroke_width(width); // line_width used for stick line width.

	  path_element& path = g_ptr.path();
	  //path.fill_on(sri.line_style_.fill_.is_visible());
	  //if (path.fill_on() == true) 	 path.style().fill_color(sri.line_style_.fill_);
	  //else
	  path.style().fill_color(blank);

	  auto last = sri.series_.cend(); // Final pair with first the last bin end, and value zero or NaN.

	  if(sri.pointScale_.size())
	  {	 auto stdvp = sri.pointScale_.begin();
		  for(auto j = sri.series_.cbegin(); j != last && stdvp!=sri.pointScale_.end(); ++j)
		  {
			double x = j->a;
			double y1 = j->b-*stdvp;
			double y2 = j->b+*stdvp;

			transform_x(x); // SVG x coordinate of start of bin,
			transform_y(y1); // SVG y coordinate of height of bin.
			transform_y(y2); // SVG y coordinate of height of bin.

			path.M(x-width*2, y1).L(x+width*2, y1).L(x, y1).L(x, y2).L(x+width*2, y2).L(x-width*2, y2); // Draw a line from point vertically.
			++stdvp;
		  }
	  }
	  if(sri.pointScale2_.size())
	  {	 auto stdvp = sri.pointScale2_.begin();
		  for(auto j = sri.series_.cbegin(); j != last && stdvp!=sri.pointScale2_.end(); ++j)
		  {
			double x = j->a;
			double y = j->b;
			double stdev=*stdvp;

			transform_x(x); // SVG x coordinate of start of bin,
			transform_y(y); // SVG y coordinate of height of bin.

			path.M(x-stdev, y).L(x+stdev, y) .Z(); // Draw a line from point vertically.
			++stdvp;
		  }
	  }

	}

  } //  draw_series_error_bars






void draw_legend()
{ //! Draw the legend border, text header (if any) and marker lines and/or shapes.
  // Assume legend box position has already been sized and positioned by function calculate_legend_box.

	// Use whichever is the biggest of point marker and font.
	float spys = std::max(legend_header_.textstyle().font_size(), serieses_[0].marker_.size());


	g(g_LEGEND_BOX)
		.style().fill_color(legend_.fill()) //
		.stroke_color(legend_.stroke())
		.stroke_width(legend_.width());	  //.stroke_on(legend_.border_on())


	  g(g_LEGEND_BOX).rect(legend_left_, legend_top_, legend_width_, legend_height_);	  // Draw border box round legend.


	float lg_y = legend_top_ + text_margin_*spys;
	if (legend_header_.text().size())
	{ // Draw the legend header text for example: "My Plot Legend".
	  legend_header_.x(legend_left_ + legend_width_/2.); // / 2. to center in legend box.
	  // Might be better to use center_align here because will fail if legend contains symbols in Unicode.
	  legend_header_.y(lg_y);
	  g(g_LEGEND_TEXT).push_back(new  text_element(legend_header_));
	  lg_y += 2*spys; // Might be 1.5? - useful if many sries makes the box too tall.
	}


	bool draw_series_lines=false;
	for(unsigned int i = 0; i < serieses_.size(); ++i)
	 if(serieses_[i].title_.size())
	 { // Show point marker, perhaps line, & text info for each of the data sries.

		auto& sri = serieses_[i];
		draw_series_lines |= sri.line_style_.is_on() ;


		float lg_x = legend_left_ + spys; // + space before point marker and/or line & text.
		// Use both stroke & fill colors from the points' style.		// Applies to both shape AND line.		//g_inrptr->style().stroke_color(sri.marker_.stroke_);		//g_inrptr->style().fill_color(sri.marker_.fill_);		//g_inrptr->style().stroke_width(sri.line_style_.width_);


		if(sri.line_style_.fill_.is_visible()) // draw a rectangle around legend
			g(g_LEGEND_POINTS).add_g_element(sri.line_style_).rect( lg_x-spys/2,  lg_y-spys/2, spys, spys);
		else
		if (draw_series_lines)
		{
			g(g_LEGEND_POINTS).add_g_element(sri.line_style_).
				push_back(new line_element( lg_x-spys, lg_y, lg_x+spys, lg_y, sri.line_style_)); // Draw horizontal lines with appropriate color. // Line sample is one char long.
		}

	  if(sri.marker_.is_on())
	  {
		  if(sri.seriesClr_.size()) ;
			//draw_pl_point(lg_x, lg_y, *g_inrptr, sri.marker_, sri.seriesClr_[j], sri.pointScale_[j], sri.marker2_);
		  else
			draw_pl_point(lg_x, lg_y, g(g_LEGEND_POINTS).add_g_element(sri.marker_), sri.marker_);

	  }

	  if(draw_series_lines)   lg_x += 1.5 * spys;// Total is short line & a space.
	  else if(sri.marker_.is_on()  || sri.line_style_.fill_.is_visible())    lg_x += 1.0 * spys;

	  // Legend text for each Data Series added to the plot.
	  g(g_LEGEND_TEXT).text(
		 lg_x, lg_y+0.25*spys, // allow space for the marker.
		 sri.title_, // Text for this data sries.
		 legend_header_.textstyle(),  left_align);
	  lg_y += 2 * spys;


	  if(sri.seriesClr_.size()>1)
		draw_colorbar({float(legend_left_)+0.5f*spys,lg_y}, {0.8f*spys,6*spys},4,2,sri);

	 } // for if
} // void draw_legend()


void draw_colorbar(float2 xy0, float2 Dxy , const unsigned int nmajor , const unsigned int nminor,  const svplot_series& sri)
{ //! Draw the legend border, text header (if any) and marker lines and/or shapes.

	const point_style& sty = sri.marker_;  const point_style& sty2 = sri.marker2_;
	int font_size  = legend_header_.textstyle().font_size();
	int point_size =  sty.size();

	float spys = std::max(font_size, point_size);

	float dy=Dxy[1]/(nminor*nmajor);


	g_element& g_leg_ptr = g(g_LEGEND_POINTS);
	double y1 = xy0.b + Dxy[1];

	for(unsigned int i = 0; i <= nminor*nmajor; ++i)
	{ // Show point marker, perhaps line, & text info for each of the data sries.

		g_element& g_ref = g_leg_ptr.add_g_element();
		float clrf = float(i)/(nminor*nmajor);
		g_ref.style()
			.fill_color(sri.clrinterf_(clrf, sty.fill_color(), sty2.fill_color()))
			.stroke_color(sri.clrinterf_(clrf, sty.stroke_color(), sty2.stroke_color()));


		g_ref.push_back(new rect_element(  xy0.a - spys/2 ,  y1- dy/2,  spys,  dy));

		if(i%nminor==0)
		{

			std::stringstream label;
			  label.precision(y_ticks_.value_precision_+2);
			  label.flags(y_ticks_.value_ioflags_); // set ALL IOflags.
			  label << sri.seriesClr_.scalefrom01(clrf);
			  if (y_ticks_.strip_0es_)
			  { // Remove unecessary e, +, leadings 0s.
				std::string v = strip_e0s(label.str());
				label.str(v);
			  }
			  if (y_ticks_.e_to_x10_)
			  { // Remove unecessary e, +, leadings 0s.
				 std::string v = e_to_x10(label.str());
				 label.str(v);
			  }

			g(g_LEGEND_TEXT).push_back(new text_element(
			 xy0.a + spys,	 y1+ 0.5*dy,
			 label.str() , // Text for this data sries.
			 legend_header_.textstyle(),
			 left_align));
		 }
		y1 -= dy;
	} // for
} // void draw_colourbar()





// point markers ********************************************************

void draw_pl_point(double x, double y, g_element& g_ptr, const point_style& sty	)
{ //aq209: Warning this does not apply style, apply it in the group if using this, to save space
	//TODO: move symbol (text characters) to g too. !?
	//! Draw a plot data point marker shape
	//!  whose size and stroke and fill colors are specified in point_style,

	float size = sty.size();
	float half_size = size / 2.f;
	float quarter = size / 4.f; // use baseline once it is supported instead of adhoc mess here

	switch(sty.shape_) // from enum point_shape no_point, round, square, point, egg
	{
	case no_point:    break;
	case circlet: 	  g_ptr.circle(x, y, half_size);  	  break;
	case point: 	  g_ptr.circle(x, y, 1); 	  break;// Fixed size round.
	case squaret:	  g_ptr.rect(x - half_size, y -half_size, size*0.9, size*0.9);	  break;
	case egg:    	  g_ptr.ellipse(x, y, half_size, size * 2.);	  break; // Tall thin egg!

	case unc_ellipse:	  break;

	 // Offset from center is not an issue with vertical or horizontal ticks.

	case vertical_tick:   g_ptr.line(x, y, x , y - size); 	  break; // Especially neat for 1-D points. // tick up from axis.
	case vertical_line:   g_ptr.line(x, y + size, x , y - size);	  break; // line up & down from axis.

	case horizontal_tick: 	  g_ptr.line(x, y, x + size, y ); 	  break;	  // horizontal_tick is pretty useless for 1-D because the horizontal line is on the X-axis. // tick right from axis.
	case horizontal_line: 	  g_ptr.line(x, y - size, x + size, y );	  break; // line left & right from axis.	  // horizontal_line is pretty useless for 1-D because the horizontal line is on the X-axis.

	case symbol:	  g_ptr.text(x, y+quarter, sty.symbols(), sty.txtStyle(), center_align, horizontal);   break;// symbol(s), size and center.
	case diamond: 	  g_ptr.text(x, y+quarter, "&#x2666;", sty.symbols_style_, center_align, horizontal);  break;	  // size / 4. puts bottom tip on the X-axis,	  // size / 2. put center above the X-axis	  // x, y, put on the X-axis - probably what is needed for 2-D plots.	  // diamond, spades, clubs & hearts fill with expected fill_color.
	case asterisk: 	  g_ptr.text(x, y+quarter, "&#x2217;", sty.symbols_style_, center_align, horizontal);  break;	  // asterisk is black filled.	  // size /3 puts the bottom tip on the X-axis.
	case lozenge:	  g_ptr.text(x, y+quarter, "&#x25CA;", sty.symbols_style_, center_align, horizontal);  break;	  // size / 3 to get tip of lozenge just on the X-axis.	  // lozenge seems not to fill?
	case club:   	  g_ptr.text(x, y+quarter, "&#x2663;", sty.symbols_style_, center_align, horizontal);  break;	  // x, y, puts club just on the X-axis

	case spade: 	  g_ptr.text(x, y+quarter, "&#x2660;", sty.symbols_style_, center_align, horizontal);  break;
	case heart: 	  g_ptr.text(x, y+quarter, "&#x2665;", sty.symbols_style_, center_align, horizontal);  break;
	case cone:  	  g_ptr.text(x, y+quarter, "&#x25BC;", sty.symbols_style_, center_align, horizontal);  break;
	//case cone:  	  g_ptr.triangle(x - half_size, y - size, x + half_size, y - size, x, y, sty.fill_.is_visible() );	  break; // Pointing down triangle.  // Last point puts the bottom tip of the triangle on the X-axis (may not be wanted for 2-D).

	case trianglet:	  g_ptr.text(x, y+quarter, "&#x25B2;", sty.symbols_style_, center_align, horizontal);  break; // Pointing up triangle.
		  // Also could use &#x25BC for pointing down triangle, and &#x25B4 for small up-pointing triangle and &#x25BE for small down triangle.

	case star:  	g_ptr.text(x, y+quarter, "&#x2605;", sty.symbols_style_, center_align, horizontal);   break;
	case cross:   	g_ptr.text(x, y+quarter, "&#x2715;", sty.symbols_style_, center_align, horizontal);   break;	  // Cross is pretty useless for 1-D because the horizontal line is on the X-axis.

	default:   break;   // TODO Other point_shapes do nothing yet.
	}
} // void draw_pl_point


void draw_pl_point(double x, double y, float clrf, float sizef, g_element& g_ptr, const svplot_series& sri )	//, const point_style& sty 2 // X and Y values (in SVG coordinates).//const point_style& sty,
{

	switch (sri.marker2_.shape_)
	{
		case y_error_bar:
		case x_error_bar:
			draw_pl_point(x,y,g_ptr,sri.marker_);
			break; //handelled separately
		default:
		{
			point_style styl = sri.marker_;
			styl.size(sri.marker_.size() *(1.0f-sizef)+sizef* sri.marker2_.size());
			g_element& g_ref = g_ptr.add_g_element();
			g_ref.style()
				.fill_color(sri.clrinterf_(clrf,  sri.marker_.fill_color(),   sri.marker2_.fill_color()))
				.stroke_color(sri.clrinterf_(clrf,sri.marker_.stroke_color(), sri.marker2_.stroke_color()));
			draw_pl_point(x,y,g_ref,styl);
		} break;
	}
}



void draw_series_points()
{ //! Draw normal 'good' non-limit points, and then any 'at limits' points.
	double x(0.);
	double y(0.);
	for(unsigned int i = 0; i < serieses_.size(); ++i)
	{ const auto& sri = serieses_[i];

		if(!sri.marker_.is_on()) continue;
		//if(sri.marker_.fill_.is_blank()) cout<<"Warning blank fill colour get rendered black sometimes"
	  g_element& g_ptr = Self_.g(g_DATA_POINTS).add_g_element(sri.marker_);

		for (size_t j = 0; j<sri.series_.size(); ++j) {
			x = sri.series_[j].a;
			y = sri.series_[j].b;
			transform_point(x, y); // Note x and y are now SVG coordinates.
			if((x > pl_left_-0.5) && (x < pl_right_+0.5) && (y > pl_top_-0.5) && (y < pl_bottom_+0.5))
			{ // Data point is inside plot window, so draw a point.
				if(sri.pointScale_.size() || sri.seriesClr_.size())
					draw_pl_point(x, y, sri.seriesClr_[j], sri.pointScale_[j], g_ptr, sri);
				else
					draw_pl_point(x, y, g_ptr, sri.marker_);
				// TODO might refactor so that only pass ux, and uy.
				if (x_values_on_ && y_values_on_)
				{ // Show the two values of the X & Y data as a pair.
						//g_element& g_ptr_vx = Self_.g(g_X_POINT_VALUES).add_g_element();
						//g_element& g_ptr_vy = Self_.g(g_Y_POINT_VALUES).add_g_element();
						//draw_point_values(x, y, g_ptr_vx, g_ptr_vy, x_values_style_, y_values_style_,0,0);//, ux, uy
				}
				else if (x_values_on_)
				{ // Show the value of the X data point too.
				// void draw_point_value(double x, double y, g_element& g_ptr, value_style& val_style, point_style& point_style, double value)
						//g_element& g_ptr_vx = Self_.g(g_X_POINT_VALUES).add_g_element();
				//draw_point_value(x, y, g_ptr_vx, x_values_style_, sri.marker_);//, ux
				}
				else if (y_values_on_) // show the value of the Y data point too.
				{	  //g_element& g_ptr_vy = Self_.g(g_Y_POINT_VALUES).add_g_element();
				//draw_point_value(x, y, g_ptr_vy, y_values_style_,sri.marker_);//, uy
				}

			}
	  }

	}

	//! Draw the abnormal 'at_limit' points (if any)  -- deactivated


} //  void draw_series_points()
