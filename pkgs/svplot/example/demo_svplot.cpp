/*! \file demo_2d_plot.cpp
    \brief Demonstration of some 2D plot features.
    \details Uses some simple math functions to generate curves.
    The detailed output shows the plot settings for each plot.

    \author Ali Q. Raeini
*/

#include <memory>
#include <iostream>
#include <svplot.hpp>


#include <map>
#include <vector>

#include <cmath>
#include <iostream>

using namespace svg;


void plot(svplot& my_plot,	std::vector<std::vector<dbl2> >& data,
	const std::string& title = std::string(),
    const std::string& x_label = std::string(), // default is ""
    double xmin = 0., double xmax = 0.,
    const std::string& y_label = std::string(),
    double ymin = 0., double ymax = 0.,
    int x_tick_values = -1, int y_tick_values = -1, // None
    int x_rotation = horizontal, int y_rotation = upward)
{
	my_plot.title(title).legend_on(true);
	my_plot.x_label(x_label).y_label(y_label)
          .y_tick_values_side(left_side).x_tick_values_side(down_side);
			//.x_range(xmin, xmax).x_num_major_ticks(2.).x_num_minor_ticks(1)  // plus 1 major = 5 ticks per major step.
			//.y_range(ymin, ymax).y_num_major_ticks(2.).y_num_minor_ticks(1); // plus 1 major = 2 ticks per major step.
	//my_plot.x_axis().axis_line_on_ = true;
	my_plot.x_ticks_location(on_axis);
	my_plot.x_axis_location(at_zero);
	my_plot.y_axis_location(at_zero);
	//!;my_plot.legend_on(true).legend_place(outside_right) .legend_title_font_size(16); // Legend settings.
          //.legend_title("Unicode &#x3A9;&#x3A6;") // Omega Phi
           //.plot_window_on(true).plot_border_color(red)


  //{my_plot.background_color(white).plot_border_width(1);  } // image
	//.legend_background_color(white) .legend_border_color(white) .plot_background_color(svg_color(white)) .plot_border_color(svg_color(green)). title_color(red)
         
  // axis settings.

  // Very pale blue grid - like old fashioned graph paper.
  //{my_plot.x_major_grid_color(svg_color(200, 220, 255))
         //.y_major_grid_width(2)
         //.y_minor_grid_width(1)  // But nothing shows - until you make .major_grid_on(true)!
         //.x_major_grid_on(true);  }

  //my_plot.x_ticks_down_on(true); // X-axis.
  //my_plot.y_ticks_right_on(true); // Y-axis.

  /*// Options for x and/or y num_minor_ticks.
  // .y_num_minor_ticks(4)  // 0 major, 2,4,6,7,8 minor, 10 major ...
  // .y_num_minor_ticks(1) // 0 major, 5, minor, 10 major ...

  // Where the ticks (and labels if any) go, left/right, on axis, or bottom/top.
  // Default x_ticks_on_window_or_axis == -1 left or bottom, +1 right to top, 0 = on axis.
  // my_plot.x_ticks_on_window_or_axis(+1); //
  // my_plot.y_ticks_on_window_or_axis(+1);
  // my_plot.x_ticks_on_window_or_axis(-1); // right or top.
  // my_plot.y_ticks_on_window_or_axis(-1);
  // x_ticks_on_window_or_axis == 0 : on axes, if possible.

  // Which side of axis line or plot window the value labels go.
  //my_plot.x_tick_values_side(0); // NO value labels.
  //my_plot.y_tick_values_side(0); // NO value labels.

  //my_plot.x_tick_values_side(top_of_plot); // Top side value labels.
  //my_plot.x_tick_values_side(bottom_of_plot); // Bottom side value labels (default).
  //my_plot.y_tick_values_side(no_labels); // NO value labels.
  //my_plot.y_tick_values_side(right_of_plot); // Right side of axis value labels.
  //my_plot.y_tick_values_side(left_of_plot); // Left side value labels (default).

  //my_plot.x_tick_values_rotation(rotate_style(x_rotation));// Use this plot function's defaults.
  //my_plot.x_tick_values_rotation(horizontal); // svplot default.
  //my_plot.x_tick_values_rotation(upward);
  //my_plot.x_tick_values_rotation(downward);
  //my_plot.x_tick_values_rotation(uphill);
  //my_plot.x_tick_values_rotation(downhill);
  */
  //my_plot.x_ticks_on_window_or_axis(0); // ticks on axes.
  //my_plot.y_ticks_on_window_or_axis(0); // ticks on axes.

  //my_plot.y_data_ioflags(ios::dec | ios::fixed).y_tick_values_precision(1);
  //my_plot.x_data_ioflags(ios::dec | ios::scientific).x_tick_values_precision(2);
  //  my_plot.x_data_ioflags(ios::dec).x_tick_values_precision(2);

  my_plot.plot(data[0], "Sqrt(x)").fill_color(red);
  my_plot.plot(data[0], "Sqrt(x)");
  my_plot.plot(data[0], "Sqrt(x)");
  my_plot.plot(data[1], "-2 + x"+superscript("2")).fill_color(orange).marker_size(5).bezier_curve(true);
  my_plot.plot(data[2], "-1 + 2x").fill_color(yellow).line_color(blue).line_width(3).shape(squaret);
} // plot


 // Several simple math functions to demonstrate:
double f(double x){  return sqrt(x);  }
double g(double x){  return -2 + x * x;  }
double h(double x){  return 10 + 4 * x;  }
double sq(double x){  return x * x;  }
double recip(double x){  return 1. / x;  }

int main()
{

 	std::vector<std::vector<dbl2 > > datas(3,std::vector<dbl2 >(21));

	for(double i = 0; i <= 20.; i += 1.)  {
		double x=i-10;
		datas[0][i].a = x;
		datas[1][i].a = x;
		datas[2][i].a = x;
		datas[0][i].b = f(std::max(x,0.0));
		datas[1][i].b = g(x);
		datas[2][i].b = h(x);
	}

   // Demonstrate/test plots with various range of x and y, some *not* including zero.


	svgraphic my_svg(2,2); 
	my_svg.size(1000, 680);// Size/scale settings. order is important

	{
		svplot& my_plot = my_svg.subplot<svplot>(1);
		plot(my_plot, datas, "ddd", "X-axis &#x00B1;&#x3A9;", -10., +10., "Y-axis &#x221E;&#x221A;", -10., +10.); // Both X & Y include zero.
		/////////////"Title with Unicode  -&#945;   &#x3A9; &#x3A6; &#x221A; &#x221E; &#x3B6; &#x00B1; &#x2080; &#x2081; &#x2082;&#x2083;",
	}
	{
		svplot& my_plot = my_svg.subplot<svplot>(0).legend_on(true);
		my_plot.plot((piece<dbl2>(&datas[2][0]+5,5)*=0.1)+=dbl2(0.0,5.0), "Sqrt(x/10+5)")
			.fill_color(none).line_color(blue).dasharray({5});
		//my_plot.plot(piece<dbl2>(&datas[2][0]+5,5)*=0.1, "Sqrt(x)").fill_color(blank).line_color(blue);
		my_plot.x_label("x (m)").y_label("y (km)");
	}
	{
		
		svplot& my_plot = my_svg.subplot<svplot>(0,1);
		my_plot.plot(piece<dbl2>(&datas[1][0]+5,5)*=0.1, "Sqrt(x)").fill_color(red);
		my_plot.x_label("x (m)").y_label("y (km)");
	}
	{
		dbls errBar(datas[1].size(),0.0001);
		svplot& my_plot = my_svg.subplot<svplot>(1,1);
		my_plot.plot((piece<dbl2>(&datas[1][0]+5,5)*0.0001)+=dbl2(0.0,0.000001), "Sqrt(x)").fill_color(red).y_error_bars(errBar);
		my_plot.x_label("x (m)").y_label("y (km)");
		std::vector<double > aa;
		my_plot.plot(aa, aa, "Sqrt(x)").fill_color(red).y_error_bars(errBar);
	}

  std::string file = "./demo_svplot_XYPM.svg";
  std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> " << file << std::endl;
  my_svg.write(file);
  //show_2d_plot_settings(my_svg->subplot<svplot>(1));

  return 0;
} // int main()

