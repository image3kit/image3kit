/*!
  \file svg_style.hpp
  \brief Styles for SVG specifying font, sizes, shape, color etc for text, values, lines, axes etc.
  \details SVG style information is fill, stroke, width, line & bezier curve.
   This module provides struct point_style & struct line_style
   and class svg_style holding the styles.
   See http://www.w3.org/TR/SVG11/styling.html
  \date Mar 2009, (Modified in 2019 by Ali Q Raeini)
  \author Jacob Voytko and Paul A. Bristow and Ali Q raeini
*/

// Copyright Jacob Voytko 2007
// Copyright Paul A. Bristow 2008, 2009, 2013

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "svg_color.hpp"
#include "svg_elements.hpp"
#define FOLD

#ifndef BOOST_SVG_SVG_STYLE_HPP
#define BOOST_SVG_SVG_STYLE_HPP

#ifdef _MSC_VER
#  pragma warning (disable : 4512) // assignment operator could not be generated.
#endif

#include <iostream>
#include <sstream>
#include <string>

//namespace boost {
namespace svg {


// Forward declarations of classes in svg_style.hpp
class svg_style; // Holds the basic stroke, fill colors and width, and their switches.
class text_style; // Text and tspan element's font family, size ...
class value_style; // Data series point value information, text, color, orientation.
class point_style; // Shape, color, of data point markers.
class line_style; // Style of line joining data series values.
class axis_line_style; // Style of the x and/or y axes lines. But NOT the ticks and value labels.
class tick_style; // Style of the x and y axes ticks, grids and tick value labels.
class box_element; //  Box colors, size and border switches.

// Ugly hack to remove unwanted sign and leading zero(s) in exponent.
const std::string strip_e0s(std::string s);
// Estimate length of string when appears as svg units.
int string_svg_length(const std::string& s, const text_style& style, bool etoX10 = false);



class value_style
{ //! \class svg::value_style
  //! \brief Data series point value label information, text, color, orientation
  //! name ID string, order in sequence, time and date.
  //! \details For example, to output: 5.123 +- 0.01 (19).
  //! Uncertainty and degrees of freedom estimate.
  //! Prefix, separator and suffix allow X and Y values to be together on one line, for example\n
  //! [1.23+- 0.01 (3), 4.56 +-0.2 (10)]\n
  //! Used in draw_pl_point_values (note plural - not used in singular draw_pl_point_value)
  //! where X value_style is used to provide the prefix and separator, and Y value_style to provide the suffix.
  //! Prefix, separator and suffix are ignored when X or Y are shown separately using draw_pl_point_value.
  //! "4.5+- 0.01 (3) Second #2, 2012-Mar-13 13:01:00"


public:
  rotate_style value_label_rotation_; //!< Direction point value labels written.
  int value_precision_; //!< Decimal digits of precision of value.
  std::ios_base::fmtflags value_ioflags_; //!< Control of scientific, fixed, hex etc.
  bool strip_0es_; //!< If true, then unnecessary zeros and + sign will be stripped to reduce length.
  text_style val_style_; //!< Font etc used for data point value marking.
  // svg_style
  svg_color stroke_; //!< Stroke color for value.
  svg_color fill_; //!< Fill color for value.
  bool plusminus_on_; //!< If an uncertainty estimate is to be appended (as + or - value).  \details See http://en.wikipedia.org/wiki/Plus-minus_sign
  svg_color plusminus_color_; //!< Color for uncertainty, for example: 0.02 in "1.23 +-0.02".
  bool addlimits_on_; //!< If an confidence interval is to be added, for example <4.5, 4.8>.
  svg_color addlimits_color_; //!< Color for confidence interval.
  bool df_on_; //!< If a degrees of freedom estimate is to be appended.
  svg_color df_color_; //!< Color for degrees for freedom, for example: 99 in "1.23 +-0.02 (99)".
  bool id_on_;  //!< If an id or name string to be appended. default == false,
  svg_color id_color_; //!< Color for id or name string".
  bool datetime_on_; //!< If an time and/or date string to be appended. default == false,
  svg_color datetime_color_; //!< Color for time and date string".
  bool order_on_;  //!< If an order in sequence number # to be appended. default == false,
  svg_color order_color_; //!< Color for sequence number #".
  std::string prefix_; //!< Prefix to data point value, default none, but typically "[".
  std::string separator_; //!< Separator between x and y values, if both on same line (none if only X or only Y, or Y below X).
  std::string suffix_; //!< Suffix to data point value, default none, but typically "]".

  // Constructors declarations.

 //!<Default Constructor Data point value label style (provides default color and font).
 value_style()
    :
    value_label_rotation_(horizontal), //!< Label orientation, default horizontal.
    value_precision_(4), //!< Precision, reduced from default of 6 which is usually too long.
    value_ioflags_(std::ios::dec), //!< Any std::ios::ioflags, for example, hex, fixed, scientific. eg. std::ios::scientific | std::ios::dec) & ~std::ios::fixed
    strip_0es_(true), //!< If true, then unnecessary zeros will be stripped to reduce length.
    val_style_(no_style),  //!< All defaults, black etc.
    stroke_(black), //!< == black.
    fill_(svg_color(0, 0, 0)), //!< no fill.
    plusminus_on_(false), //!< If uncertainty estimate to be appended.
    plusminus_color_(black), //!< Default color for uncertainty of value.
    addlimits_on_(false), //!< If uncertainty estimate to be appended.
    addlimits_color_(black), //!< Default color for uncertainty of value.
    df_on_(false), //!< If a degrees of freedom estimate to be appended.
    df_color_(black), //!< Default color for degrees of freedom is black.
    id_on_(false), //!< If an id or name string to be appended.
    id_color_(black), //!< Default color for an id or name string is black.
    datetime_on_(false), //!< If a date and date to be appended.
    datetime_color_(black), //!< Default color for date and date is black.
    order_on_(false), //!< If a order #  to be appended.
    order_color_(black), //!< Default color for order # is black.
    prefix_(""),
    separator_(","),
    suffix_("")
    {  } //! Default constructor initialises all private data.


    //!< Constructor Data point value label style (provides default color and font).
    value_style(rotate_style r, //!< Label orientation, default horizontal.
      int p, //!< Reduced from default of 6 which is usually too long.
      std::ios_base::fmtflags f = std::ios::dec, //std::ios_base::fmtflags(), //!< Any std::ios::ioflags, for example, hex, fixed, scientific.
      bool s = true, //!< If true, then unnecessary zeros will be stripped to reduce length.
      text_style ts = no_style, //!< All defaults, black etc.
      const svg_color& scol = black, //!< == black.
      const svg_color& fcol = black,  //!< no fill.
      bool pm = false, //!< If uncertainty estimate to be appended.
      const svg_color& plusminus_color = black, //!< Default color for uncertainty of value.
      bool lim = false, //!< If confidence limits to be appended.
      const svg_color& addlimits_color = black, //!< Default color for confidence limits.
      bool df = false,  //!< If a degrees of freedom estimate to be appended.
      const svg_color& df_color = black,//!< Default color for uncertainty of value.
      bool id = false,  //!< If a degrees of freedom estimate to be appended.
      const svg_color& id_color = black,//!< Default color for uncertainty of value.
      bool dt = false,  //!< If a degrees of freedom estimate to be appended.
      const svg_color& dt_color = black,//!< Default color for uncertainty of value.
      bool ordno = false,  //!< If a degrees of freedom estimate to be appended.
      const svg_color& order_color = black,//!< Default color for uncertainty of value.
      // Separators [,] provide, for example: [1.23+-0.01 (3), 4.56 +-0.2 (10)]
      std::string pre = "", //!< Prefix, for example: "[",
      std::string sep  = "", //!< separator, for example: ,\&\#x00A0;", // If put ", " the trailing space seems to be ignored, so add Unicode explicit space.
      std::string suf  = "") //!< suffix, for example: "]")
    :
    value_label_rotation_(r), value_precision_(p), value_ioflags_(f), strip_0es_(s),
    val_style_(ts), stroke_(scol), fill_(fcol),
    plusminus_on_(pm), plusminus_color_(plusminus_color),
    addlimits_on_(lim), addlimits_color_(addlimits_color),
    df_on_(df), df_color_(df_color),
    id_on_(id), id_color_(id_color),
    datetime_on_(dt), datetime_color_(dt_color),
    order_on_(ordno), order_color_(order_color),
    prefix_(pre), separator_(sep), suffix_(suf)
    {  } //! Constructor setting parameters with some defaults.



}; // class value_style



enum point_shape
{ //! \enum point_shape used for marking a data point.
  no_point = 0, //!< No marker for data point.
  circlet, /*!< Circle. Name was changed to round to avoid clash with function named circle,
  but was then found to clash with C++ Standard numeric function round.
  Full qualification `point_shape::round` requires C++11 support to compile, so then changed to circlet.
  */
  squaret, //!< Square.
  star, //!< Star (using polygon).
  diamond, //!< Diamond card shape.
  cross, //!< cross
  cone, //!< Cone pointing up - 'rightwayup'.
  trianglet, //!< Triangle pointing down 'upsidedown'.
  asterisk, //!< Asterix as * symbol
  vertical_line,  //!< Vertical line up & down from axis.
  horizontal_line, //!< Horizontal line left & right from axis.
  vertical_tick, //!< Vertical tick up from axis.
  horizontal_tick, //!< Horizontal line right from axis.
  unc_ellipse, //!< Ellipse sized using uncertainty estimate of x and y, typically about twice standard deviation or 95% confidence interval.
  point, //!< Small solid point.
  egg, //!< Ellipse.
  lozenge, //!< Lozenge or square with corners pointing up and down..
  heart, //!< Heart playing card shape.
  club, //!< Club playing card shape.
  spade, //!< Spade playing card shape.
  symbol /*!< Unicode symbol including letters, digits, greek & 'squiggles'.
  \verbatim
    Default letter "X".\n
    Other examples: "&#x3A9;"= greek omega, "&#x2721;" = Star of David hexagram
    &#2720 Maltese cross & other dingbats. \n
    See also http://en.wikipedia.org/wiki/List_of_Unicode_characters#Basic_Latin geometric shapes
    that may be a better way to make these symbols: &#25A0 black square ...to &#25FF
    But unclear how many browsers implement these properly.
  \endverbatim
  */

  ,x_error_bar //!< horizental error bar, use as marker shape 2 only
  ,y_error_bar //!< vertical error bar, use as marker shape 2 only
  ,column_histogram,//!< Stick or column line (stroke width) vertically to/from X-axis.
  bar_y_block, //!< Rectangular (optionally filled) block style horizontal to Y-axis,
  bar_y_stick, //!< Bar or row line (stroke width) horizontal to Y-axis.
  bar_x_stick, //!< Stick or column line (stroke width) vertical to X-axis.
  bar_x_block  //!< Rectangular (optionally filled) block style vertical to X-axis,

}; // enum point_shape

class point_style : public svg_style
{ /*! \class svg::point_style
    \brief Shape, color, of data point markers.
    \details (optional value & uncertainty) not implemented yet.
  */
  //friend std::ostream& operator<< (std::ostream&, point_style);

public:

  point_shape shape_; //!< shape: round, square, point...
  std::string symbols_; //!< Unicode symbol(s) (letters, digits, squiggles etc).\n
  //! Caution: not all Unicode symbols are output by all browsers!\n
  //! see http://en.wikipedia.org/wiki/Hexagram, symbols("&#x2721;")
  //! Positioning of symbols (especially > 1 symbols) may be imprecise.

  text_style symbols_style_; //!< font, size, decoration of symbols.
  bool show_x_data_; //!< If true, show the X value like "1.2" near the point. (If both true, then show both X and Y as a pair like "1.2, 3.4".)
  bool show_y_data_; //!< If true, show the Y value like "3.4" near the point. (If both true, then show both X and Y as a pair like "1.2, 3.4".)


// Constructor.

  point_style( //!< Constructor with defaults for data members. .
    const svg_color& stroke = black,  //!< Color of circumference of shape.
    const svg_color& fill = blank, //!< Fill color of the centre of the shape.
    int size = 5, //!< Diameter of circle, height of square, font_size  ...
    point_shape shape = circlet, //!< shape: circlet, square, point...
    const std::string& symbols = "X") //!< Unicode symbol(s) (letters, digits, squiggles etc).
  : svg_style(stroke,fill),
    shape_(shape), symbols_(symbols),
    show_x_data_(false), show_y_data_(false)
  { // Best to have a fixed-width font for symbols? But there are always problems centering a symbol at the right point.
    symbols_style_.font_family("Lucida Sans Unicode").font_size(size);
    //.baseline("central"); need to manually shift the symbols for now, TODO revert back once this supported
  }

// Member Function Definitions.

	point_style& size(float i)  {    symbols_style_.font_size(i);     return *this;  } //! Set size of shape or symbol used to mark data value plot point(s).// Font size, in case using a symbol as marker.Diameter of circle, height of square, font_size


	float size() const  {    return symbols_style_.font_size();  } //! \return size of shape or symbol used to mark data value plot point(s).

	point_style& fill_color(const svg_color& f)  { fill_ = f;   return *this;	 } //! Set fill color of shape or symbol used to mark data value plot point(s).

	const svg_color& fill_color() const  {    return fill_;	 } //! \return  fill color of shape or symbol used to mark data value plot point(s).


	point_style& stroke_color(const svg_color& f)  {    stroke_ = f;    return *this;   } //! Set stroke color of shape or symbol used to mark data value plot point(s).


	svg_color stroke_color() const    {     return stroke_;	 }//! \return  stroke color of shape or symbol used to mark data value plot point(s).

	point_style& shape(point_shape s)  {     shape_ = s;    return *this;   } //! Set shape used to mark data value plot point(s).

	point_shape shape() const          {    return shape_;	 } //! \return  shape used to mark data value plot point(s).


	point_style& symbols(const std::string s)  {    symbols_ = s;    return *this;   } //! Override default symbol "X" - only effective if .shape(symbol) used.

	bool is_on() const                { return  shape_ && (stroke_on() || fill_on());  } //! \return True if line(s) will join data points.

	const std::string& symbols() const  {    return symbols_;	 } //! \return plot point marking symbol (only effective if .shape(symbol) used).


	point_style& txtStyle(text_style ts)	 { 	symbols_style_ = ts;    return *this;   }

	text_style& txtStyle() const  {     return const_cast<text_style&>(symbols_style_);  }//! \return text_style& To allow control of symbol font, size, decoration etc.


	bool is_bar_style() const	{ 	return bar_y_block<=shape_ && shape_<=bar_x_block; 	}

}; // struct point_style



class line_style : public svg_style
{ //! \class svg::line_style Style of line joining data series values.
  // TODO dotted and dashed line style would be useful for monochrome plots.
public:
  bool bezier_curve_; //!< If true, data points will be joined by bezier curved line(s), otherwise by straight lines (default if line colour is set).
  std::string dasharray_; //!< svg stroke-dasharray style

  //! Constructor to set plot line style, but providing defaults for all member data.
  line_style(const svg_color& col = black, const svg_color& fill_col = blank, double width = 2, bool bezier = false)
    : svg_style(col,fill_col,width),  bezier_curve_(bezier)    {  } // Provides defaults for all private data.


// Member Functions.

  line_style& width(double w)  {  width_ = w;    return *this; } //! Set width of line(s) joining data points. // return *this to Make chainable

  double width()  {    return width_;	 } //! \return  width of line(s) joining data points.

  line_style& color(const svg_color& f)  {    stroke_ = f;    return *this;	 } //! Set color of line(s) joining data points.
  line_style& dasharray(const std::string& str)  {  dasharray_ = str;   return *this;	 } //! Set dasharray_str of line(s) joining data points.


  line_style& area_fill(const svg_color& f)  {   fill_ = f;    return *this;	 }//! Set if area under line joining data points is to be color filled.



  bool is_on() const  {  return stroke_on() || fill_on() || bezier_curve_;  } //! \return True if line(s) will join data points.

  line_style& line_on(bool is)  {  stroke_on(is);    return *this;	 } //! Set true if line(s) are to join data points.

  bool bezier_curve() const  {   return bezier_curve_;	 } //! \return true if bezier curved line(s) are to join data points.

  line_style& bezier_curve(bool is)  {    bezier_curve_ = is;    return *this;	 } //! Set true if bezier curved line(s) are to join data points.



}; // class line_style



/*! Number of standard deviations used for text_plusminus text display.\n
 Nominal factor of 2 (strictly 1.96) corresponds to 95% confidence limit. 	*/
static const double text_plusminus = 2.;

static const double sin45 = 0.707; //!< Used to calculate 'length' if axis value labels are sloping.
static const double reducer = 0.9; //!< To make uncertainty and degrees of freedom estimates a bit smaller font to help distinguish from value.

 enum legend_places
 { //! \enum legend_places Placing of legend box, if requested by legend_on == true.
	nowhere = 0, //!< Placing of legend box not requested or not calculated yet.
	inside_left = -1,  //!< Default place for inside is top left of plot window, (exact location controlled by legend_top_left()).
	inside_right = -2,  //!< Default place for inside is top left of plot window, (exact location controlled by legend_top_left()).
	outside_left = 1, //!< Outside on the left of the graph.
	outside_right = 2, //!< Outside right (Default).
	outside_top = 3, //!< Outside at top.
	outside_bottom = 4, //!< Outside at bottom.
	somewhere = 5 //!< legend_top_left(x, y)
 };


 enum axis_intersect
 { //! \enum (x and) y_ax_intersect  If and how the Y axes intersects X axis.
	at_auto = -3, //!< default value to be rest automatically if not set by user
	at_min  = -1, //!< Y-axis free to left of end of X-axis (case of all X definitely < 0).
	at_zero = 0, //!< y_intersects_x when X values include zero, so intersects the X axis.
	at_max  = +1 //!< Y-axis free to left of end of X-axis (case of all X definitely > 0).
  };


enum dim
{ //! \enum dim dimension of plot. (Used so that an axis knows what type it is, or none = N).
  N = 0, X = 1, Y = 2
};

class axis_line_style
{ //! \class svg::axis_line_style
  //! \brief Style of the X or Y-axes lines.
  //! \details (But NOT the ticks and value labels because different styles for X and Y-axes are possible).
public:
  //dim dim_; //!< None, X or Y.
  //double min_; //!< minimum X value (Cartesian units).
  //double max_; //!< maximum Y value (Cartesian units).
  svg_color color_; //!< Axis line (stroke) color.
  double axis_width_; //!< Axis line width.
  double position_; //!< Depending on value of dim, either X-axis (y = 0) transformed into SVG Y coordinates or Y-axis (x = 0) transformed into SVG X coordinates (-1 if not calculated yet).

  axis_intersect  location_; //!< Intersection with other axis, or not.



  bool axis_line_on() {return color_.o_;}; //!< Draw an X horizontal or Y vertical axis line.

  // Default constructor.
  axis_line_style( //!< Default constructor.  Sets all member data items with defaults for all.
    const svg_color col = black, //!< Axis line color.
    double width = 1 //, //!< Axis line width.
    // bool label_on = true, //!< Label axis with text - example: "length".
    // bool axis_lines_on = true  //!< Draw an X horizontal or a Y vertical axis line.
    )
   :   color_(col), axis_width_(width),//, min_(min), max_(max)  dim_(d),
    //axis_position_(axis_position),
    position_(-1.0), // -1 means not calculated yet.   //int axis_position = 0, //!< How the axes intersect with values as below:\n
	location_(at_auto)
									 //! enum x_ax_intersect {bottom = -1, x_intersects_y = 0, top = +1};
									 //! enum y_ax_intersect {left = -1, y_intersects_x = 0, right = +1};
									 //! If axes look like an L, then is bottom left.
									 //! If a T then y intersects and X is at bottom.
  {  }



  // Set and get member functions.


  axis_line_style& color(const svg_color& color)
  {  color_ = color;    return *this;	 }//! Set color of an axis line.

  svg_color color()
  { //! \return  color of an axis line.
    return color_;	 }

  axis_line_style& width(double w)
  { //! Set width of an axis line.
    axis_width_ = w;    return *this;	 }

  double width()
  { //! \return  width of an axis line.
    return axis_width_;	 }

  //axis_line_style& position(int pos)
  //{ //! How the axes intersect.
    //axis_position_ = pos;
    //return *this;
  //}

  //double position()
  //{ //! \return How the axes intersect.
  //  return axis_position_;
  //}

  //bool axis_line_on() const
  //{ //! If returns true, then either an X or a Y axis line to be drawn.
    //return axis_line_on_;
  //}


}; // class axis_line_style

//! The place for ticks value labels on the axis.
enum tick_place : int
{//! TODO make everything positive
  auto_place = +4,  left_of_plot = (-1)^127,  on_axis = 0,  right_of_plot = +2,  bottom_of_plot = (-1)^255,  top_of_plot = +2
};
enum tick_side : int
{//! TODO make everything positive
   left_side = (-1)^127,  down_side = (-1)^127,  no_tick = 0,  right_side = +2,  up_side = +2,  both_sides = left_side|right_side,
    auto_side = 4 ,
    auto_outside = left_side | down_side | 4,
    auto_inside = right_side | up_side | 4,
    auto_right = right_side | 4,
    auto_left = left_side | 4,
};

class tick_style
{ /*! \class svg::tick_style
   \brief Style of the X and Y axes ticks, grids and their tick value labels.
   \details
   But NOT the X and Y axes lines.
   These can be either on the axis lines or on the plot window edge(s),
   (because different styles for x and y are possible).
  */
  friend class plot2d;

public:
    //dim dim_; //!< X, Y, or None.
    double min_; //!< Minimum x value (Cartesian units).
    double max_; //!< Maximum x value (Cartesian units).
    //bool up_ticks_on_; //!< Draw ticks up from horizontal X-axis line.
    //bool down_ticks_on_; //!< Draw ticks down from horizontal X-axis line.
    //bool left_ticks_on_; //!< Draw ticks left from vertical Y-axis line.
    //bool right_ticks_on_; //!< Draw ticks right from vertical Y-axis line.
    tick_place ticks_loc_; //!< Value labels & ticks on a plot window border (rather than on X or Y-axis).
    // Simplest to have all of these although only one pair (up or down) or (left or right) is used. Unused are always false.
    tick_side values_side_; //!< Which side of axis for label values for major ticks. < 0 means to left (for Y) or down (for X) (default), 0 (false) means no ticks value labels (just ticks), > 0 means to right (for Y) or top(for X). values_side_|=2 means auto
    tick_side ticks_sides_; //!< Which side of axis for label values for major ticks. < 0 means to left (for Y) or down (for X) (default), 0 (false) means no ticks value labels (just ticks), > 0 means to right (for Y) or top(for X). values_side_|=2 means auto
    rotate_style label_rotation_; //!< Direction axis value labels written.
    bool major_grid_on_;  //!< Draw X grid at major ticks.
    bool minor_grid_on_; //!< Draw X grid at minor ticks.
    unsigned int num_minor_ticks_; //!< number of minor ticks, eg 4 gives major 0, minor 1,2,3,4, major 5 (All units in svg units, default pixels).
    double  major_ticks_frac_;  //!< Number of major ticks, can also be a decimal number (= (max-min)/interval). No set function because x_num_minor_ticks_ used to determine this instead, but one could calculate x_minor_interval_.
    svg_color major_tick_color_; //!< Color (stroke) of tick lines.
    double major_tick_width_; //!< Width of major tick lines.
    double major_tick_length_;//!< Length of major tick lines.
    svg_color minor_tick_color_; //!< Color (stroke) of tick lines.
    double minor_tick_width_; //!< Width of minor tick lines.
    double minor_tick_length_; //!< Length of minor tick lines.
    svg_color major_grid_color_; //!< Color of major grid lines.
    double major_grid_width_; //!< Width of major grid lines.
    svg_color minor_grid_color_; //!< color of minor grid lines.
    double minor_grid_width_; //!< Wdith of minor grid lines.

    svg_color values_color_; //!< Color of tick values labels. (just fill_color for now (stroke makes characters fuzzy.)
    int value_precision_; //!< Precision for tick value labels, usually 3 will suffice.
    std::ios_base::fmtflags value_ioflags_;  //!< IO formatting flags for the axis default std::ios::dec.
    bool strip_0es_; //!< If redundant zero, + and e are to be stripped, for example "+1.000e3" to "1e3".
    bool e_to_x10_; //!< If redundant zero, + and e are to be stripped, for example "+1.000e3" to "1e3".
    double label_max_space_;  //!< Space (SVG units, pixels) needed for value label adjusted for rotation.
    //! For Y-axis -1 = left, 0 = false = on X-axis, +1 = right. Default -1 to left of plot window.
    //! For X-axis -1 = bottom, 0 = false = on Y-axis, +1 = top. Default -1 below bottom of plot window.
    //! 0 = false puts the ticks and their labels on the X or Y axis line which may be in the middle of the plot.
    //! For 1D the default overrides the constructor default of -1 below, to tick and value label the X-axis.
    //! For 2D the default is left at -1, to use bottom and left of plot window to tick and value label X and Y-axis.

    text_style value_label_style_; //!< text style (font, size...) for value labels.

     //! Constructor setting several parameters, but providing default values for all member data.
    tick_style(
    // dim d = X, //!< X or Y axis (-1 if not assigned yet).
    const text_style& txtstyle = no_style, //!< Default text font style.
    double max = -10.,  //!< Maximum x value,  <min means we need to update the results.
    double min = +10., //!< Minimum x value., >max means we need to update the results
    unsigned int num_minor_ticks = 4) //!< Number of minor ticks between major ticks.
    : // Constructor.
		//dim_(d), // 1 or 2 D
		min_(min),
		max_(max),
		ticks_loc_(auto_place), // Value labels & ticks on the plot window, rather than on X or Y-axis only. Default -1 means left or bottom of plot window.
		values_side_(auto_outside), // Label values side for major ticks left (right or none).  2 /-2 means auto
		ticks_sides_(auto_inside), // Label values side for major ticks left (right or none).  2 /-2 means auto
		label_rotation_(horizontal), // Direction axis value labels written.
		major_grid_on_(false),  // Draw grid at major ticks.
		minor_grid_on_(false),// Draw grid at minor ticks.
		num_minor_ticks_(num_minor_ticks),
		major_ticks_frac_(-1),
		major_tick_color_(black), // line stroke color.
		major_tick_width_(1.5),
		major_tick_length_(4),
		minor_tick_color_(black), // line stroke color.
		minor_tick_width_(1),
		minor_tick_length_(2),
		major_grid_color_(svg_color(200, 220, 255)),
		major_grid_width_(1.), // pixels.
		minor_grid_color_(svg_color(200, 220, 255)),
		minor_grid_width_(0.5), // pixels.

		values_color_(black),
		value_precision_(3), // precision for tick value labels, usually 3 will suffice.
		value_ioflags_(std::ios::dec),  // IO formatting flags for the axis, eg. std::ios::scientific | std::ios::dec) & ~std::ios::fixed
		// Note that ALL the flags are set, overwriting any defaults, so std::dec is wise.
		// This should give the default 'normal' iosflags with neither fixed, scientific nor showpoint set.
		strip_0es_(true), // strip superflous zeros and signs.
		e_to_x10_(true), //
		label_max_space_(0.), // Space (estimated in SVG units) of longest label on axis adjusted for rotation.
		value_label_style_(txtstyle)
	{} // tick_style constructor.

  double label_length(double value)
  { //! Find the length of label (like "1.23E-5") for a value.
    // Needs to know the IO precision & flags for the axis, and if zeros are to be stripped, so can't be a free function.
    std::stringstream label;
    label.precision(value_precision_);
    label.flags(value_ioflags_);
    label << value; // "1.2" or "3.4e+000"...
    double r;
    if (strip_0es_) // Do want to strip unecessary e, +, & leading exponent zeros.
    {
      std::string stripped = strip_e0s(label.str());
      r = string_svg_length(stripped, value_label_style_, e_to_x10_);      // want x_or y_tick_values_style_ here!
      return r;
    }
    r = string_svg_length(label.str(), value_label_style_, e_to_x10_);
    return r;
  } // double label_length

  double longest_label(bool include_zero)
  { //! Calculate label_max_length_ with the longest value label as pixels,
    //! return the count of digits etc.
    if(values_side_ != 0) // ! none
    { // Show values by the tick as "1.2" or "3.4e+000"...
      double longest = 0;

      // Check length of label for the ticks on the positive side (right or above zero).
      double major_interval = (max_-min_)*major_ticks_frac_;
      major_interval = std::max(major_interval,std::abs(max_)*1e-5); // avoid inf loop
      for(double v = include_zero ? 0.0 : std::max(0.0,min_); v <= max_; v += major_interval)
      {
        if (v != 0. || ticks_loc_ != 0)
        { // Avoid a major tick at x == 0 where there *is* a vertical Y-axis line,
          // or avoid a major tick at y == 0 where there *is* a horizontal X-axis line.
          // (won't be a Y-axis line for 1-D, where both the zero tick & value label is always wanted).
          double l = label_length(v);
          if (l > longest)           longest = l;
        }
      }
      // Check length of label of the ticks on the negative side (left of zero).
      for(double v = include_zero ? 0.0 : std::min(0.0,max_); v >= min_; v -= major_interval)
      {
        if (v != 0. || ticks_loc_ != 0)
        { // Avoid a major tick at x == 0 where there *is* a vertical Y-axis line.
          // (won't be Y-axis line for 1-D where the zero tick is always wanted).
          // But no tick means no value label 0 either unless on_pl_window.
          double l = label_length(v);
          if (l > longest)          longest = l;
        }
      }
      return longest;
    }
    else
    {
      return 0;
    }
  } // longest_label()


  bool autoscale() const	 {	    return major_ticks_frac_<0.0 || max_<min_;   }//! \return true if to draw ticks up from horizontal X-axis line.

  tick_side ticks_side() const	 {	    return ticks_sides_;   }//! \return true if to draw ticks up from horizontal X-axis line.


  tick_style& ticks_side(tick_side t)
  { //! Set true to draw ticks up from horizontal X-axis line.
    ticks_sides_=t;    return *this; //! \return tick_style& to make chainable.
  }



  //int major_value_labels_side() const  {     return values_side_;	 }//! \return side for tick value labels: left (<0), none (==0) or right (>0).


}; // class tick_style

class box_element : public svg_style
{ //! \class svg::box_element Style of a rectangular box. (Used for boxplot image and plot window).
public:
    double margin_; //!< Marginal (pixels) space around the box (inside or out).
	//rect_element  rect_;//TODO, transfer from plot2d to here, remove style from parent

	//! Constructor to set parameters but provides defaults for all variables.
	box_element(
		const svg_color& scolor = black, //!< stroke color
		const svg_color& fcolor = blank, //!< fill color (white = no fill), blank => not drawn.
		double width = 1, //!< of border.
		double margin = 4. //!< Margin around box (SVG units, default pixels).
	) :
		svg_style(scolor, fcolor, width),
		margin_(margin)  //,rect_(-1,-1,-1,-1)
	{  }

// Member Functions definitions.

	box_element& stroke(const svg_color& color)  {   stroke_ = color;    return *this;   } //! Set (stroke) color for box outline.

	svg_color stroke()	 {	    return stroke_;   } //! \return (stroke) color for box outline.

	box_element& fill(const svg_color& color) {   fill_ = color;    return *this;  } //! Set fill color for box.

	svg_color fill()	 {	    return fill_;   }//! \return Fill color for box.

	box_element& width(double w) {  width_ = w;    return *this;   }//! Set width for box.

	double width()	 {	    return width_;   } //! \return width for box.

	box_element& margin(double w)  {    margin_ = w;    return *this;   } //! Set marginal (default pixels) space around the box (inside or out).


	double margin()	 {	    return margin_;   }//! \return marginal (default pixels) space around the box (inside or out).


}; // class box_element



inline const std::string strip_e0s(std::string s)
{ /* To remove redundant sign and leading zero(s) in exponent, for example, "1.2e+000" becomes "1.2"
    \details Used to work out the longest value label before update_internals.
    Should also be useful for values that spill over into exponent format
    'by accident' - when leading zeros are likely.
  */

  // An ugly hack but works...  (Could also do the same for uppercase E cases). (Considered doing a repeated strip but complicated).

  using std::string;
  size_t j = s.find("e+000");  if (j != string::npos)  s.erase(j, 5); // remove "e+000" completely, leaving no e... at all.

  j = s.find("e-000");  if (j != string::npos)    s.erase(j, 5); // remove entire "e-000".

  j = s.find("e+00");  if (j != string::npos)   s.erase(j + 1, 3); // remove "+00" from "e+009", leave d, so becomes e9.

  j = s.find("e-00");  if (j != string::npos)   s.erase(j+2, 2); // remove "00", leave "-" and any trailing d.

  j = s.find("e+0");  if (j != string::npos)    s.erase(j + 1, 2); // remove "+0", leave "dd"

  j = s.find("e-0");  if (j != string::npos)    s.erase(j+2, 1); // remove "-0", leave "-dd"



  j = s.find("000e");  if (j != string::npos)     s.erase(j, 3);

  j = s.find("00e");  if (j != string::npos)    s.erase(j, 2);

  j = s.find("0e");  if (j != string::npos)    s.erase(j, 1);

  j = s.find(".e");  if (j != string::npos)    s.erase(j, 1);
  if (s.back()=='e')    s.erase(s.size()-1, 1);



  return s; //! \return length of trimmed string (perhaps unchanged).
} // const std::string strip(double d)



inline int string_svg_length(const std::string& s, const text_style& style, bool eto10)
{
 int d = 0.; // Estimated or actual width of resulting svg string.
 for (std::string::const_iterator i = s.begin(); i != s.end(); i++)
 {
    if (*i == '&') // Start of Unicode 'escape sequence'
    {
       while ((*i != ';') && (i != s.end())) // In case mistakenly not null terminated.
          ++i; // Only count a Unicode string like &#x3A9; as 1 character (omega) wide.
    }
    if (*i == '<') // Embedded xml like <sub> or <super>.
    {
       while ((*i != '>') && (i != s.end())) // In case mistakenly not terminated.
           ++i; // Only count <...>; as NO characters wide.
       --d;
    }
    if (eto10 && *i == 'e') d+=2;
    ++d;
 }
 const double wh = 0.7; //!< font text width/height ratio. AQR: Orig was 0.7
  // http://www.w3.org/TR/SVG/text.html#FontSizeProperty
  // Font size is the height of the text's font, so width = wh * font_size.
  // Even after reading http://www.w3.org/TR/SVG/fonts.html,\n
  // unclear how to determine the exact width of digits, so an
  // arbitrary average width height ratio wh = 0.7 is used as a good approximation.
 return d * style.font_size() * wh;
}



inline const std::string e_to_x10(std::string s)
{
  using std::string;
  size_t j = s.find("e");
  if(s[1]=='e'&&s[0]=='1')
    s="10"+superscript(s.substr(2)); // remove "-0", leave "-dd"
  else
    if (j != string::npos)
      s=s.substr(0,j)+"⨯10"+superscript(s.substr(j+1)); // remove "-0", leave "-dd"
  return s; //! \return length of trimmed string (perhaps unchanged).
} // const std::string strip(double d)















// Caution: these two enum and ids must match because
// the enum value is used to index the array of id strings.
// void set_ids() copies all strings to matching image.get_g_element(i).id()
// Add any new id items to both!

// Order determines the painting order.
enum pl_doc_structure
{ //! \enum pl_doc_structure Plot document structure whose order controls the painting order, later ones overwriting earlier layers.
	g_BACKGROUND = 0, //! Must be zero to index array document_ids_[]
	g_AREA,          //! the smaller plot area (if used).
	g_Y_MINOR_GRID, //! Y minor grid.
	g_Y_MAJOR_GRID, //! Y major grid.
	g_X_MINOR_GRID, //! X minor grid.
	g_X_MAJOR_GRID, //! X major grid.
	g_DATA_LINES, //! Lines joining data points.
	g_DATA_POINTS, //! Normal data point markers.
	g_LIMIT_POINTS, //! 'At limit or NaN' data point markers.
	g_LEGEND_BOX, //! Legend box.
	g_LEGEND_POINTS, //! Legend data series point markers, circle, cross...
	g_LEGEND_TEXT, //! Legend text describing each data series.
	g_TITLE, //! Title of the whole plot.
	g_X_POINT_VALUES, //! X Data point value labels.
	g_Y_POINT_VALUES, //! Y Data point value labels.
	g_FUNCTIONS, //! Lines and curves, often to show a fit to the data.
	g_NOTES, //! Free text and shapes to annotate a plot.
	g_Y_AXIS, //! X axis line.
	g_X_AXIS, //! Y axis line.
	g_Y_MINOR_TICKS, //! Y minor ticks.
	g_X_MINOR_TICKS, //! X minor ticks
	g_Y_MAJOR_TICKS, //! Y major ticks.
	g_X_MAJOR_TICKS, //! X major ticks.
	g_X_TICKS_VALUES, //! X-axis tick values labels, for example 10, 20, 30 ...
	g_Y_TICKS_VALUES, //! Y-axis tick values labels, for example 1, 2, 3 ...
	g_Y_LABEL, //! Y axis text labels "length (cm)".
	g_X_LABEL, //! X axis text labels "height (m)".
	g_PLOT_ITEMS, //! Last enum value used as count of children (22).
};




} // namespace svg

#endif // BOOST_SVG_SVG_STYLE_HPP
