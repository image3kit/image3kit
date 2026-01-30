/*! \file

   \brief Boost.Plot SVG plot Implemention details.
   \details See svg.hpp etc for user functions.
      svg_tag.hpp defines all classes that can occur in the SVG parse tree.

   \author Jacob Voytko and Paul A. Bristow and Ali Q. Raeini
// Copyright Jacob Voytko 2007, 2008
// Copyright Paul A Bristow 2007, 2008, 2009, 2012
// Copyright Ali Q. Raeini 2018, 2019
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)
*/

 
class text_string; // Ancestor to both tspan and strings for the text_element class.
class tspan_element; // Within a text_element, adjust text and font properties.




#ifdef svg_EXTRA_SHAPES
#undef svg_EXTRA_SHAPES

class text_string //! An ancestor to both tspan and strings for the text_element class.
{                //!  This allows an array of both types to be stored in text_element. 
  protected:
	 std::string text_; //!< Actual text string for SVG text.
  public:
	 virtual ~text_string() {};  //! write functions output SVG commands.
	 virtual void write(std::ostream& o_str)  //! write functions output SVG commands.
	 {	 o_str << text_;	 }
	 text_string(const std::string& text): text_(text) 	{} //! Construct from text.
	 text_string(const text_string& rhs): text_(rhs.text_) 	{} //! Copy construct.
};
class tspan_element : public text_string, public svg_element // TODO also add public text_style rm style_
{ //! \class svg::tspan_element
	//! \brief tspan (not text) to be stored in text_string.
	//! \details See 10.5 tspan element http://www.w3.org/TR/SVG/text.html#TSpanElement
	private:
	  double x_;  //!< Absolute X position.
	  double y_;  //!< Absolute Y position.
	  double dx_; //!< Relative X position of a 1st single character of text.
	  double dy_; //!< Relative Y position of a 1st single character of text.
	  int rotate_; //!< Rotation of a 1st single character of text.
	  // A list of shifts or rotations for several characters is not yet implemented.
	  double text_length_;  //!< Allows the author to provide exact alignment.
	  //! dx_, dy_, and rotate_ can all be omitted, usually meaning no shift or rotation,
	  //! but see http://www.w3.org/TR/SVG/text.html#TSpanElement for ancestor rules.
	  //! but x_, y_, and text_length need a flag.
	  bool use_x_; //!> If true then use X absolute position.
	  bool use_y_; //!> If true then use Y absolute position.
	  bool use_text_length_; //!< If true then use user calculated length rather than SVG (not used).
	  text_style style_; //!< font variants.

 public:

	  tspan_element(const std::string& text, //!< Text string to display.
		      const text_style& style = no_style) //!< Text style (font).
	  :  text_string(text),x_(0), y_(0), dx_(0), dy_(0), rotate_(0),
		 text_length_(0), use_x_(false), use_y_(false),
		 use_text_length_(false), style_(style)
	  {  } //! Construct tspan element (with all defaults except text string).

	tspan_element(const tspan_element& rhs) //! copy constructor.
		 :
		 text_string(rhs),
		 x_(rhs.x_), y_(rhs.y_), dx_(rhs.dx_), dy_(rhs.dy_), rotate_(rhs.rotate_),  text_length_(rhs.text_length_), 
		 use_x_(rhs.use_x_), use_y_(rhs.use_y_),
		 use_text_length_(rhs.use_text_length_), style_(rhs.style_)
	  {  }
	 
		 // TODO all may need refactoring to separate declaration from definition

	  tspan_element& text(const std::string& text)	 { //! Set text string to use with SVG tspan command.
		 text_=text;	 return *this;	 }

	  tspan_element& dx(double dx)	 { dx_ = dx;	 return *this;	 } //! Relative X position of a 1st single character of text string to use with SVG tspan command.
		

	  tspan_element& dy(double dy)	 { dy_ = dy;	 return *this;	 } //! Relative Y position of a 1st single character of text string to use with SVG tspan command.
		

	  tspan_element& rotation(int rotation)	 { //!< Note implementation so far only rotates the 1st character in string.
		 //!< text_element rotation rotates the whole text string, so it *much* more useful.
		 rotate_ = rotation;	 return *this;	 }

	  tspan_element& x(double x)	 { //! Absolute X position of a 1st single character of text string to use with SVG tspan command.
		 x_ = x; 	 use_x_ = true;	 return *this;	 }

	  tspan_element& y(double y)
	  {//! Absolute Y position of a 1st single character of text string to use with SVG tspan command.
		 y_ = y; 	 use_y_ = true; 	 return *this;	 }

	  tspan_element& text_length(double text_length)	 { //! Set user estimate of text length (see http://www.w3.org/TR/SVG/text.html#TSpanElement TSPAN SVG Specification).
		 text_length_ = text_length;	 use_text_length_ = true;	 return *this;	 }

	  tspan_element& font_size(unsigned int size)	 { //! font size of 1st single character of text string to use with SVG tspan command.
		 style_.font_size(size); 		 return *this;	 }

	  tspan_element& font_family(const std::string& family)
	  {//! font family of 1st single character of text string to use with SVG tspan command.
		 style_.font_family(family);	 return *this;	 }

	  tspan_element& font_style(const std::string& style)	 { //! font style of 1st single character of text string to use with SVG tspan command.
		 //! font-style: normal | bold | italic | oblique
		 //! Examples: "italic"
		 //! http://www.croczilla.com/~alex/conformance_suite/svg/text-fonts-02-t.svg
		 style_.font_style(style);	  return *this;	 }

	  tspan_element& font_weight(const std::string& w)	 { //! font weight of 1st single character of text string to use with SVG tspan command.
	  //! svg font-weight: normal | bold | bolder | lighter | 100 | 200 .. 900
		 //! Examples: "bold", "normal"
		 //! http://www.croczilla.com/~alex/conformance_suite/svg/text-fonts-02-t.svg
		 //! tests conformance.  Only two weights are supported by Firefox, Opera, Inkscape.
		 style_.font_weight(w);	 return *this;	 }

	  tspan_element& fill_color(const svg_color& color)	 { //! Set fill color for a tspan element.
		 style_info_.fill_color(color);		 style_info_.fill_on(true); return *this;	 }

	  tspan_element& textstyle(const text_style& style)	 { //! Set text style (font) for a tspan element.
		 style_ = style;	 return *this;	 }


	  const text_style& textstyle()	 { //! \return text_style& to permit access to font family, size ...
		 return style_;	 }

	  const text_style& textstyle() const  {//! \return text_style& to permit access to font family, size (const version).
		 return style_;	 }

	  std::string text()	 {	 return text_;	 } //! Get text from a tspan element.
	
	  double x() 	{ 	return x_; 	} //! Get absolute X position for tspan element.

	  double y() 	{ 	return y_; 	} //! Get absolute Y position for tspan element.

	  double dx() 	{ 	return dx_; 	} //! Get relative X position for tspan element.


	  double dy() 	{ 	return dy_; 	} //! Get relative Y position for tspan element.

	  int rotation()
	  {  //! Get rotation for the next character for tspan element.
		 return rotate_;	 }

	  double text_length()	 { //! Get user estimated length for a text string.
		 return text_length_;	 }

	  bool use_style()	 { //! Get true if to use the estimated text string length.
		 return use_text_length_;	 }

	  unsigned int font_size()	 { //! Get the font size for tspan element (from its text style).
		 return style_.font_size();	 }

	  const std::string& font_family()	 { //! Get the font family for tspan element (from its text style).
		 return style_.font_family();	 }

	  const std::string& font_weight() const	  { //! Get the font weight for tspan element (from its text style).
		 return style_.font_weight();	 }

	  const std::string&  font_style()	 { //! Get the font style for tspan element (from its text_style).
		 return style_.font_style();	 }

	  const std::string&  font_style() const
	  { //! Get the font style for tspan element (from its text_style). const version.
		 return style_.font_style();	 }

	  svg_color fill_color()	 { //! Get the fill color for tspan element .
		 return style_info_.fill_color();	 }

	  bool fill_on()	 { //! Get true if to use fill color for tspan element .
		 return style_info_.fill_on();	 }

	  void write(std::ostream& os)	 { //! Output SVG XML for tspan_element
		 os << "<tspan";
		 write_attributes(os); // id & clip_path
		 style_info_.write_out(os); // fill, stroke, width...

		 // All of the conditional writes within tspan_element.

		 // First, all elements that can be tested based on their value.
		 if(rotate_ != 0)
			os << " rotate=\"" << rotate_ << "\"";

		 if(dx_!= 0)
			os << " dx=\"" << dx_ << "\"";

		 if(dy_!= 0)
			os << " dy=\"" << dy_ << "\"";


		 // Now, add all elements that can be tested with the flags.
		 if(use_x_ == true)
			os << " x=\"" << x_ << "\"";

		 if(use_y_  == true)
			os << " y=\"" << y_ << "\"";

		 if(use_text_length_ == true)
			os << " textLength=\"" << text_length_ << "\"";

		 if (style_.font_size() != 0)
			os << " font-size=\"" << style_.font_size() << "\"";

		 if (style_.font_family() != "") // Example: Arial.
			os << " font-family=\"" << style_.font_family() << "\"";

		 if (style_.font_style().size() != 0) // Example: italic.
			os << " font-style=\"" << style_.font_style() << "\"";

		 if (style_.font_weight().size() != 0) // Example: bold.
			os << " font-weight=\"" << style_.font_weight() << "\"";

		 if (style_.font_stretch().size() != 0)
			 os << " font-stretch=\"" << style_.font_stretch() << "\"";

		 if (style_.font_decoration().size() != 0)
			os << " font-decoration=\"" << style_.font_decoration() << "\"";

		 os << ">" << text_ << "</tspan>";	 } //   void write(std::ostream& os)
	}; // class tspan_element





  struct h_path: public path_point
  { /*! \struct svg::h_path
      \brief  Draws a horizontal line from the current point (cpx, cpy) to (x, cpy).
       which becomes the new current point. No y needed, start from current point y.
    */
    double x; //!< x horizontal SVG coordinate.
    void write(std::ostream& o_str)
    { //! Write horizontal line SVG command.
        o_str << (relative ? "h" : "H" ) << x << " ";
    }

    h_path(double x, bool relative = false)
    :  path_point(relative), x(x)
    { }
  }; // struct h_path


  struct v_path: public path_point
  { /*! \struct svg::v_path
        \brief Draws a vertical line from the current point (cpx, cpy) to (cpx, y).
        No x coordinate needed - use current point x.
    */
    double y; //!< Y vertical line SVG coordinate.
    void write(std::ostream& o_str)
    { //! Write vertical line SVG command.
        o_str << ( relative ? "v" :  "V" ) << y << " ";
    }

    v_path(double y, bool relative = false)
      :  path_point(relative), y(y)
    { //!< Constructor (defines all member variables).
    }
  }; // struct v_path


  struct c_path: public path_point
  { /*! \struct svg::c_path
     \brief Draws a cubic Bezier curve from the current point to (x, y) using (x1, y1).
     \details 8.3.5 The curve commands: C, Q & A.
    */
    double x1; //!< Middle of curve.
    double y1; //!< Middle of curve.
    double x2; //!< End point.
    double y2; //!< End point.
    double x; //!< Current (start point).
    double y; //!< Current (start point).

    void write(std::ostream& o_str)
    { //!< Write a cubic Bezier curve SVG XML to ostream.
        o_str << (relative ? "c" : "C" ) << x1 << "," << y1 << " " << x2 << "," << y2 << " "
        << x << "," << y << " ";
    } // void write(ostream&)

    c_path(double x1, double y1, double x2, double y2,
            double x, double y, bool relative = false)
      : path_point(relative), x1(x1), y1(y1), x2(x2), y2(y2), x(x), y(y)
    { }
  }; // struct c_path


  struct q_path: public path_point
  { /*! \struct svg::q_path
      \brief Draws a quadratic Bezier curve from the current point to (x,y).
       using (x1,y1) as the control point.
    */
    double x1; //!< quadratic Bezier curve control X coordinate.
    double y1; //!< quadratic Bezier curve control Y coordinate.
    double x; //!< quadratic Bezier curve end X coordinate.
    double y; //!< quadratic Bezier curve end Y coordinate.

    void write(std::ostream& o_str)
    { //!< Write a quadratic Bezier curve SVG XML to ostream.
      o_str << (relative ? "q" : "Q" ) << x1 << " " << y1 << " " << x << " " << y << " ";
    }

    q_path(double x1, double y1, double x, double y, bool relative = false)
      : path_point(relative), x1(x1), y1(y1), x(x), y(y)  {  } //! Constructor quadratic Bezier curve.

  }; //struct q_path




  struct t_path: public path_point
  { /*! \struct svg::t_path
      \brief Draws a quadratic Bezier curve from the current point to (x,y).
      \details see also s_path for a cubic Bezier curve.
    */
    double x; //!< SVG X coordinate.
    double y; //!< SVG Y coordinate.

    void write(std::ostream& o_str)
    { //! Write SVG command for a cubic Bezier curve.
        o_str << (relative ? "t" : "T") << x << "," << y << " ";
    }

    t_path(double x, double y, bool relative = false)
      : path_point(relative), x(x), y(y)
    {    } //! Constructor of path that draws a quadratic Bezier curve from the current point to (x,y)

  }; // struct t_path


  struct a_path : public path_point
  { /*! \struct svg::a_path
      \brief Draws a elliptical arc from the current point to (x,y),
        using two radii, axis rotation, and two control flags.
      \details See 8.3.8 The elliptical arc curve commands.!
        Useful for pie charts, etc.
     */
    double x; //!< X End of arc from current point.
    double y; //!< Y End of arc from current point.
    double rx; //!< X radius
    double ry; //!< Y radius
    double x_ax_rotation; //!< Any rotation of the X axis.
    bool large_arc; //!< true if arc >= 180 degrees wanted.
    bool sweep; //!< true if to draw in positive-angle direction

    void write(std::ostream& o_str)
    { //! Write elliptical arc path XML to ostream.
        o_str << (relative ? "a" : "A")	<< rx << "," << ry << " " << x_ax_rotation
			<< ((large_arc) ? 1 : 0) << "," << ((sweep) ? 1 : 0) << " " << x << "," << y << " ";
    }

    //! Construct elliptic arc path.
    a_path(double x, double y, double rx, double ry, double x_ax_rotation, bool large_arc = false, bool sweep = false, bool relative = false)
      : path_point(relative), x(x), y(y), rx(rx), ry(ry), x_ax_rotation(x_ax_rotation), large_arc(large_arc), sweep(sweep)
    {  }
  }; // struct a_path



  struct poly_path_point
  { /*! \struct svg::poly_path_point
      \brief polyline or polygon point coordinates (x, y)
      \details  9.6 polyline & 9.7 The 'polygon' element.
      */
    double x; //!< polygon or polyline path point X SVG coordinate.
    double y; //!< polygon or polyline path point Y SVG coordinate.
    // Polygon & polyline points are always absolute, never relative, and values have no preceeding letter like M or L,  So NOT derived from path_point.

    void write(std::ostream& o_str)
    { //! Output SVG XML,   //! Leading space is redundant for 1st after "points= ", but others are separators, and arkward to know which is 1st.
      o_str << " " << x << "," << y; // x, y separator comma for clarity.
    } // void write(std::ostream& o_str)

    poly_path_point(double x, double y)//! Construct a polygon or polyline path point from X and Y SVG coordinate.
      : x(x), y(y)     {    } 

    poly_path_point()//! Default constructor.
      : x(0.), y(0.)   {    } 


  }; // struct poly_path_point
 

  class polygon_element: public svg_element
  {  /*! \struct svg::polygon_element
     \brief The 'polygon' element defines a closed shape
     consisting of a set of connected straight line segments.

     \details http://www.w3.org/TR/SVG/shapes.html#PolygonElement
     The 'polygon'  element 9.9.7.
     A polygon is defined by including a 'path'  element
     which contains a points="(path data)"  attribute,
     where the d attribute contains the x, y coordinate pairs.
    */

    std::vector<std::unique_ptr<poly_path_point> > poly_points; //!< All the x, y coordinate pairs,
    //!< emplace_backed by calls of p_path(x, y).
  public:
    bool fill; //!< polygon to have fill color.

    //polygon_element(const polygon_element& rhs)
    //{ //! Copy constructor.
      //poly_points = (const_cast<polygon_element&>(rhs)).poly_points.release();    // 'empty' the vector of points.
    //}

    polygon_element() : fill(true)
    {    } //! Default constructor empty polygon (with fill on).


    polygon_element (double x, double y, bool f = true) : fill(f)
    { //! Constructor - One absolute (x, y) point only.
      //! Can add more path points using member function P.
      poly_points.emplace_back(new poly_path_point(x, y));
    }

    polygon_element (double x1, double y1, double x2, double y2, double x3, double y3, bool f = true)
      :      fill(f)
    { //! Constructor - Absolute (x, y) only. Used by triangle.
      poly_points.emplace_back(new poly_path_point(x1, y1));
      poly_points.emplace_back(new poly_path_point(x2, y2));
      poly_points.emplace_back(new poly_path_point(x3, y3));
    }

    polygon_element (double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, bool f = true)
      :      fill(f)
    { //! Constructor - Absolute (x, y) only. Used by rhombus.
      poly_points.emplace_back(new poly_path_point(x1, y1));
      poly_points.emplace_back(new poly_path_point(x2, y2));
      poly_points.emplace_back(new poly_path_point(x3, y3));
      poly_points.emplace_back(new poly_path_point(x4, y4));
    }

    polygon_element (double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double x5, double y5, bool f = true)
      :
      fill(f)
    { //! Constructor - Absolute (x, y) only. Used by pentagon.
      poly_points.emplace_back(new poly_path_point(x1, y1));
      poly_points.emplace_back(new poly_path_point(x2, y2));
      poly_points.emplace_back(new poly_path_point(x3, y3));
      poly_points.emplace_back(new poly_path_point(x4, y4));
      poly_points.emplace_back(new poly_path_point(x5, y5));
    }

    polygon_element (double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double x5, double y5, double x6, double y6, bool f = true)
      :
      fill(f)
    { //! Constructor - Six absolute (x, y) only. Used by hexagon.
      // Might be done more efficiently with fixed size std::array?
      poly_points.emplace_back(new poly_path_point(x1, y1));
      poly_points.emplace_back(new poly_path_point(x2, y2));
      poly_points.emplace_back(new poly_path_point(x3, y3));
      poly_points.emplace_back(new poly_path_point(x4, y4));
      poly_points.emplace_back(new poly_path_point(x5, y5));
      poly_points.emplace_back(new poly_path_point(x6, y6));
    }

    polygon_element (std::vector<poly_path_point>& points, bool f = true)
      :
      fill(f)
    { //! Constructor from vector of path points.
      poly_points.reserve(points.size()); // Since we know how many will be pushed.
      for(std::vector<poly_path_point>::iterator i = points.begin(); i != points.end(); ++i)
      {
        poly_path_point p = (*i);
        poly_points.emplace_back(new poly_path_point(p.x, p.y));
      }
    }

    // Member function to add more points to polygon.
    polygon_element& P(double x, double y)
    { //! Add another point (x, y) - absolute only.
      poly_points.emplace_back(new poly_path_point(x, y));     return *this;
    }

    void write(std::ostream& o_str)
    {  //! \verbatim SVG XML:
       //!     Example:   <polygon fill="lime" stroke="blue" stroke-width="10"  points="850,75 958,137.5 958,262.5 742,137.5" />
       //! \endverbatim

      o_str << "<polygon points=\"";
      for(auto i = poly_points.begin(); i != poly_points.end(); ++i)
         (*i)->write(o_str); //  x, y coordinates as " 1,2"
      o_str << "\"";
      write_attributes(o_str);
      style_info_.write_out(o_str);
      if(!fill) o_str << " fill = \"none\"";
      o_str<<"/>";
    }


  }; // class polygon_element


  class polyline_element: public svg_element
  { /*! \class svg::polyline_element
     \brief The 'polyline'  element: defines a set of connected straight line segments.
     \details
      http://www.w3.org/TR/SVG/shapes.html#PolylineElement
     9.6 The polyline  element: defines a set of connected straight line segments.
     Typically, polyline elements define open shapes.
     A polyline is defined by including a 'path'  element
     which contains a points="(path data)"  attribute,
     where the points attribute contains the x, y coordinate pairs.
     * perform an absolute moveto operation
       to the first coordinate pair in the list of points
     * for each subsequent coordinate pair,
       perform an absolute lineto operation to that coordinate pair.
     The advantage of polyline is in reducing file size,
     avoiding M and repeated L before x & y coordinate pairs.
     */

   public: //temporary for experimental gil
// private:
    std::vector<std::unique_ptr<poly_path_point> > poly_points; //!< All the (x, y) coordinate pairs,
    // emplace_back by calls of p_path(x, y).
  public:
    //bool fill; // not needed for polyline, unlike polygon.

    //polyline_element(const polyline_element& rhs)
    //{ //! copy constructor.
      //poly_points = (const_cast<polyline_element&>(rhs)).poly_points.release();
    //}

    polyline_element()
    { //! Construct an 'empty' line.
      //! Can new line path points add using polyline_element member function P.
    }

    polyline_element (double x1, double y1)
    { //! One (x, y) path point, absolute only.
      poly_points.emplace_back(new poly_path_point(x1, y1));
    }

    polyline_element (double x1, double y1, double x2, double y2)
    { //! Two (x, y) path points, absolute only.
      poly_points.emplace_back(new poly_path_point(x1, y1));
      poly_points.emplace_back(new poly_path_point(x2, y2));
    }

    polyline_element (std::vector<poly_path_point>& points)
    { //! Constructor from vector of path points.
      for(std::vector<poly_path_point>::iterator i = points.begin(); i != points.end(); ++i)
      {
        poly_path_point p = (*i);
        poly_points.emplace_back(new poly_path_point(p.x, p.y));
      }
    }

    // Member function to add new points to existing line.
    polyline_element& P(double x, double y)
    { //! Absolute (x, y) only, so Capital letter P.
      poly_points.emplace_back(new poly_path_point(x, y));
      return *this; //! \return polyline_element& to make chainable.
    }

    void write(std::ostream& o_str)
    { /*! \verbatim
          Output polyline info (useful for Boost.Test).
          Example: <polyline points=" 100,100 200,100 300,200 400,400"/>
          \endverbatim
      */
      o_str << "<polyline points=\"";
      for(std::vector<std::unique_ptr<poly_path_point> >::iterator i = poly_points.begin(); i!= poly_points.end(); ++i)
      {
        (*i)->write(o_str); //  x, y coordinates as " 1,2"
      }
      o_str << "\"";
      write_attributes(o_str);
      style_info_.write_out(o_str);
      o_str<<"/>";
    } // void write(std::ostream& o_str)

  }; // class polyline_element


#endif



#ifdef path_element_EXTRA_SHAPES
#undef path_element_EXTRA_SHAPES

    path_element& h(double x)
    { //! Line horizontal (relative).
      path.emplace_back(new h_path(x, true));	 return *this;
    }

    path_element& H(double x)
    { //! Line horizontal (absolute).
      path.emplace_back(new h_path(x, false));	 return *this;
    }

    path_element& v(double y)
    { //! Line vertical (relative).
      path.emplace_back(new v_path(y, true));	 return *this;
    }

    path_element& V(double y)
    {//! Line vertical (absolute).
      path.emplace_back(new v_path(y, false));	 return *this;
    }


    path_element& c(double x1, double y1, double x2, double y2, double x, double y)
    { //! Draws a cubic Bezier curve from the current point to (x,y) using (x1,y1).(relative).
      path.emplace_back(new c_path(x1, y1, x2, y2, x, y, true));	 return *this;
    }

    path_element& C(double x1, double y1, double x2, double y2, double x, double y)
    { //! Draws a cubic Bezier curve from the current point to (x,y) using (x1,y1).(absolute).
      path.emplace_back(new c_path(x1, y1, x2, y2, x, y, false));	 return *this;
    }

    path_element& q(double x1, double y1, double x, double y)
    {  //! Quadratic Curve Bezier (relative).
      path.emplace_back(new q_path(x1, y1, x, y, true));	 return *this;
    }

    path_element& Q(double x1, double y1, double x, double y)
    { //! Quadratic Curve Bezier (absolute).
      path.emplace_back(new q_path(x1, y1, x, y, false));	 return *this;
    }
    path_element& t(double x, double y)
    { //! Draws a quadratic Bezier curve from the current point to (x,y)(relative).
      path.emplace_back(new t_path(x, y, true));	 return *this;
    }

    path_element& T(double x, double y)
    { //! Draws a quadratic Bezier curve from the current point to (x,y)(absolute).
      path.emplace_back(new t_path(x, y, false));	 return *this;
    }
#endif


#ifdef g_element_EXTRA_SHAPES
#undef g_element_EXTRA_SHAPES
    //text_element& text(text_element& txt )
    //{ //! Add a new text element.
    //  //! \return A reference to the new child node just created.
    //  children.emplace_back(new text_element(txt));
    //  return *(static_cast<text_element*>(children.back().get()));
    //}

     emplace_back info about polygon shapes:

     Polygon for shapes with many vertices.
    polygon_element& polygon(double x, double y, bool f = true)
    { //! Add a new polygon element.
      //! \return A reference to the new child node just created.
      children.emplace_back(new polygon_element(x, y, f));
      return *(static_cast<polygon_element*>(children.back().get()));
    }

    //JVTODO: Replace with template version
    polygon_element& polygon(std::vector<poly_path_point>& v, bool f = true)
    { //! Add a new complete polygon element.
      //! \return A reference to the new child node just created.// emplace_back a complete many-sided polygon to the document.
      children.emplace_back(new polygon_element(v, f));
      return *(static_cast<polygon_element*>(children.back().get()));
    }
    //JVTODO: Replace with template version
    polyline_element& polyline(std::vector<poly_path_point>& v)
    {  //! Add a new complete polyline.
       //! \return A reference to the new child node just created.
      children.emplace_back(new polyline_element(v));
      return *(static_cast<polyline_element*>(children.back().get()));
    }

    polyline_element& polyline(double x, double y)
    { //! Add a new polyline element, but 1st point only, add others later with .P(x, y)...
      //! \return A reference to the new child node just created.
      children.emplace_back(new polyline_element(x, y));
      return *(static_cast<polyline_element*>(children.back().get()));
    }
    polygon_element& rhombus(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, bool f = true)
    { //! Add a new rhombus element.
      //! \return A reference to the new child node just created.
      children.emplace_back(new polygon_element(x1, y1, x2, y2, x3, y3, x4, y4, f = true));
      return *(static_cast<polygon_element*>(children.back().get()));
    }

    polygon_element& pentagon(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double x5, double y5, bool f = true)
    { //! Add a new pentagon element.
      //! \return A reference to the new child node just created.
      // emplace_back a complete pentagon to the document.
      children.emplace_back(new polygon_element(x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, f));
      return *(static_cast<polygon_element*>(children.back().get()));
    }

    polygon_element& hexagon(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double x5, double y5, double x6, double y6, bool f = true)
    { //! Add a new hexagon element.
      //! \return A reference to the new child node just created.
      // emplace_back a complete 6-sided star to the document.
      children.emplace_back(new polygon_element(x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6, f));
      return *(static_cast<polygon_element*>(children.back().get()));
    }

    polygon_element& polygon()
    { //! Add a new polygon element.
      //! \return A reference to the new polygon element just created.
      children.emplace_back(new polygon_element()); // Empty polygon,
      // to which poly_path_points can be added later using member function P.
      return *(static_cast<polygon_element*>(children.back().get()));
    }

    polyline_element& polyline()
    { //! Add a new polyline element.
      //! \return A reference to the new polyline element just created.
      children.emplace_back(new polyline_element()); // Empty polyline.
      return *(static_cast<polyline_element*>(children.back().get()));
    }

 
    polygon_element& triangle(double x1, double y1, double x2, double y2, double x3, double y3, bool f = true)
    { //! Add a new triangle element.
      //! \return A reference to the new child node just created.
      children.emplace_back(new polygon_element(x1, y1, x2, y2, x3, y3, f));
      return *(static_cast<polygon_element*>(children.back().get()));
    }
#endif

 


//Annotations:
#ifdef svplot_EXTRA_SHAPES
#undef svplot_EXTRA_SHAPES

	Derived& draw_note(double x, double y, std::string note, rotate_style rot /*= horizontal*/, align_style al/* = center_align*/, const svg_color& col /* black */, text_style& tsty/* = no_style*/)
	{ //! \brief Annotate plot with a  text string (perhaps including Unicode), putting note at SVG Coordinates X, Y.
		//! \details Defaults color black, rotation horizontal and align = center_align   Using center_align is recommended as it will ensure that will center correctly (even if original string is made much longer because it contains Unicode, for example Greek or math symbols, taking about 6 characters per symbol) because the render engine does the centering.
		g_element g = add_g_element(); // New group.
		g.style().fill_color(col); // Set its color
		g.push_back(new text_element(x, y, note, tsty, al, rot)); // Add to document image.
		// Could warn if X and/or Y outside - but even if OK, then text could still stray outside image.
		return *this;
	} // void draw_note()


	Derived& draw_svg_line(double x1, double y1, double x2, double y2, const svg_color& col)//!  \details Default color black.
	{ //! \brief Annotate plot with a line from SVG Coordinates X1, Y1 to X2, Y2.
		g_element gr = add_g_element(); // New group.
		gr.style().stroke_color(col);	 //g->style().width(w); // todo
		gr.push_back(new line_element(x1, y1, x2, y2));
		// No checks on X or Y - leave to SVG to not draw outside image. Could warn if X and/or Y outside ?
		return *this;
	} // void draw_line()



	Derived& draw_plot_line(double x1, double y1, double x2, double y2, const svg_color& col  )
	{ //! \brief Annotate plot with a line from user's Cartesian Coordinates X1, Y1 to X2, Y2.
		//!  \details For example, -10, -10, +10, +10, Default color black.
		//update_internals(); // To ensure the scale and shift are setup for transform.
		// It would be better to store the line (and curve and text) information like plot data sries to ensure that transform works correctly. This assumes that the notes, lines and curves are the last item before the write.
		transform_point(x1, y1);
		transform_point(x2, y2);
		g_element* g = &add_g_element(); // New group.
		g->style().stroke_color(col);
		g->push_back(new line_element(x1, y1, x2, y2));
		// No checks on X or Y - leave to SVG to not draw outside image. Actually we want to use clip_path for the plot area. Could warn if X and/or Y outside ?
		return *this;
	} // void draw_pl_line()


	Derived& draw_transform_curve(double x1, double y1, double x2, double y2, double x3, double y3, const svg_color& col)
	{ //! \brief Annotate plot with a line from user's Cartesian Coordinates X1, Y1 via X2, Y2 to X3, Y3.
		///	  \details For example, -10, -10, +10, +10, Default color black.
		//update_internals(); // To ensure the scale and shift are setup for transform.
		// It would be better to store the line (and curve and text) information like plot data sries to ensure that transform works correctly.
		// This assumes that the notes, lines and curves are the last item before the write.
		transform_point(x1, y1);
		transform_point(x2, y2);
		transform_point(x3, y3);
		g_element* g = &add_g_element(); // New group.
		g->style().stroke_color(col);
		g->push_back(new qurve_element(x1, y1, x2, y2, x3, y3));
		// No checks on X or Y - leave to SVG to not draw outside image.
		// Actually we want to use clip_path for the plot area.	 // Could warn if X and/or Y outside ?
		return *this;
	} // void draw_pl_curve 

	//! \endcond // DETAIL


#endif //ANOTATE_ENABLE


