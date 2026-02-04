/*! \file

   \brief Boost.Plot SVG plot Implemention details.
   \details See svg.hpp etc for user functions.
      svg_tag.hpp defines all classes that can occur in the SVG parse tree.

   \author Jacob Voytko and Paul A. Bristow and  Ali Q. Raeini
// Copyright Jacob Voytko 2007, 2008
// Copyright Paul A Bristow 2007, 2008, 2009, 2012
// Copyright Ali Q. Raeini 2018, 2019
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)
*/


#ifndef BOOST_SVG_TAG_HPP
#define BOOST_SVG_TAG_HPP


#include "svg_color.hpp"

#include <ostream>
#include <string>
#include <vector>


namespace svg {




const double dblepsilon = 2*std::numeric_limits<double>::min();
//double value_precision_ = 2*std::numeric_limits<double>::min();

inline std::string superscript(const std::string& text)
{ return "<tspan style=\"line-height:200%;font-size:80%\" dy=\"-5\">"+text+"</tspan><tspan dy=\"5\"></tspan>"; }
inline std::string  subscript(const std::string& text)
{ return "<tspan style=\"line-height:200%;font-size:80%\" dy=\"5\">"+text+"</tspan><tspan dy=\"-5\"></tspan>"; }

inline std::string  italic(const std::string& text)
{ return "<tspan style=\"font-style:italic;\">"+text+"</tspan>"; }
inline std::string  overline(const std::string& text)
{ return "<tspan text-decoration=\"overline\">"+text+"</tspan>"; }
inline std::string  underline(const std::string& text)
{ return "<tspan text-decoration=\"underline\">"+text+"</tspan>"; }


//! Default font chosen is a Unicode font like ['Lucida Sans Unicode] that
//! has the best chance of [symbols] being rendered corrrectly.
//! Used for title, legend, axes ... unless overridden by an explicit font specification.
const static char* default_font("Serif"); // Linux software replace this with a suitable (Free)serif font

 //! All assignement class members return references to this to make chainable.

/*!
 This is the style information for any group (g) tag.
 This could be expanded to include more data from the SVG standard.

 There are some strange effects for text on some browsers
 (Firefox especially) when only stroke is specified.
 fill is interpreted as black, and the font outline is fuzzy and bolder.
\verbatim
*   <g id="title" stroke="rgb(255,0,0)"> .. is red border and black fill.
*   (because created as a graphic not a builtin font?)
*   <g id="title" fill="rgb(255,0,0)"> .. is red sharp font.
*   <g id="title" stroke="rgb(255,0,0)" fill="rgb(255,0,0)"> red and red fill also fuzzy.
*   So for text, only specific the fill unless a different outline is really wanted.
*   Defaults for text provide a built-in glyph, for example for title:
*   <g id="title">
*     <text x="250" y="36" text-anchor="middle" font-size="18" font-family="Verdana">
*       Plot of data
*     </text>
*   </g>
*   and this is not a graphic.
 \endverbatim
*/

class svg_style
{ //! \class svg::svg_style Holds the basic SVG stroke, fill colors and width, and their switches.
  friend std::ostream& operator<< (std::ostream&, const svg_style&);
  //std::ostream& operator<< (std::ostream& os, const svg_style& s)
	char* extra_css_; //working with char* is dangrous, so keep private, std::string would be easier but uses much more memory

public: // Accesses only by set and get member functions below.
  // Private data member variables names end with _,
  // to permit use of names for set & get member functions.
  svg_color stroke_; //!< Color of SVG stroke (line or outline).
  svg_color fill_; //!< Color of SVG fill.
  float width_; //!< Width of line.  Only valid if > 0 
public:
  // Constructors declarations:
// Constructors.
  svg_style(const svg_color& stroke, const svg_color& fill, float width=1)
  : extra_css_(0),   stroke_(stroke), fill_(fill), width_(width)
  {  } //! Construct svg_style with specified fill and stroke colors, and width.
  svg_style(const svg_style& st)
  : extra_css_(0),   stroke_(st.stroke_), fill_(st.fill_), width_(st.width_)
  {  if(st.extra_css_) set_css(std::string(st.extra_css_));} //! Construct svg_style with specified fill and stroke colors, and width.

  
  
  svg_style()
  : //! Default svg_style has everything off.
  extra_css_(0), stroke_(svg_color(0, 0, 0, 0)), //! Stroke default is black. not shwon (opacity=0)
  fill_(blank), //! No fill color.
  width_(0) //! No width specified.
  {  } // Default constructor initialises all private data.


	~svg_style()  { if(extra_css_) delete[] extra_css_; extra_css_=0; }

 // Set svg_style member functions to set fill, stroke & width.
  svg_style& set_css(const std::string& str)
  { //! example set_css("stroke-dasharray=\"4\" ")
		if(extra_css_) delete[] extra_css_; //a bit low level, but this is a low level class anyway
		extra_css_ = new char[ str.size()+1 ];
		extra_css_[str.copy(extra_css_, str.size())]='\0';
		return *this;  //! \return svg_style& to make chainable.
  }
  const char* css() {return extra_css_;};

  // Set svg_style member functions
  // to set fill color and stroke color & width.
  svg_color fill_color() const
  { //! \return SVG fill color.
    return svg_color(fill_);
  }

  svg_color stroke_color() const  {   return svg_color(stroke_);  } //! \return SVG stroke color.

  float stroke_width() const  {    return width_;  } //! \return SVG stroke width.

  bool fill_on() const  { //! \return true if fill wanted.
    return fill_.o_>1; //0 is not shown, 1 is "none" (hence fully transparent)
  }

  svg_style& fill_on(bool is)  { //! Set fill is wanted.
    fill_.o_ = is*255;    return *this;   }

  bool stroke_on() const  { //! \return true if SVG stroke is on.
    return stroke_.o_;  }

  svg_style& stroke_on(bool is)  { //! Set true if SVG stroke is wanted.
    stroke_.o_ = is*255;    return *this;   }

  bool width_on() const  { //! \return true if to use SVG stroke width.
    return width_>dblepsilon;  }

  svg_style& width_on(bool is)  { //! Set true to use SVG stroke width.
    width_ *= is;    return *this;  }

  // Set svg_style member functions to set fill, stroke & width.
  svg_style& stroke_color(const svg_color& col)  { //! Set stroke color (and set stroke on).
    stroke_ = col;    return *this;   }

 
  svg_style& fill_color(const svg_color& col)  { //! Set fill color (and set fill on true, unless color is blank).
    //fill_on_ = ! is_blank(col); // If blank fill is off or "none".
    fill_ = col;    return *this;   }

  svg_style& stroke_width(double width)  { //! Set stroke width (and set width on).
      width_ = width;      return *this;  }


  void write_out(std::ostream& os)
  { //! Write any stroke, fill colors and/or width info to SVG XML document.
    if(stroke_.o_)
    { // (Note: start with space but no terminating space)
        os << " stroke=\"";
        stroke_.write(os);
        os << "\"";
    }
    if(fill_.o_==255) //  && (fill_ != blank))
    { // Don't add fill info if color is blank.
        os << " fill=\"";
        fill_.write(os);
        os << "\"";
    }
    else if(fill_.o_>1) //  && (fill_ != blank))
    { // Don't add fill info if color is blank.
        os << " fill=\"";
        fill_.write(os);
        os << "\"";
        os << " fill-opacity=\"";
        os << (fill_.o_-1.0)/254.0f;//! accepted numbers are 1-255, but then mapped to 0-1.0
        os << "\"";
    }
    else if(fill_.o_) //==1
    {
        os << " fill=\"none\"";
	}

    if(width_ > dblepsilon) // We never want a 0 (or <0) width output?
    {   os << " stroke-width=\"" << width_<< "\"";    }
    
    if(extra_css_) // We never want a 0 (or <0) width output?
    {   os <<" "<< extra_css_;    }

   /*! \details Example output: \<g id="yMinorTicks" stroke="rgb(0,0,0)" stroke-width="1"\>
   */
 } // void write

}; // class svg_style



//! \class svg::text_style
class text_style//!  \brief Font family, font size, weight, style, stretch & decoration.
{
  //friend std::ostream& operator<< (std::ostream&, const text_style&);
  //friend bool operator== (const text_style&, const text_style&);
  //friend bool operator!= (const text_style&, const text_style&);

 public: // Or private?
  float font_size_; //!< Font size (SVG units, default pixels).
  std::string font_family_; //!< Font family, examples: "Arial", "Times New Roman", "Verdana", "Lucida Sans Unicode".
  std::string weight_; //!< Font style, examples: "bold", "normal".
  std::string font_style_; //!< Font weight, examples: normal | italic | oblique.
  std::string stretch_; //!< Font stretch, examples: normal | wider | narrower.
  std::string decoration_; //!< Font decoration, examples: "underline" | "overline" | "line-through".
  //std::string baseline_; //!< dominant-baseline, Example middle, central bottom , top, not supported in inkskape

 public:

//! Default constructor only sets font size = 12, and leaves other font details as SVG defaults.
  text_style( //!< Constructor to allow all text style (font etc) to be set.
    float size = 16, //!<  font size (default: 12 pixels).
    const std::string& font = default_font, //!< Examples: "Arial", "Times New Roman", "Verdana", "Lucida Sans Unicode"
    const std::string& weight = "", //!< font weight: Examples: "bold", "normal"
    const std::string& font_style = "", //!< font-style: normal | italic | oblique
    const std::string& stretch = "", //!< font-stretch: normal | wider | narrower ...
    const std::string& decoration = ""
    //,const std::string& baseline = ""
    ) //!< No decoration.
  : // Constructor.
	  font_size_(size),
	  font_family_(font),
	  weight_(weight),
	  font_style_(font_style),
	  stretch_(stretch),
	  decoration_(decoration)
	  //,baseline_(baseline)
  {  } // text_style default constructor, defines defaults for all private members.


  float font_size() const   {//! \return  font size (svg units, usually pixels).
						return font_size_;  }

  text_style& font_size(float i)  { //! Set font size (svg units usually pixels) default 10.
						font_size_ = i;    return *this;   }

  //const std::string& font_family() const  { //! \return  font family as string.
						//return font_family_;  }

  text_style& font_family(const std::string& s)  { //! Set font family, for example: "Arial", "Times New Roman", "Verdana", "Lucida Sans Unicode".
						font_family_ = s;    return *this;  }

  //const std::string& font_style() const  { //! \return  font style. font-style: normal | bold | italic | oblique.   Example "normal" is default.
						//return font_style_; }

  text_style& font_style(const std::string& s)  { //! Set font style.  Example: my_text_style.font_style("italic");\n
						font_style_ = s;    return *this;  }

  //const std::string& font_weight() const  {  //! Set font weight.  Example: my_text_style.font_style("bold");\n
						//return weight_;  }

  text_style& font_weight(const std::string& s)  { //! svg font-weight: normal | bold | bolder | lighter | 100 | 200 .. 900.  Examples: "bold", "normal"
						weight_ = s;    return *this;  }

  //const std::string& font_stretch() const  { //! \return font stretch, for example: normal | wider | narrower .
    //return stretch_;  }

  text_style& font_stretch(const std::string& s)  { //! Examples: "wider" but implementation by browsers varies.  font-stretch: normal | wider | narrower ...
    stretch_ = s;    return *this;  }

  //const std::string& font_decoration() const  { //! \return  font decoration.
    //return decoration_;  }

  text_style& font_decoration(const std::string& s)  { //! Set font decoration.      Examples: "underline" | "overline" | "line-through"
    decoration_ = s;    return *this;   }

  //text_style& baseline(const std::string& s)  { //! Set font decoration.      Examples: "underline" | "overline" | "line-through"
    //baseline_ = s;    return *this;   }//! \return assignement functions return reference to text_style to make chainable.




}; //   end class text_style


const static text_style no_style = text_style(12); //!< Text style that uses all constructor defaults.



// Forward declarations of classes defined in this module.

class g_element;  // (group element)
class svg_element; // svg_element is base class for:
class text_element; // text with position, size, font, (& styles) & orientation.
class rect_element; // clipping path restricts the region to which paint can be applied.
class circle_element; // Represents a single circle.
class ellipse_element; // Represents a single ellipse.
class line_element; // Represents a single line.
struct path_point; // Base class for m_path, z_path, q_path, h_path, v_path, c_path, s_path.
struct poly_path_point; // for polyline & polygon
class polygon_element; // closed shape consisting of a set of connected straight line segments.
class polyline_element; // a set of connected straight line segments.
class path_element; // d= moveto, lineto...
struct m_path; // moveto coordinates (x, y), outputs "M1.2,3.4"
struct l_path; // lineto coordinates (x, y).
struct z_path; // z indicates a closepath.
struct h_path; // Draws a horizontal line from the current point (cpx, cpy) to (x, cpy).
struct v_path; // Draws a vertical line from the current point (cpx, cpy) to (cpx, y).
struct c_path; // Draws a cubic Bezier curve from the current point to (x,y) using (x1,y1).
struct q_path; // Draws a quadratic Bezier curve from the current point to (x,y).
struct s_path; // Draws a cubic Bezier curve from the current point to (x,y).
struct t_path; // Draws a quadratic Bezier curve from the current point to (x,y).
struct a_path; // Draws a elliptical arc from the current point to (x,y).


class svg_element
{ /*! \class svg::svg_element
	   \brief svg_element is base class for all the leaf elements.
	   \details
	   rect_element, circle_element, line_element, text_element,
	   polygon_element, polyline_element, path_element, clip_path_element,
	   g_element.\n

	   g_element ('g' element is a container element
	   for grouping together related graphics elements).\n
	   See http://www.w3.org/TR/SVG/struct.html#NewDocument 5.2.1 Overview.
	*/

	protected:
	svg_style   style_info_; //!< Colors fill, stroke, width, get by function style.
	std::string id_name_; //!< SVG id name, set & get by function id.
	std::string clip_id_; //!< SVG clip path name, set & get by function clip_id.

	void write_attributes(std::ostream& s_out)
	{ //! Output group_element id and clip-path.
	  if(id_name_.size() != 0)
	  { // Example: id="imageBackground"
		s_out << " id=\"" << id_name_ << "\""; // Prefix with space.
	  }

	  if(clip_id_.size() != 0)
	  { // Example: clip-path="url(#pl_window)"
		s_out << " clip-path=\"url(#" << clip_id_ << ")\""; // Prefix with space.
		// Might be nicer to suffix with newline - but after the >
	  }
	  // should transform be here allow translate and rotate?
	  /*! \details
		Classes inherited from svg_element add other references, 5.3.1, like color, fill, stroke, gradients...
		*/
	  /*
		\verbatim
		  Example id: <g id="yMinorGrid" ></g>
		  Example class: <g class="grid_style"></g>
		  Example URI: fill="url(#Gradient01) // local URL
		\endverbatim
	  */
	} // void write_attributes(std::ostream& s_out)

	public:

	svg_element(const svg_style& style_info)
				: style_info_(style_info)   {    } //! Constructor with some defaults.


	svg_element()    {     }//! Default constructor.


	virtual void write(std::ostream& rhs) = 0; //!< write functions output SVG commands.

	virtual ~svg_element()
	{ //! destructor.
	}

	bool operator==(const svg_element& lhs)
	{ //! Compare svg_elements, useful for Boost.Test.
	  return lhs.id_name_ == id_name_;
	}

	 //! Compare svg_elements for inequality, useful for Boost.Test.
	bool operator!=(const svg_element& lhs) 	{ 	return lhs.id_name_ != id_name_; 	}

	// Set and get member functions.
	svg_style& style()
	{ //! \return  reference to svg_style to provide indirect access to colors & width via style().stroke_color(), fill_color(), width()
	  return style_info_;
	}

	const svg_style& style() const
	{ //! \return  reference to const svg_style to provide indirect access to colors & width via style().stroke_color(), fill_color(), width() (const version).
	  return style_info_;
	}

	void id(const std::string& id_name)
	{ //! Provide a unique name for an element. Example: id="plotBackground"

	  /*! \details
		See http://www.w3.org/TR/SVG/struct.html#IDAttribute
		5.10.1 Attributes common to all elements: id and xml:base
		The id and xml:base attributes are available on all SVG elements:
		Attribute definitions:
		id = "name"
	  */
	  id_name_ = id_name;
	}

	std::string id()
	{ //! \return  the unique name for an element, for example id() ="plotBackground".
	  return id_name_;
	}

	//void class_id(const std::string& class_id)
	//{ //! Class class id, non-unique identifier for an element.
	  ///*! \details
		//http://www.w3.org/TR/2001/REC-SVG-20010904/styling.html#ClassAttribute
		//6.12 Attributes common to all elements: id and xml:base
		//Example: class="info"
	  //*/
	  //class_name_ = class_id;
	//}

	void clip_id(const std::string& id)
	{ //! Set name of a clip path, for example: g_ptr.clip_id(pl_window_clip_);
	  clip_id_ = id;
	}

}; // class svg_element


// Derived elements whose write member functions
// output SVG XML for line, circle, rectangle, text...
// Reminder: Within a literal C string, \"  is needed to output a " ;-)


// Represents a straight line
class line_element: public svg_element
{ /*! \class svg::line_element
	\brief Line from (x1, y1) to (x2, y2).
	/details Straight line from SVG location (x1, y1) to (x2, y2).

*/
public:
//  private:
double x1_; //!< Line from (x1_, x2_) to (y1_, y2_)
double y1_; //!< Line from (x1_, x2_) to (y1_, y2_)
double x2_; //!< Line from (x1_, x2_) to (y1_, y2_)
double y2_; //!< Line from (x1_, x2_) to (y1_, y2_)

public:
line_element(double x1, double y1, double x2,  double y2)
  :   x1_(x1),  y1_(y1), x2_(x2), y2_(y2)
{ //! Constructor assigning all line_element private data.
}

line_element(double x1, double y1,
			 double x2, double y2,
			 const svg_style& style_info       ) 
			: svg_element(style_info),
			x1_(x1), y1_(y1), x2_(x2),y2_(y2)
{  } //! Constructor assigning all line_element private data, and also inherited svg_element data.



void write(std::ostream& rhs)
{ //! output line from (x1_, y1_) to (x2_, y2_) by  writing XML SVG command to draw a straight line.
	//! \verbatim Example: <line x1="5" y1="185" x2="340" y2="185"/> \endverbatim 
  rhs << "<line x1=\"" << x1_ << "\" y1=\"" << y1_  << "\" x2=\"" << x2_ << "\" y2=\"" << y2_ << "\"/>";
}
}; // class line_element


// Represents a curve (quadratic)
class qurve_element: public svg_element
{ /*! \class svg::qurve_element
	\brief Quadratic Bezier curved line from (x1, y1) control point (x2, y2) to (x3, y3).
	\details Note x2 is the Bezier control point - the curve will \b not pass thru this point.
*/
 public: //temporary for experimental gil
	//  private:
	double x1_; //!< Quadratic curved line from (x1_, y1_) control point (x2_, y2_) to (y3_, y3_).
	double y1_; //!< Quadratic curved line from (x1_, y1_) control point (x2_, y2_) to (y3_, y3_).
	double x2_; //!< Quadratic curved line from (x1_, y1_) control point (x2_, y2_) to (y3_, y3_).
	double y2_; //!< Quadratic curved line from (x1_, y1_) control point (x2_, y2_) to (y3_, y3_).
	double x3_; //!< Quadratic curved line from (x1_, y1_) control point (x2_, y2_) to (y3_, y3_).
	double y3_; //!< Quadratic curved line from (x1_, y1_) control point (x2_, y2_) to (y3_, y3_).

 public:
	// bool fill; now inherited from parent svg class.
	qurve_element(double x1, double y1, double x2,  double y2, double x3,  double y3)
	  :   x1_(x1), y1_(y1), x2_(x2), y2_(y2),  x3_(x3), y3_(y3)
	{ //!< Quadratic curved line constructor (info inherited from parent svg_element class).
	}

	qurve_element(double x1, double y1,
				 double x2, double y2, // Control point - will not pass thru this point.
				 double x3, double y3,
				 // Inherited from svg_element:
				 const svg_style& style_info)
				:  svg_element(style_info),
				x1_(x1),y1_(y1), x2_(x2),y2_(y2), x3_(x3),y3_(y3)
			   
	{ //!< Quadratic curved line constructor, including svg_element info.
	}

	void write(std::ostream& o_str)
	{ /*! output quadratic curved line from (x1_, y1_) control point (x2_, y2_) to (x3_, y3_)
	   \details
		  \verbatim Example:
		  \endverbatim
	  */
	  o_str << "<path d=\"M" << x1_ << "," << y1_
		  << " Q" << x2_ << "," << y2_ << " " // Control point - will not pass thru this point.
		  << x3_ << "," << y3_ <<"\"";
	  if(style_info_.fill_on() == false)
		  o_str << " fill = \"none\"";
	  o_str<<"/>";
	}
}; // class qurve_element


class rect_element : public svg_element
{ /*! \class svg::rect_element
	\brief Rectangle from top left coordinate, width and height.
	\details
	 Represents a single rectangle.
	 http://www.w3.org/TR/SVG/shapes.html#RectElement
*/

 public: //temporary for experimental gil

	//  private:
	double x_; //!< X-axis coordinate of the side of the rectangle which has the smaller x-axis coordinate value.
	double y_; //!< Y-axis coordinate of the side of the rectangle which has the smaller y-axis coordinate value.
	//!< So (0, 0) is top left corner of rectangle.
	double width_; //!< x + width is top right.
	double height_; //!< y + height is bottom left.
	//!< x + width and y + height is bottom right.
	public:

	rect_element(double x, double y, double w, double h)
	  : x_(x), y_(y), width_(w), height_(h)
	{ //! Constructor defines all private data (no defaults).
	}

	rect_element(double x, double y, double w, double h,
				 const svg_style& style_info)                // Inherited from svg_element:
	  :        svg_element(style_info),     x_(x), y_(y), width_(w), height_(h)
	{    } //! Constructor defines all private data (inherites info from svg_element).

	rect_element(rect_element & rect)
				: svg_element(rect),
			 x_(rect.x_), y_(rect.y_), width_(rect.width_), height_(rect.height_)
	{ //! Constructor defines all private data (inherites info from svg_element).
	}
	void operator=(const rect_element & rect)
	{     //! Constructor defines all private data (inherites info from svg_element).
	   x_=(rect.x_); y_=(rect.y_); width_=(rect.width_); height_=(rect.height_);
	}

	double x() const
	{ //! x-axis coordinate of the side of the rectangle which has the smaller x-axis coordinate value.
	  return x_;
	}

	double y() const
	{ //! y-axis coordinate of the side of the rectangle which has the smaller y-axis coordinate value.
	  return y_;
	}

	double width() const
	{ //! x + width is top right.
	  return width_;
	}

	double height() const
	{ //! y + height is bottom left.
	  return height_;
	}

	void write(std::ostream& rhs)
	{ /*! \verbatim
		Output SVG XML for rectangle.
	   For example: <rect  x="0" y="0"  width="500"  height="350"/>
	   \endverbatim
	   */
	  rhs << "<rect";
	  write_attributes(rhs); // id (& clip_path)
	  rhs << " x=\"" << x_ << "\"  y=\"" << y_ << "\""
		<< " width=\"" << width_ << "\"  height=\""<< height_<< "\"/>";
	}

	bool operator==(const rect_element& lhs) const
	{ //! Comparison (useful for Boost.Test).
	  return (lhs.x() == x_) && (lhs.y() == y_) &&  (lhs.width() == width_) && (lhs.height() == height_);
	}
	bool operator!=(const rect_element& lhs) const
	{ //!< Comparison rect_elements (useful for Boost.Test).
	  return (lhs.x() != x_) || (lhs.y() != y_) ||  (lhs.width() != width_) || (lhs.height() != height_);
	}
}; // class rect_element



// Represents a single circle
class circle_element : public svg_element
{/*! \class svg::circle_element
	\brief Circle from center coordinate, and radius.
	\details
	 Represents a single circle.
	 http://www.w3.org/TR/SVG/shapes.html#CircleElement
*/
 private:
		double x;
		double y;
		double radius;
 public:
	circle_element(double x,  //!< coordinate X of center of ellipse.
	  double y, //!< coordinate Y of center of ellipse.
	  double radius = 5) //!< radius of ellipse.
	  : x(x), y(y), radius(radius)
	{ //! Constructor defines all private data (default radius only).
	}

	circle_element(double x, double y, double radius,//! Constructor defines all private data.
				 const svg_style& style_info        )
	  : svg_element(style_info),  x(x), y(y), radius(radius)
	{ }


	void write(std::ostream& rhs)
	{ /*! Output SVG XML
	\verbatim
	   Example: <circle cx="9.78571" cy="185" r="5"/>
	\endverbatim
	*/
	  rhs << "<circle";
	  write_attributes(rhs);
	  rhs << " cx=\"" << x << "\" cy=\"" << y << "\" r=\"" << radius << "\"/>";
	}
}; // class circle_element


class ellipse_element : public svg_element
{ /*! \class svg::ellipse_element
  \brief Ellipse from center coordinate, and radius.
  \details
  Represents a single ellipse.
  http://www.w3.org/TR/SVG/shapes.html#EllipseElement
  9.4 The 'ellipse'  element.
  Default is 'horizontal' but can be rotated.
  */
	private:
		double cx_; //!< coordinate x of center of ellipse, default 0.
		double cy_; //!< coordinate y, default 0.
		double rx_; //!< radius x, default 4 pixels.
		double ry_; //!< radius y, default 8 pixels.
		double rotate_; //! rotation in degrees from horizontal (default 0.).
		// Only hacked in - should be in attributes?
public:
	ellipse_element(double cx, //!< coordinate X of center of ellipse.
		double cy, //!< coordinate Y  of center of ellipse.
		double rx = 4, //!< X radius default.
		double ry = 8) //!< Y radius default.
	:	cx_(cx), cy_(cy), rx_(rx), ry_(ry), rotate_(0.)
	{ //!< Constructor defines all private data (with default radii).
	}

 ellipse_element(double cx, //!< coordinate X of center of ellipse.
	double cy, //!< coordinate Y of center of ellipse.
	double rx, //!< radius X.
	double ry, //!< radius Y.
	const svg_style& style_info) //!< name of clip path.
	:
		svg_element(style_info),
		cx_(cx), cy_(cy), rx_(rx), ry_(ry), rotate_(0.)
 { //!< Constructor sets ellipse and its style (defaults define all private data).
 }
 ellipse_element(
	  double cx, //!< coordinate X of center of ellipse.
	  double cy, //!< coordinate Y of center of ellipse.
	  const svg_style& style_info)
	:	svg_element(style_info),	cx_(cx), cy_(cy), rx_(4), ry_(8), // 4 and 8 are the same defaults used above.
		rotate_(0.) // Horizontal.
 {  //!< Constructor that also includes style, id, class and clip.
 }

 void write(std::ostream& os)
 { //!  Output SVG XML for ellipse.  Example: \<ellipse rx="250" ry="100" fill="red"  /\>
	os << "<ellipse";
	write_attributes(os);
	if(rotate_ != 0)
		  os << " transform= \"" << " rotate=(" << rotate_ << ")\"";// Should this be in atttributes?
	
	os	<< " cx=\"" << cx_ << "\" cy=\"" << cy_ << "\""
		<< " rx=\"" << rx_ << "\" ry=\"" << ry_  << "\"/>";
 }
}; // class ellipse_element


enum align_style
{ //! \enum align_style Represents a single block of text, with font & alignment.
 left_align, //!< Align text to left.
 right_align, //!< Align text to right.
 center_align //!< Center text.
};


//! \enum rotate_style Rotation of text (in degrees clockwise from horizontal).
enum rotate_style
{
  // Might include more steps for others too.
  horizontal = 0, //!< normal horizontal left to right, centered.
  //uphill = -45, //!< slope up. seems bit steep
  slopeup = -30, //!< slope up.
  uphill = -45, //!< 45 up.
  steepup = -60, //!< up near vertical.
  upward = -90, //!< vertical writing up.
  backup = -135, //!< slope up backwards - upside down!
  righttoleft= -180, //!< horizontal to left.
  lefttoright = 360, //!< horizontal to right.
  slopedownhill = 30, //!< 30 gentle slope down.
  downhill = 45, //!< 45 down.
  steepdown = 60, //!<  60 steeply down.
  downward = 90,  //!< vertical writing down.
  backdown = 135, //!< slope down backwards.
  upsidedown = 180 //!< upside down!  (== -180)
};



class text_element : public svg_element
{ /*! \class svg::text_element
		\brief Holds text with position, size, font, (& styles) & orientation.
		\details
		Not necessarily shown correctly (or nicely) by all browsers, alas.
		SVG Coordinates of 1st character EM box, see http://www.w3.org/TR/SVG/text.html#TextElement 10.2

		So any text with y coordinate = 0 shows only any roman lower case descenders!\n

		\verbatim
		  (Text may contain embedded xml Unicode characters for Greek, math etc, for example: &#x3A9;).
		\endverbatim
		\n
		int size; // " font-size = 12"
		http://www.w3.org/TR/SVG/text.html#CharactersAndGlyphs
		std::string font;  // font-family: "Arial" | "Times New Roman" | "Verdana" | "Lucida Sans Unicode"
		"sans", "serif", "times"
		http://www.w3.org/TR/SVG/text.html#FontFamilyProperty
		10.10 Font selection properties\n

		std::string style_; // font-style: normal | bold | italic | oblique
		std::string weight; // font-weight: normal | bold | bolder | lighter | 100 | 200 .. 900
		std::string stretch; // font-stretch: normal | wider | narrower ...
		std::string decoration; // // "underline" | "overline" | "line-through"
		Example:
		\verbatim
		  <text x="250" y="219.5" text-anchor="middle"  font-family="verdana" font-size="12">0 </text>
		\endverbatim

  */
private: // Access only via member functions below:
  double x_; //!< Left edge.
  double y_; //!< Bottom of roman capital character.
  std::vector<std::string> data_; //!< Stores all of the containing tests and tspans. tspans should be converted to std::string.
  text_style style_; //!< font variants.
  align_style align_; //!< left_align, right_align, center_align
  rotate_style rotate_; //!< horizontal, upward, downward, upsidedown

  void generate_text(std::ostream& os)
  {
	 for(auto i = data_.begin(); i != data_.end();  ++i)
		os << *i;
 }
public:

	text_element(
		 //!< Coordinates of 1st character EM box, see
		 //!< http://www.w3.org/TR/SVG/text.html#TextElement 10.2
		 double x = 0., //!< X = Left edge.
		 double y = 0., //!< Y =  Bottom left of (western) character (roman capital).
		 //!< So any text with Y coordinate = 0 shows only the roman lower case descenders!
		 //!< One must increase Y to allow for the height (font size) of the character.
		 const std::string text = "", //!< Text string to output (may include Unicode string like "&#x221A;" for square root symbol.
		 text_style ts = no_style, //!< Text font style,default left to SVG defaults.
		 align_style align = left_align, //!< Alighment of text, left, center or right, default left_align.
		 rotate_style rotate = horizontal) //!< orientation of text, default horizontal.
	: // Constructor.
		x_(x), y_(y), // location.
		data_(),
		style_(ts), // Simpler to include all these as members?
		align_(align),
		rotate_(rotate)	 
	{ //! text_element Default Constructor, defines defaults for all private members.
		if(text.size()) data_.emplace_back(text); // Adds new text string.
	}

  //text_element(const text_element& rhs)
  //: 	svg_element(rhs), x_(rhs.x_), y_(rhs.y_), style_(rhs.style_), align_(rhs.align_), rotate_(rhs.rotate_)	 { //! Copy constructor.
		//for(auto& dta:rhs.data_) data_.emplace_back(dta);
	//}

  //text_element& operator=(text_element&& rhs)	 { //! Assignment operator.
	 //x_ = rhs.x_;
	 //y_ = rhs.y_;
	 //data_.clear(); // Copy data_
	 //for(auto& dta:rhs.data_) data_.emplace_back(dta);
	 //style_ = rhs.style_;
	 //align_ = rhs.align_;
	 //rotate_ = rhs.rotate_;
	 //return *this;	 }

  text_style& textstyle()	 { //! Get text style for font size, family, decoration ...
	 return style_;	 }

  const text_style& textstyle() const  { //! Get text style for font size, family, decoration ...
						return style_;	 }

  text_element& textstyle(text_style& ts)	 { //! Set text style for font size, family, decoration ...
						style_ = ts;	 return *this;	 }

  text_element&  alignment(align_style a)  { //! left_align, right_align, center_align
							align_ = a;	 return *this;	 } // TODO Change name to align????

  align_style alignment()	 { //! left_align, right_align, center_align
							return align_;	 }

  text_element&  rotation(rotate_style rot)// TODO Change name to rotate???
  { //! Degrees: horizontal  = 0, upward = -90, downward, upsidedown \n Generates: transform = "rotate(-45 100 100 )"
	 rotate_ = rot;	 return *this;	 }

  rotate_style rotation() const  { //! \return rotation of text element.
	 return rotate_;	 }

  // set functions return *this to be chainable, for example:  // my_text_element.style(no_style).x(999).y(555).alignment(right_align).rotation(vertical);

  text_element& x(double x)	 { //! x coordinate of text to write.
	 x_ = x;	 return *this;	 }

  double x() const  { //! x coordinate of text to write.
	 return x_;	 }

  text_element& y(double y)	 { //! y coordinate of text to write.
	 y_ = y;	 return *this;	 }

  double y() const  { //! y coordinate of text to write.
	 return y_;	 }

  void text(const std::string& t)	 { //! Get tspan text string to write.
	 data_.emplace_back(t);	 }

  //tspan_element& tspan(const std::string& t)	 { //! Add text span element.
	 //data_.emplace_back(new tspan_element(t, style_));
	 //return *(static_cast<tspan_element*>(data_.back().get()));	 }

  //tspan_element& tspan(const std::string& t, const text_style& style)	 { //! Add text span element (with specified text style).
	 //data_.emplace_back(new tspan_element(t, style));
	 //return *(static_cast<tspan_element*>(data_.back().get()));	 }


  std::string text()	 { //! \return text string of a text_element.
	 std::stringstream os;
	 generate_text(os);
	 return os.str();	 }

  size_t data_size()	 { 	 return data_.size();	 }//! \return number of text sub elements


  void write(std::ostream& os)	 { //! Output text_element, style & attributes to stream.
	 // Changed to new convention on spaces:
	 // NO trailing space, but *start* each item with a space.

	 if(!data_.size()) return;

	 os << "\n<text x=\"" << x_ << "\" y=\"" << y_ << "\"";
	 std::string anchor;
	 switch(align_)
	 {
	  case left_align: 	 anchor = ""; 	 break;	// anchor = "start"; // This is the initial == default. so should be possible to reduce file size this by:
	  case right_align: 	 anchor = "end"; 	 break;
	  case center_align: 	 anchor = "middle"; 	 break;
	  default: 	 anchor = ""; 	 break;
	 }
	 if(anchor.size())   os << " text-anchor=\"" << anchor << "\"";
	 if(rotate_)         os << " transform = \"rotate(" << rotate_ << " " << x_ << " " << y_ << " )\"";

	 if (style_.font_size_)          os << " font-size=\"" << style_.font_size_ << "\"";

	 if (style_.font_family_.size())  os << " font-family=\"" << style_.font_family_ << "\"";  // Example: Arial.

	 if (style_.font_style_.size())  os << " font-style=\"" << style_.font_style_ << "\"";  // Example: italic.
	 
	 if (style_.weight_.size())      os << " font-weight=\"" << style_.weight_ << "\"";	  // Example: bold.
	 
	 if (style_.stretch_.size())     os << " font-stretch=\"" << style_.stretch_ << "\"";
	 
	 if (style_.decoration_.size())  os << " font-decoration=\"" << style_.decoration_ << "\"";
	 //if (style_.baseline_.size())  os << " dominant-baseline=\"" << style_.baseline_ << "\"";
	 os << '>' ;
	 generate_text(os);
	 os << "</text>";
	 // Example:
  } // void write(std::ostream& os)
}; // class text_element_



struct path_point	//! Base class for m_path, z_path, q_path, h_path, v_path, c_path, s_path.
{	/*! \struct svg::path_point
	\details Paths represent the outline of a shape which can be
	filled, stroked, used as a clipping path, or any combination of the three.
	*/
	bool relative; //!< If true relative else absolute.

	virtual void write(std::ostream& rhs) = 0; //!< write functions output SVG commands like "M1.2, 3.4",
	virtual ~path_point() 	 {  } //! Destructor.

	path_point(bool relative) : relative(relative) 	 {  }
}; // struct path_point



  struct m_path: public path_point
  { /*! \struct svg::m_path
      \brief move to coordinates (x, y)
     \details See 8.3.2 The "moveto" commands.
     */
    double x; //!< End of move SVG X coordinate.
    double y; //!< End of move SVG Y coordinate.

    void write(std::ostream& o_str)    { //! write moveto X and Y coordinates to stream, for example: "M52.8571,180 "
		o_str << ( relative ? "m" : "M" ) << x << "," << y << " ";    }// separator changed to comma for clarity when reading XML source.
 

    m_path(double x, double y, bool relative = false) //! Construct a move to
      : path_point(relative),  x(x), y(y)
    {   }
  }; // struct m_path


  struct z_path: public path_point
  { /*! \struct svg::z_path
       \brief Close current path.
       \details
       http://www.w3.org/TR/SVG/paths.html#PathElement
    */
    void write(std::ostream& o_str)    {      o_str << "Z";    } //! Write close current path SVG command.

    z_path() : path_point(false)   {    } //! Constructor defines all member variables.

  }; // struct z_path


  // 8.3.4 The "lineto" commands L, H & V.
  struct l_path: public path_point
  {  /*! \struct svg::l_path
      \brief Draw a line from the current point to the given (x,y) coordinate
       which becomes the new current point.
    */
    double x; //!< End of line SVG X coordinate.
    double y; //!< End of line SVG Y coordinate.
    void write(std::ostream& o_str)
    { //! Write line to SVG command.
        o_str << (relative ? "l" : "L" ) << x << "," << y << " ";
    }

    l_path(double x, double y, bool relative = false)
      : path_point(relative), x(x), y(y)
    { //! Constructor defines all member variables.
    }
  }; // struct l_path


  struct s_path : public path_point
  { /*! \struct svg::s_path
      \brief Draws a cubic Bezier curve from the current point to (x,y).
      \details see also t_path for a quadratic Bezier curve.
    */
    double x1; //!< Cubic Bezier curve control X coordinate.
    double y1; //!< Cubic Bezier curve control Y coordinate.
    double x; //!< Cubic Bezier curve end X coordinate.
    double y; //!< Cubic Bezier curve end Y coordinate.
    void write(std::ostream& o_str)
    {  o_str <<(relative ? "s" : "S" ) << x1 << "," << y1 << " " << x << "," << y << " "; } //! Write SVG Cubic Bezier curve command.

    s_path(double x1, double y1, double x, double y, bool relative = false)
      : path_point(relative), x1(x1), y1(y1), x(x), y(y)   {  } //! Constructor Cubic Bezier curve.
  
  }; // struct s_path

//#define svg_EXTRA_SHAPES
//#include <svg_elements_extra.hpp>


  class path_element: public svg_element
  {  /*! \class svg::path_element
     \brief Path element holds places on a path used by move, line ...
     \details
     http://www.w3.org/TR/SVG/paths.html#PathElement
     8.3.1 General information about path data.
     A path is defined by including a 'path'  element
     which contains a d="(path data)"  attribute,
     where the d attribute contains the moveto, line, curve
     (both cubic and quadratic Beziers), arc and closepath instructions.
     */
  public: //temporary for experimental gil
//  private:
    std::vector<std::unique_ptr<path_point> > path; //!< All the (x, y) coordinate pairs,
    //!< filled by calls of m, M, l , L... that emplace_back.
    // bool fill; now inherited from parent svg class.

    //path_element(const path_element& rhs)
    //{ //! Copy constructor.
      //path = (const_cast<path_element&>(rhs)).path.release();
    //}

    path_element(const svg_style& style_info)
      :
      svg_element(style_info)
    { //! Construct empty path element.
    }

    path_element()
    { //! Construct an empty path element.
      // fill now got from the parent svg fill color.
    }

    path_element& fill_on(bool on_)
    { //! Set area fill, on or off.
      style_info_.fill_on(on_);	 return *this;
    }

    bool fill_on()
    { //! \return area fill, on or off.
      return style_info_.fill_on();
    }
    // Note 1: return of path_element& permits chaining calls like
    // my_path.M(3, 3).l(150, 150).l(200, 200)...;

    // Note 2: By convention:
    // lower case (like m call m_path(x, y, true) for relative
    // but upper case (like M) calls m_path(x, y, false) for absolute.

    path_element& m(double x, double y)
    { //! Move relative by x and y.
      path.emplace_back(new m_path(x, y, true));	 return *this;
    }

    path_element& M(double x, double y)
    { //! Move to absolute x and y.
      path.emplace_back(new m_path(x, y, false));	 return *this;
    }

    path_element& l(double x, double y)
    { //! Line to (relative).
      path.emplace_back(new l_path(x, y, true));	 return *this;
    }

    path_element& L(double x, double y)
    { //! Line to (absolute).
      path.emplace_back(new l_path(x, y, false));	 return *this;
    }
    path_element& z()
    { //! Path end. Note lower case z, see path_element& Z() below.
      path.emplace_back(new z_path());	 return *this;
    }

    path_element& Z()
    { //! Path end. Note Upper case Z also provided for compatibility with
      //! http://www.w3.org/TR/SVG/paths.html#PathDataClosePathCommand 8.3.3 which allows either case.
      path.emplace_back(new z_path());	 return *this;
    }

    path_element& s(double x1, double y1, double x, double y)
    { //! Draws a cubic Bezier curve from the current point to (x,y) (relative).
      path.emplace_back(new s_path(x1, y1, x, y, true));	 return *this;
    }

    path_element& S(double x1, double y1, double x, double y)
    { //! Draws a cubic Bezier curve from the current point to (x,y) (absolute).
      path.emplace_back(new s_path(x1, y1, x, y, false));	 return *this;
    }

    void write(std::ostream& o_str)
    { //! Write SVG path command.
      //! Example: \verbatim <path d="M5,175 L5,195 M148.571,175" /> \endverbatim
      if (path.begin() != path.end() )
      { // Is some path info (trying to avoid useless <path d=""/>"
        // TODO or would this omit useful style & attributes?
        o_str << "<path d=\"";
        for(auto i = path.begin(); i != path.end(); ++i)
          (*i)->write(o_str); // M1,2
        o_str << "\"";

        write_attributes(o_str); // id & clip_path
        style_info_.write_out(o_str); // fill, stroke, width...

        // line above should write fill = "none" that seems to be needed for reasons unclear.  Even when g_element does not specify a fill, it seems to be interpreted as black fill.
        if(!fill_on()) 	  o_str << " fill=\"none\"";
        o_str<<"/>"; // closing to match <path d=
      }

    } // void write(std::ostream& o_str)

	//#define path_element_EXTRA_SHAPES
	//#include <svg_elements_extra.hpp>

  }; // class path_element



  class g_element: public svg_element
  { /*! \class svg::g_element
      \brief g_element (group element) is the node element of our document tree.
      'g' element is a container element for grouping together:  \verbatim  <g> ... </g>  \endverbatim
      \details g_element ('g' element is a container element
      for grouping together related graphics elements).\n
      See http://www.w3.org/TR/SVG/struct.html#NewDocument 5.2.1 Overview.


      'g' element is a container element for grouping together \verbatim <g /> </g>\endverbatim
      related graphics elements, for example, an image background rectangle with border and fill:
      \verbatim <g id="background" fill="rgb(255,255,255)"><rect width="500" height="350"/></g>
      \endverbatim
   */
  public: //temporary for experimental gil

//  private:
    std::vector<std::unique_ptr<svg_element> > children; /*!< Children of this group element node,
      containg graphics elements like text, circle line, polyline...
    */
    std::string clip_name;  //!< Name of clip path.
    bool clip_on; //!< true if to clip anything outside the clip path.
  public:

    g_element() : clip_on(false)   { }//! Construct g_element (with no clipping).

    g_element(const svg_style& styl) :svg_element(styl), clip_on(false)  { }//! Construct g_element (with no clipping).

    svg_element& operator[](unsigned int i)  { return *children[i];  } //! \return child svg_element node.

    size_t size()    {       return children.size();    }//! \return Number of child nodes.

    void write(std::ostream& os)
    { /*! Output all children of a group element.
         \verbatim
         Example:
         <g fill="rgb(255,255,255)" id="background"><rect x="0" y="0" width="500" height="350"/></g>
         \endverbatim
      */

      if (children.size() > 0)        //   Avoid useless output like: <g id="legendBox"></g>
      {
        os << "<g"; // Do NOT need space if convention is to start following item with space.
        write_attributes(os); // id="background" (or clip_path)
        style_info_.write_out(os); // stroke="rgb(0,0,0)" fill= "rgb(255,0,0)" ...
        os << ">" ;
        for(unsigned int i = 0; i < children.size(); ++i)
        {
          children[i]->write(os);
        }
        os << "</g>" << std::endl;
      }
    } // void write(std::ostream& rhs)

    g_element& g(int i)
    { //! i is index of children nodes.
      return *(static_cast<g_element*>(children[i].get())); // \return reference to child node g_element.
    }

    g_element& add_g_element(svg_style styl = svg_style())
    { //! Add a new group element.
      //! \return A reference to the new child node just created.
      children.emplace_back(new g_element(styl));
      return *(static_cast<g_element*>(children.back().get()));
    }

    line_element& line(double x1, double y1, double x2, double y2)
    { //! Add a new line element.
      //! \return A reference to the new child node just created.
      children.emplace_back(new line_element(x1, y1, x2, y2));
      return *(static_cast<line_element*>(children.back().get()));
    }

    rect_element& rect(double x1, double y1, double dx, double dy)
    { //! Add a new rect element.
      //! \return A reference to the new child node just created.
      children.emplace_back(new rect_element(x1, y1, dx, dy));
      return *(static_cast<rect_element*>(children.back().get()));
    }

    circle_element& circle(double x, double y, unsigned int radius = 5)
    { //! Add a new circle element.
      //! \return A reference to the new child node just created.
      children.emplace_back(new circle_element(x, y, radius));
      return *(static_cast<circle_element*>(children.back().get()));
    }

    ellipse_element& ellipse(double rx, double ry, double cx, double cy)
    { //! Add a new ellipse element.
      //! \return A reference to the new child node just created.
      children.emplace_back(new ellipse_element(rx, ry, cx, cy));
      return *(static_cast<ellipse_element*>(children.back().get()));
    }

    text_element& text(double x = 0., double y = 0.,
    const std::string& text = "",
    const text_style& style = no_style, // Use svg implementation's defaults.
    const align_style& align = left_align,
    const rotate_style& rotate = horizontal)
    { //! Add a new text element.
      //! \return A reference to the new child node just created.
      children.emplace_back(new text_element(x, y, text, style, align, rotate));
      return *(static_cast<text_element*>(children.back().get()));
    }

	//#define g_element_EXTRA_SHAPES
	//#include <svg_elements_extra.hpp>

    path_element& path()
    { //! Add a new path element.
      //! \return A reference to the new path just created.
      children.emplace_back(new path_element()); // Empty path.
      return *(static_cast<path_element*>(children.back().get()));
    }

    void push_back(svg_element* g)
    { //! Add a new child node g_element.
      children.emplace_back(g);
    }

    void clear()    {      children.clear();    } //! Remove all the child nodes.
  }; // class g_element

} // namespace svg

#endif // BOOST_SVG_TAG_HPP
