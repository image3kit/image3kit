/*!
  \file svg_color.hpp
  \brief SVG standard names of colors, and functions to create and output colors.
  \date 9 Feb 2009
  \author Jacob Voytko & Paul A. Bristow
*/
// Copyright Jacob Voytko 2007
// Copyright Paul A. Bristow 2007, 2009
// Copyright Ali Q. Raeini 2019
//
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_SVG_SVG_COLOR_HPP
#define BOOST_SVG_SVG_COLOR_HPP


#include <ostream>

//namespace boost {
							namespace svg _begins_

  /*!
    \brief Colors that have SVG standard special names.
    \details  The reason that the underscore separator convention does not match
    the normal Boost format is that these names that are specified by the SVG standard.
    http://www.w3.org/TR/SVG/types.html#ColorKeywords
    color "tan" is also renamed to "tanned" to avoid clash with global function name tan in math.h.
  */
 enum svg_color_constant
 { //! \enum svg_color_constant SVG standard names for some colors.
    //! See http://www.w3.org/TR/SVG/types.html#ColorKeywords
    blank=0,
    none=1, transparent = 1,
    black, darkblue, darkred, darkgreen, darkmagenta, darkcyan, darkgoldenrod, aliceblue, brown, lightgreen1, aquamarine, azure, beige,
    bisque, blanchedalmond, blue, blueviolet,
    burlywood, cadetblue, chartreuse, chocolate, coral,
    cornflowerblue, cornsilk, crimson, cyan,
    darkgray, darkgrey, darkkhaki,
    darkolivegreen, darkorange, darkorchid,
    darksalmon, darkseagreen, darkslateblue, darkslategray,
    darkslategrey, darkturquoise, darkviolet,
    deepskyblue, deeppink, dimgray, dimgrey, dodgerblue, firebrick,
    floralwhite, forestgreen, fuchsia, gainsboro, ghostwhite, gold,
    goldenrod, gray, grey, green, greenyellow, honeydew, hotpink,
    indianred, indigo, ivory, khaki, lavender, lavenderblush,
    lawngreen, lemonchiffon, lightblue, lightcoral, lightcyan,
    lightgoldenrodyellow, lightgray, lightgreen, lightgrey,
    lightpink, lightsalmon, lightseagreen, lightskyblue,
    lightslategray, lightslategrey, lightsteelblue, lightyellow,
    lime, limegreen, linen, magenta, maroon, mediumaquamarine,
    mediumblue, mediumorchid, mediumpurple, mediumseagreen,
    mediumslateblue, mediumspringgreen, mediumturquoise,
    mediumvioletred, midnightblue, mintcream, mistyrose, moccasin,
    navajowhite, navy, oldlace, olive, olivedrab, orange,
    orangered, orchid, palegoldenrod, palegreen, paleturquoise,
    palevioletred, papayawhip, peachpuff, peru, pink, plum,
    powderblue, purple, red, rosybrown, royalblue, saddlebrown,
    salmon, sandybrown, seagreen, seashell, sienna, silver,
    skyblue, slateblue, slategray, slategrey, snow, springgreen,
    steelblue, tanned,
    // tan, // Note that tan would clash with geometric tan in math.h!
    teal, thistle, tomato, turquoise, violet,
    wheat, whitesmoke, antiquewhite, yellow, yellowgreen,
    white,
    dark
  }; // enum svg_color_constant

  void constant_to_rgb(svg_color_constant c, unsigned char& r, unsigned char& g, unsigned char& b);

  /*! \brief SVG standard colors, see also enum svg_color_constant
      \details svg_color is the struct that contains information about RGB colors.
      For the constructor, the SVG standard specifies that numbers
      outside the normal rgb range are to be accepted,
      but are constrained to acceptable range of integer values [0, 255].
  */
 class svg_color
 {
    //friend bool operator== (const svg_color& lhs, const svg_color& rhs);
    //friend bool operator!= (const svg_color& lhs, const svg_color& rhs);
    friend bool is_blank(const svg_color& col);

  public: //temporary for experimental gil
	svg_color() {};
//  private:
    unsigned char r_; //!< red unsigned char provides range [0 to 255].
    unsigned char g_; //!< green unsigned char provides range [0 to 255].
    unsigned char b_; //!< blue unsigned char provides range [0 to 255].
    unsigned char o_; //!< opacity unsigned char provides range [0 to 255]. 0 is reserved, it means not to be drawn
    //bool is_blank_; //!< true means "Not to be displayed" a 'pseudo-color'. If is_blank_ == true should write output to SVG XML file as "none".

  public:

    svg_color(int red, int green, int blue, int opacity=255)
    {  //! \brief Construct an SVG color from RGB values.
       //!\details Constrain rgb to [0 .. 255].
      //!  Default is to construct a 'pseudo-color' blank.
      opacity = ( opacity < 0 ) ? 0 : opacity;
      red = ( red < 0 ) ? 0 : red;
      green = ( green < 0 ) ? 0 : green;
      blue = ( blue < 0 ) ? 0 : blue;
      o_ = (unsigned char)(( opacity > 255 ) ? 255 : opacity);
      r_ = (unsigned char)(( red > 255 ) ? 255 : red);
      g_ = (unsigned char)(( green > 255 ) ? 255 : green);
      b_ = (unsigned char)(( blue > 255 ) ? 255 : blue);
    } // svg_color(int red, int green, int blue)


    svg_color(svg_color_constant col)
    { //! Set a color (including blank) using the SVG 'standard' colors defined in enum svg::svg_color_constant
      if (col == blank)
      { // NotAColor.
        r_ = 255; // Safer to assign *some* value to rgb: zero, or 255 or something,
        g_ = 255; // rather than leaving them random.
        b_ = 255; // Default 'blank' color here is white.
        o_ = 0;
      }
      else
      { // Proper color.
        o_ = 255;
        constant_to_rgb(col, r_, g_, b_);
      }
    }


	bool is_visible() const 	{ 	return o_; 	} //! \return true if color is blank.


	svg_color operator*(float frac) const 	{ 	return svg_color(r_,g_,b_,o_*frac); 	}    // obselete
	svg_color operator/(float frac) const 	{ 	return svg_color(r_/frac,g_/frac,b_/frac,o_); 	} // obselete
	svg_color rgbX(float frac) const 	{ 	return svg_color(r_*frac,g_*frac,b_*frac,o_); 	}
	svg_color opaX(float frac) const 	{ 	return svg_color(r_,g_,b_,o_*frac); 	}

	//unsigned int red() const 	{ 	return r_; 	}   //! \return red component of color [0, 255]
	//unsigned int green() const 	{ 	return g_; 	} //! \return green component of color [0, 255]
	//unsigned int blue() const 	{ 	return b_; 	} //! return blue component of color [0, 255]

	void write(std::ostream& os) const
	{ //! Write to ostream a color in svg format.
		//! \details Usage: my_color.write(cout); Outputs: rgb(127,255,212)
		if(o_) os << "rgb(" << (unsigned int)r_ << "," << (unsigned int) g_ << "," << (unsigned int)b_ << ")" ;
		else   os << "none";
	}
 }; // class svg_color


 inline bool is_visible(const svg_color& col)
 { //! \return true if color is blank.
	return col.o_;
 }

 inline std::ostream& operator<<(std::ostream& os, const svg_color& color)
 {	//!\brief Output color to stream as RGB. or none, compatible with old (version 1) svg specifications
 	color.write(os);
 	return os;
 }


  //! SVG standard colors, \see svg_color_constant
  const static svg_color color_array[] =
  {
    svg_color(255,255,255,0),          // blank - "Not to be displayed" pseudo-color.
    svg_color(255,255,255,1),          // "none" - fully  transparent color.
    svg_color(0  , 0  , 0  ), // black
    svg_color(60 , 120 , 255), // darkblue
    svg_color(245, 110, 110  ), // darkred
    svg_color(50 , 180, 0  ), // darkgreen
    svg_color(120, 0  , 210), // darkmagenta
    svg_color(0  , 200, 240), // darkcyan
    svg_color(204, 174, 22 ), // darkgoldenrod
    svg_color(200, 208, 215), // aliceblue
    svg_color(145, 72 , 42 ), // brown
    svg_color(200, 240, 0  ), // lightgreen1
    svg_color(127, 255, 212), // aquamarine [4]
    svg_color(220, 235, 235), // azure
    svg_color(245, 245, 220), // beige
    svg_color(255, 228, 196), // bisque
    svg_color(255, 235, 205), // blanchedalmond
    svg_color(0  , 0  , 255), // blue
    svg_color(138, 43 , 226), // blueviolet
    svg_color(222, 184, 135), // burlywood
    svg_color(95 , 158, 160), // cadetblue
    svg_color(127, 255, 0  ), // chartreuse
    svg_color(210, 105, 30 ), // chocolate
    svg_color(255, 127, 80 ), // coral
    svg_color(100, 149, 237), // cornflowerblue
    svg_color(255, 248, 220), // cornsilk
    svg_color(220, 20 , 60 ), // crimson
    svg_color(0  , 255, 255), // cyan
    svg_color(169, 169, 169), // darkgray
    svg_color(169, 169, 169), // darkgrey
    svg_color(189, 183, 107), // darkkhaki
    svg_color(85 , 107, 47 ), // darkolivegreen
    svg_color(255, 140, 0  ), // darkorange
    svg_color(153, 50 , 204), // darkorchid
    svg_color(233, 150, 122), // darksalmon
    svg_color(143, 188, 143), // darkseagreen //#008BFF #FF0000  #00FF00
    svg_color(72 , 61 , 139), // darkslateblue
    svg_color(47 , 79 , 79 ), // darkslategray
    svg_color(47 , 79 , 79 ), // darkslategrey
    svg_color(0  , 206, 209), // darkturquoise
    svg_color(148, 0  , 211), // darkviolet
    svg_color(0  , 191, 255), // deepskyblue
    svg_color(255, 20 , 147), // deeppink
    svg_color(105, 105, 105), // dimgray
    svg_color(105, 105, 105), // dimgrey
    svg_color(30 , 144, 255), // dodgerblue
    svg_color(178, 34 , 34 ), // firebrick
    svg_color(255, 250, 240), // floralwhite
    svg_color(34 , 139, 34 ), // forestgreen
    svg_color(255, 0  , 255), // fuchsia
    svg_color(220, 220, 220), // gainsboro
    svg_color(248, 248, 255), // ghostwhite
    svg_color(255, 215, 0  ), // gold
    svg_color(218, 165, 32 ), // goldenrod
    svg_color(128, 128, 128), // gray
    svg_color(128, 128, 128), // grey
    svg_color(0  , 255, 0  ), // green
    svg_color(173, 255, 47 ), // greenyellow
    svg_color(240, 255, 240), // honeydew
    svg_color(255, 105, 180), // hotpink
    svg_color(205, 92 , 92 ), // indianred
    svg_color(75 , 0  , 130), // indigo
    svg_color(255, 255, 240), // ivory
    svg_color(240, 230, 140), // khaki
    svg_color(230, 230, 250), // lavender
    svg_color(255, 240, 245), // lavenderblush
    svg_color(124, 252, 0  ), // lawngreen
    svg_color(255, 250, 205), // lemonchiffon
    svg_color(173, 216, 230), // lightblue
    svg_color(240, 128, 128), // lightcoral
    svg_color(224, 255, 255), // lightcyan
    svg_color(250, 250, 210), // lightgoldenrodyellow
    svg_color(211, 211, 211), // lightgray
    svg_color(144, 238, 144), // lightgreen
    svg_color(211, 211, 211), // lightgrey
    svg_color(255, 182, 193), // lightpink
    svg_color(255, 160, 122), // lightsalmon
    svg_color(32 , 178, 170), // lightseagreen
    svg_color(135, 206, 250), // lightskyblue
    svg_color(119, 136, 153), // lightslategray
    svg_color(119, 136, 153), // lightslategrey
    svg_color(176, 196, 222), // lightsteelblue
    svg_color(255, 255, 224), // lightyellow
    svg_color(0  , 255, 0  ), // lime
    svg_color(50 , 205, 50 ), // limegreen
    svg_color(250, 240, 230), // linen
    svg_color(255, 0  , 255), // magenta
    svg_color(128, 0  , 0  ), // maroon
    svg_color(102, 205, 170), // mediumaquamarine
    svg_color(0  , 0  , 205), // mediumblue
    svg_color(186, 85 , 211), // mediumorchid
    svg_color(147, 112, 219), // mediumpurple
    svg_color(60 , 179, 113), // mediumseagreen
    svg_color(123, 104, 238), // mediumslateblue
    svg_color(0  , 250, 154), // mediumspringgreen
    svg_color(72 , 209, 204), // mediumturquoise
    svg_color(199, 21 , 133), // mediumvioletred
    svg_color(25 , 25 , 112), // midnightblue
    svg_color(245, 255, 250), // mintcream
    svg_color(255, 228, 225), // mistyrose
    svg_color(255, 228, 181), // moccasin
    svg_color(255, 222, 173), // navajowhite
    svg_color(0  , 0  , 128), // navy
    svg_color(253, 245, 230), // oldlace
    svg_color(128, 128, 0  ), // olive
    svg_color(107, 142, 35 ), // olivedrab
    svg_color(255, 165, 0  ), // orange
    svg_color(255, 69 , 0  ), // orangered
    svg_color(218, 112, 214), // orchid
    svg_color(238, 232, 170), // palegoldenrod
    svg_color(152, 251, 152), // palegreen
    svg_color(175, 238, 238), // paleturquose
    svg_color(219, 112, 147), // palevioletred
    svg_color(255, 239, 213), // papayawhip
    svg_color(255, 218, 185), // peachpuff
    svg_color(205, 133, 63 ), // peru
    svg_color(255, 192, 203), // pink
    svg_color(221, 160, 221), // plum
    svg_color(176, 224, 230), // powderblue
    svg_color(128, 0  , 128), // purple
    svg_color(255, 0  , 0  ), // red
    svg_color(188, 143, 143), // rosybrown
    svg_color(65 , 105, 225), // royalblue
    svg_color(139, 69 , 19 ), // saddlebrown
    svg_color(250, 128, 114), // salmon
    svg_color(244, 164, 96 ), // sandybrown
    svg_color(46 , 139, 87 ), // seagreen
    svg_color(255, 245, 238), // seashell
    svg_color(160, 82 , 45 ), // sienna
    svg_color(192, 192, 192), // silver
    svg_color(135, 206, 235), // skyblue
    svg_color(106, 90 , 205), // slateblue
    svg_color(112, 128, 144), // slategray
    svg_color(112, 128, 144), // slategrey
    svg_color(255, 250, 250), // snow
    svg_color(0  , 255, 127), // springgreen
    svg_color(70 , 130, 180), // steelblue
    svg_color(210, 180, 140), // tanned
    svg_color(0  , 128, 128), // teal
    svg_color(216, 191, 216), // thistle
    svg_color(255, 99 , 71 ), // tomato
    svg_color(64 , 224, 208), // turquoise
    svg_color(238, 130, 238), // violet
    svg_color(245, 222, 179), // wheat
    svg_color(245, 245, 245), // whitesmoke
    svg_color(250, 235, 215), // antiquewhite
    svg_color(255, 255, 0  ), // yellow
    svg_color(154, 205, 50 ), // yellowgreen
    svg_color(255, 255, 255), // white
    svg_color(0, 0, 0,127), // dark
  }; // svg_color color_array[]


 inline void constant_to_rgb(svg_color_constant c, unsigned char& r, unsigned char& g, unsigned char& b)
 { /*! Convert a named SVG standard color, see enum svg::svg_color_constant
 		to update three variables (r, g, b) holding red, green and blue values.
 		Asserts that c NOT the blank color.
 	*/
 	dAsrt(c != blank);
 	const svg_color color = color_array[c];
 	r = color.r_;
 	g = color.g_;
 	b = color.b_;
 }



 inline void rgbtohsv(int fR, int fG, int fB, int& fH, int& fS, int& fV) {
 	using std::min;
 	using std::max;
 	fV = max(max(fR, fG), fB);
 	fS = max(fV - min(min(fR, fG), fB), 1);

 	if(fS > 0) {
 		if(fV == fR)        {      fH =  ((40 * (fG - fB) / fS) );//% 240
 		} else if(fV == fG) {      fH =  ((40 * (fB - fR) / fS) + 80);
 		} else              {      fH =  ((40 * (fR - fG) / fS) + 160);    }

 		if(fH < 0) {    fH = (240 + fH);  }

 	} else {    fH=0;  }
 }


 inline void hsvtorgb(unsigned char& fR, unsigned char& fG, unsigned char& fB,  const int& fH, const int& fS, const int& fV) {
 	//int fC = std::max(fS,1);*fV // Chroma
 	int fX =  fS - abs((fH*fS/40)%(2*fS+1) - fS);
 	int fM = fV - fS;
 	switch (fH/40)
 	{
 		case 0:  fR = fS;    fG = fX;   fB=0; break;
 		case 1:  fR = fX;    fG = fS;   fB=0; break;//yellow
 		case 2:  fR=0;     fG = fS;   fB = fX; break;
 		case 3:  fR=0;     fG = fX;   fB = fS; break;//cyan
 		case 4:  fR = fX;    fG=0;    fB = fS; break;
 		case 5:  fR = fS;    fG=0;    fB = fX; break;//magenta
 		default: fR = fS;    fG=0;    fB = fX;  //case 6/0
 	}
 	fR = std::min(fR+fM,255);  fG = std::min(fG+fM,255);  fB = std::min(fB+fM,255);
 }

 inline svg_color rgbtween(const float& f, svg_color c0=svg_color(244,246,250), svg_color c1=svg_color(255,180,180))
 {	//! HSV interpolation between RGB colours,
 	//! tries to keep value to vary linearly too
 	auto mag = [](const svg_color& rgb) {return int(rgb.r_)+rgb.g_+rgb.b_;}; //syncXdsdsfkghhfd * 4 *5 *3
 	//auto mag = [](const svg_color& rgb) -> int {return 4*int(rgb.r_)*int(rgb.r_)+ 5*int(rgb.g_)*int(rgb.g_)+ 3*int(rgb.b_)*int(rgb.b_);};
 	int h0,s0,v0,h1,s1,v1;
 	rgbtohsv(c0.r_,c0.g_,c0.b_,h0,s0,v0);
 	rgbtohsv(c1.r_,c1.g_,c1.b_,h1,s1,v1);
 	int mg0 = (mag(c0)*2+v0*1);// optional, to keep magnitude constantish
 	int mg1 = (mag(c1)*2+v1*1);// optional, to keep magnitude constantish

 	h1 = int(240+h0+f*(h1-h0))%240;
 	s1 = (s0+f*(s1-s0))  ;
 	v1 = v0+f*(v1-v0);
 	mg1 = mg0+f*(mg1-mg0)+2;

 	svg_color c(0,0,0, c0.o_+f*(c1.o_-c0.o_));
 	hsvtorgb(c.r_,c.g_,c.b_,h1,s1,v1);

 	{                      // optional, to keep magnitude constantish
 		mg0 = (mag(c)*2+v1*1)+2;
 		mg1=(mg1-mg0)/8; //syncXdsdsfkghhfd
 		c.r_=std::min(std::max(0,mg1+c.r_),255);
 		c.g_=std::min(std::max(0,mg1+c.g_),255);
 		c.b_=std::min(std::max(0,mg1+c.b_),255);
 	}

 	return c;
 }
 inline svg_color rgbtween_0w(const float& f, svg_color c0, svg_color c1)
 {//! alternative colour interpolation, not used by internally
 	svg_color c01=rgbtween(f,c0,c1);
 	if(c01.g_>c01.b_ && c01.g_>c01.r_)
 	{ int del=c01.g_-std::max(c01.r_,c01.b_);
 		 c01.r_+=3*del/4;   c01.b_+=3*del/4;  }
 	//if(f<0.5f) return rgbtween(2.0f*f,c0,c01);
 	//else
 	return c01;
 }
							_end_of_(namespace svg)

#endif // BOOST_SVG_SVG_COLOR_HPP
