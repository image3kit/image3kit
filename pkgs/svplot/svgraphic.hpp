
// Copyright Jacob Voytko 2007
// Copyright Paul A Bristow 2007, 2009
// Copyright Ali Q. Raeini 2019

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)


/*!
	\file svg.hpp
	\brief Scalable Vector Graphic (SVG) format elements.
	\details Provides classes and methods to create the basic SVG graph elements.
	Graph elements, point, path, line, circle, rect and polygon, text are used by the 2D and Boxplot functions,
	but could also be used for generating other graphics in SVG format.

	\author Jacob Voytko & Paul A. Bristow & Ali Q Raeini
 */

/*!
	\mainpage Scalable Vector Graphic (SVG) Plot Package

	\section intro_sec Introduction

	Humans have a fantastic capacity for visual understanding, and merely looking
	at data organized in one, two, or three dimensions allows us to see relations
	not otherwise visible in a list of numbers. Computers, however, deal with
	information numerically, and C++ and the
	\htmlonly
		<a href="http://en.wikipedia.org/wiki/C%2B%2B#Standard_library"> Standard Template Library (STL) </a>
	\endhtmlonly
	do not currently offer a way to bridge the gap.

	This library helps the user to easily plot data stored in STL containers.

	This project is focused on using STL containers in order to graph data on
	one- or two-dimensional (2D) plots.

	The plot is currently be written to a
	\htmlonly
		<a href="http://en.wikipedia.org/wiki/Scalable_Vector_Graphics"> Scalable Vector Graphic image </a>
	\endhtmlonly
	\htmlonly
		<a href="http://www.w3.org/TR/SVG11/"> Scalable Vector Graphics (SVG) </a>
	\endhtmlonly
	\htmlonly
		<a href="http://www.w3.org/TR/REC-xml/#sec-comments"> XML specification</a>
	\endhtmlonly

	and file format for describing two-dimensional vector graphics.

	SVG files (.svg) can be viewed with any modern internet browser:

	- Mozilla Firefox
	- Google Chrome
	- Adobe Illustrator
	- Opera,
	- Microsoft Internet Explorer.
	- \htmlonly
		  <a href="http://www.inkscape.org"> Inkscape</a>
	  \endhtmlonly, which has a lot of editing features, and conversion to wide variety of other graphics formats
		like pdf and png.

	- And by Many other graphics programs, for example
	 \htmlonly
	   <a href="http://svg.software.informer.com/software/"> Most popular SVG software</a>
	 \endhtmlonly

	The goals of the project are:

	- To let users produce simple plots with minimal intervention by using sane defaults.
	- To allow users to easily customize plots but allow very fine-grained control of appearance.
	- To allow the user to talk to the plot classes using coordinate units rather than pixels or other arbitrary measures.
	- To produce and transfer plots quickly enough to be useful in real-time.
	- To represent uncertainty visually and numerically by showing only significant digits,
	- To create the backbone of a `svg` class that can be extended to fully support the standard.
	- Compliance with the  \htmlonly <a href="http://www.w3.org/TR/SVG11/"> W3C Scalable vector Graph standard.</a> \endhtmlonly
	- Validation using W3C Markup Validation Service.
	- Copyright and license conditions built into the SVG file.

	\section why_svg Why SVG format?

	SVG provides very high quality images that display well on small screens like handheld devices,
	conventional monitors of all sizes,
	and on big screens (for example, projected on lecture theatre screens),
	and especially when printed, but files are tiny
	(compared to images produced by spreadsheets like Microsoft Excel).

	SVG files are microscopic when compressed using a zip program
	to convert to types like .zip, or the specific compressed SVG file type (.svgz).

	SVG files can be added to other documents like Adobe PDF, and display well, without bloat.


	Many types of plots can be produced:


	- 2D plots allow X versus Y variables to be displayed.  You can show values,
	  and optionally their uncertainty (the deprecated term is 'error bars').

	- Boxplots can be produced, optionally with values and outliers.

	- Bar charts can also be produced.

	This page generated from file svg.hpp.

*/
#ifndef BOOST_SVGRAPHIC_HPP
#define BOOST_SVGRAPHIC_HPP


#define _begins_     {
#define _end_of_(sec) }

#include <string>
#include <ostream>
#include <fstream>
#include <vector>
#include <string>

#include "svg_elements.hpp" // element class definitions.
#include "svplot_elems.hpp"
//#include "svplot.hpp" // Could be used to check declarations and definitions match correctly.

// SVG stands for Scalable Vector Graphics,
// an XML grammar for stylable graphics, usable as an XML namespace.

// Gzip compression - can give files that are 1/10th size of gif or jpeg. TODO: implement with zlib maybe
// Use default values whenever possible rather than defining all attributes and properties explicitly.
// Take advantage of the path data compaction facilities so use relative coordinates.
// Use h and v for horizontal and vertical lines.
// Use s or t for cubic and quadratic Bezier segments whenever possible.
// Eliminate extraneous white space and separators.
// Use symbols if the same graphic appears multiple times in the document.
// Use CSS property inheritance and selectors to consolidate commonly used properties into named styles
// or to assign the properties to a parent group element.
// Use filter effects to help construct graphics via client-side graphics operations.

//namespace boost {
							namespace svg _begins_
   //! \namespace boost \brief www.Boost.org. REMOVED
	//! \namespace svg \brief Scalable Vector Graph plot functions, classes and data.


//! base class for svg_2d_plot etc
class svchart : public g_element
{
 public:
	svchart(double Dx ,double Dy ,double X0=0.0 ,double Y0=0.0)
	:  X0_(X0),Y0_(Y0),Dx_(Dx),Dy_(Dy), updated_(false), clip_rect(X0,Y0,Dx,Dy)
	{
		static int plotIndex=0;	 pl_window_clip_=("pl_window"+std::to_string(++plotIndex));
	}
	virtual ~svchart() {  }
	double X0_,Y0_,Dx_,Dy_; // plot X0, Y0, and width and height
	bool updated_;  //< track whther the transformation  is set up or not yet
 	rect_element clip_rect;//!< Points on clip path (used for plot window).

 	std::string  pl_window_clip_;//!< Points on clip path (used for plot window).

	void clip_path(const rect_element& rect, const std::string& id)	{ //! Rectangle outside which 'painting' is 'clipped' so doesn't show.
		 clip_rect = rect;
		 pl_window_clip_=id;
	}
	virtual void update_image() = 0;
	double x_size() {return Dx_;}
	double y_size() {return Dy_;}
	svchart& x_size(double d) { Dx_=d; return *this; }
	svchart& y_size(double d) { Dy_=d; return *this; }
	void write(std::ostream& os)  { //!< update plot and write to ostream.
		os<<	"<g  id=\"plot"<<int(X0_)<<int(Y0_)<<"\" transform=\"translate("<< X0_<<", "<< Y0_<<")\">\n";
		this->update_image();
		os << "<clipPath id=\"" << pl_window_clip_ << "\">";
		clip_rect.write(os);
		os  << "</clipPath>" << std::endl;

		 //! Write all visual group elements: pl_background, grids, axes ... title.
		 for(size_t i = 0; i < size(); ++i)
			(*this)[i].write(os);
		os<<	"</g>\n";
	}
 };

class svgraphic
{ /*! \class svgraphic
	\brief Class to output Scalable Vector Graph XML elements: point, path, line, circle, rect, polygon and text.
	\details  Class to add basic Scalable Vector Graph XML graph elements:
	point, path, line, circle, rect, polygon and text to SVG images,
	including metadata like author, copyright and license.
	Finally output the final image as SVG XML to a @c std::stream or file.
*/
 protected:
	unsigned int x_size_; //!< SVG image X-axis size (in SVG units (default pixels).
	unsigned int y_size_; //!< SVG image Y-axis size (in SVG units (default pixels).
	unsigned int n_x_;
	unsigned int n_y_;
	std::vector<svchart*> subplots_; //!< To hold all group elements of the svg document.

	//! Document metadata:
	std::string title_document_; //!< This is a SVG document title, not a plot title. //!< SVG document title (appears in the SVG file header as \verbatim <title> ... </title> \endverbatim).
	std::string image_desc_; //!< Information about the SVG image, for example, the program that created it. //!< SVG image description (appears in the SVG file header as \verbatim <desc> ... </desc> \endverbatim).
	std::string css_;  //!< Stylesheet filename.//!< Cascading Style Sheet.

	int coord_precision_; //!< Number of decimal digits precision for output of X and Y coordinates to SVG XML.


 public:
	svgraphic(int nI=1, int nJ=1, int dx=480, int dy=360) //! Define default constructor.
	:	x_size_(dx*nI), //!< X-axis of the whole SVG image (default SVG units, default pixels). Sync with svchart,
		y_size_(dy*nJ), //!< Y-axis of the whole SVG image (default SVG units, default pixels). Sync with svchart,
		n_x_(nI), n_y_(nJ),
		subplots_(nI*nJ,0),
		coord_precision_(3) //!< 3 decimal digits precision is enough for 1 in 1000 resolution: suits small image use. Higher precision (4, 5 or 6) will be needed for larger images, but increase the SVG XML file size, especially if there are very many data values.
	{  } // Default constructor.

	~svgraphic()
	{
		for(auto& plt: subplots_)  if(plt) { delete plt; plt=0; }
	}


	void x_size(unsigned int x)	{		x_size_ = x; 	} //! Set X-axis (horizontal) image size.

	void y_size(unsigned int y)	{ y_size_ = y; 	} //! Set Y-axis (vertical) image size.

	unsigned int x_size() 	{	return x_size_;	} //! \return  X-axis (horizontal width) SVG image size.

	unsigned int y_size() {	return y_size_; } //! \return  Y-axis (vertical height) SVG image size.

	std::pair<double, double> xy_sizes() { 		return std::pair<double, double>(x_size_,y_size_); 	}//! \return Both X and Y sizes (horizontal width and vertical height) of the SVG image.


	unsigned int nsubplots() { return static_cast<unsigned int>(subplots_.size());  } //! \return How many group elements groups have been added to the document.


	template<class subplotType>
	subplotType& subplot(int iX, int iY=0)
	{ //! \return How many group elements groups have been added to the document.
		int iPlt = iX+iY*n_x_;
		if(subplots_[iPlt]==0)
		{
			if(iPlt>=int(subplots_.size())) std::cout<<"Error: requesting too large plot index from svg "<<iPlt<<" actual size: "<<n_x_<<"*"<<n_y_<<std::endl;;
			subplots_[iPlt]=new subplotType(x_size_/n_x_-2, y_size_/n_y_-2, iX*x_size_/n_x_+1, iY*y_size_/n_y_+1);
		} //int Error_memory_leak risk;
		return *static_cast<subplotType*>(subplots_[iPlt]);
	}

	void coord_precision(int digits)
	{ //! \brief Set decimal digits to be output for X and Y coordinates.
		//! \details Default stream precision 6 decimal digits is probably excessive.\n
		//!  4.1 Basic data types, integer or float in decimal or scientific (using e format).
		//!  3 or 4 probably enough if image size is under 1000 x 1000.
		//!  This will reduce .svg file sizes significantly for curves represented with many data points.\n
		//!  Used in @c svgraphic.write below and so applies to all the entire @c svg document.
		coord_precision_ = digits;
	}

	int coord_precision() {		return coord_precision_; 	} //! \return  Decimal digits to be output for X and Y coordinates.

	void write(const std::string& filename)  {
		/*! \brief Write whole .svg 'file' contents to file.
		\details @c svgraphic.write() also has two flavors, a file and an ostream.
		The file version opens an ostream, and calls the stream version.
		The stream version first clears all unnecessary data from the graph,
		builds the document tree, and then calls the write function for the root
		document node, which calls all other nodes through the Visitor pattern.
		TODO provide a filtered-stream version that writes in zipped format type .svgz ?
		http://lists.w3.org/Archives/Public/www-svg/2005Dec/0308.html
		recommends MUST have  correct Content-Encoding headers.
		*/

		std::string file(filename); // Copy to avoid problems with const if need to append.
		if (file.find('.') == std::string::npos) { // No file type suffix, so provide the default .svg.
		  file.append(".svg");
		}
		std::ofstream f_out(file);
		if(f_out.fail())  throw std::runtime_error("Unable to open file " + file);

		write_begin(f_out);
		for(auto* splt:subplots_) if(splt) splt->write(f_out);
		f_out << "</svg>" << std::endl;
	}

	void write_begin(std::ostream& s_out)
	{ //! Write whole .svg 'file' contents to stream (perhaps a file).

		//! Output the DTD SVG 1.1 header into the svg g_element document.
		s_out << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>"<< std::endl;


		// write svg document, begin <svg tag.
		// <svg xml:space="preserve" width="5.5in" height=".5in">

		s_out << "<svg width=\"" << x_size_ << "\" height =\"" << y_size_
		  << "\" version=\"1.2\"\n" // http://www.w3.org/TR/SVG11/

		  // 1.2 draft specification at http://www.w3.org/TR/SVG12/
		  "xmlns:svg =\"http://www.w3.org/2000/svg\"\n"
		  "xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\"\n"
		  "xmlns:cc=\"http://web.resource.org/cc/\"\n"
		  "xmlns:dc=\"http://purl.org/dc/elements/1.1/\"\n"
		  "xmlns =\"http://www.w3.org/2000/svg\"\n"
		  // xml namespace containing svg shapes rect, circle...
		  // so can write rect or circle avoiding need for qualification svg:rect, svg:circle...
		  // This site isn't visited, but if missing Firefox, at least, will fail to render.
		  // Might also need xlink and ev,
		  // but Inkscape doesn't provide it, so we don't until required.
		  //   "xmlns:xlink=\"http://www.w3.org/1999/xlink\"\n"
		  // tells that elements and attributes which are prefixed by "xlink:"
		  // are a part of the xlink specification http://www.w3.org/1999/xlink.
		  // Need to use xlink:href to refer to xlink.
		  //  "xmlns:ev=\"http://www.w3.org/2001/xml-events\"\n"
		  << '>' << std::endl;

		s_out << "<!-- Written using svplot:  https://github.com/aliraeini/svplot -->\n" << std::endl;

		if (image_desc_ != "")
		{
		  s_out << "<desc>" << image_desc_ << "</desc>" << std::endl;
		}
		if (title_document_ != "")
		{
		   s_out << "<title>" << title_document_ << "</title>" << std::endl;
		}

		s_out.precision(coord_precision());
		//licence_.write(s_out); //. TODO

		// Output CSS (Cascading Style Sheet) - (not yet used or implemented).
		if (css_.size() != 0) // css != ""
		{
		  // [CDATA[ ... ]] enclosing the style information is a standard XML construct for hiding information  necessary since CSS style sheets can include characters, such as ">", which conflict with XML parsers.
		  s_out << "<defs><style type=\"text/css\"><![CDATA[" << css_ << "]]></style></defs>" << std::endl;
		  // CSS inline style can be declared within a style attribute in SVG
		  // by specifying a semicolon-separated list of property declarations,
		  // where each property declaration has the form "name: value".
		  // For example a style: style="fill:red; stroke:blue; stroke-width:3"
		  // class=
		  // Multiple class names must be separated by whitespace.
		  // Example: <defs><style type="text/css"><![CDATA[]]>
		  // .axes { fill:none;stroke:#333333;stroke-width:1.6 }
		  // .title{ font-size:20px;font-weight:bold;font-style:italic;fill:magenta }
		  // .legend_header{ font-size:16px;font-weight:bold;fill:darkblue;text-anchor:middle }
		  // .legend_item{ font-size:16px;font-weight:normal;fill:blue }
		  // .x_ax_value{ font-size:12px;font-weight:normal;fill:black }
		  //   </style></defs>
		}

	}

	void size(unsigned int x, unsigned int y)	{	x_size_ = x;	y_size_ = y;	/*! Set both X and Y image size (SVG units, default pixels).*/	}


	void description(const std::string d)	{	image_desc_ = d; 	/*! Write description to the SVG document (for header as \<desc\> ... \</desc\>).*/}



	//! \return  description of the SVG document (for header as \<desc\>).
	const std::string& description(){ 	return image_desc_;	}



	//!  Set document title for the SVG document (for header as \<title\> ... \</title\>).
	void document_title(const std::string title){  title_document_ = title;  }

	//!   \return document title for the SVG document (for header as \<title\>).
	const std::string document_title() {  return title_document_;  }


	//void load_stylesheet(const std::string& input)
	//{ // Load a stylesheet into string css from an input file.
	//  std::ifstream if_str(input.c_str());
	//  if(if_str.fail())  {    throw std::runtime_error("Error opening file " + input);  }
	//  if(!validate_stylesheet(if_str))  {   throw std::runtime_error("Error loading stylesheet!");  }
	//  if_str.clear();
	//  if_str.seekg(0);
	//  std::string tmp;
	//  css = "";
	//  while(std::getline(if_str, tmp))  {    css += tmp;  }
	//} // svg& load_stylesheet
}; // class svg





// svg_licence deactivated toue to low priority, need some re-arrangement to reactivate
#define BOOST_SVG_LICENSE_HPP

#ifndef BOOST_SVG_LICENSE_HPP
#define BOOST_SVG_LICENSE_HPP

class svg_licence{
	std::string holder_copyright_; //!< SVG info on holder of copyright (probably == author, but could be an institution).
	std::string date_copyright_; //!< SVG info on date of copyright claimed.
	std::string author_; //!< Author(s) name. (Probably == copyright holder).
	bool is_license_; //!< If true, to include Creative Commons license as metadata:
	std::string reproduction_; //!< License requirements for reproduction: "permits", "requires", or "prohibits".
	std::string attribution_; //!< License requirements for attribution: "permits", "requires", or "prohibits".
	std::string commercialuse_; //!< License requirements for commerical use: "permits", "requires", or "prohibits".
	std::string distribution_; //!< License requirements for distribution: "permits", "requires", or "prohibits".
	std::string derivative_works_; //!< License requirements for derivative: "permits", "requires", or "prohibits".


	svg_licence()
	:
		holder_copyright_(""),  //!< Name of copyright holder.
		date_copyright_(""), //!<  Date of copyright claim.
		author_(""), //!< Author of image (defaults to the copyright holder).
		is_license_(false), //!< If true, license text is written as comment in SVG XML. (Default is no license).
		reproduction_("permits"), //!< Default license permits reproduction.
		attribution_("requires"), //!< Default license requires attribution.
		commercialuse_("permits"), //<! Default license permits commerical use.
		distribution_("permits"), //!< Default license permits distribution.
		derivative_works_("permits"), //!< Default license permits derivative works.
	{}




	void license(
		const std::string reproduction = "permits",
		const std::string distribution = "permits",
		const std::string attribution = "requires",
		const std::string commercialuse = "permits",
		const std::string derivative = "permits")
	{  /*! Set several license requirements for the svg document.
		   If any are set, then a license is wanted, so @c svg#is_license is set @c true.
		   This can be changed using function @c license_on().
		 */
		reproduction_ = reproduction;
		distribution_ = distribution;
		attribution_ = attribution;
		commercialuse_ = commercialuse;
		derivative_works_ = derivative;
		is_license_ = true;  // Assume want this if set any of these requirements.
	} // void license


	void license_on(bool l)
	{ /*! Set (or not) license using all requirements (default permits).\n
		 Implicitly set by setting any license requirement using @c license function.
		 */
		is_license_ = l;
	}

	bool license_on()
	{ //! Return true if a license has been requested for @c svg header metatadata.
		return is_license_;
	}

	const std::string& reproduction()
	{ //! \return  License reproduction requirement.
		return reproduction_;
	}

	const std::string& distribution()
	{ //! \return  License distribution requirement.
		return distribution_;
	}
 //! \return  License attribution requirement.
	const std::string& attribution()	{		return attribution_;	}

 //! \return  License commercial use requirement.
 const std::string& commercialuse()	{		return commercialuse_;	}




		//!  Set author for the SVG document (default is \<copyright_holder\>).
	void author(const std::string a)		{ author_ = a;	}


		//!  \return  author of the SVG document (for header as \<author\>).
	const std::string& author()			{return author_;	}
		//!  Set document title for the SVG document (for header as  \<copyright_holder\>).
	void copyright_holder(const std::string copyright_holder)	{	holder_copyright_ = copyright_holder;	}

		//!  \return  document title for the SVG document (for header as  \<copyright_holder\> ).
	const std::string copyright_holder()	{	return holder_copyright_;	}

	//!  Set copyright date for the SVG document (for header as \<copyright_date\>).
	void copyright_date(const std::string copyright_date){ 		date_copyright_ = copyright_date;	}


	const std::string copyright_date(){		return date_copyright_;		/*!  \return copyright date for the SVG document (for header as \<copyright_date\>).*/}




	write_header(std::ostream& s_out)
	{

		if (author_ == "")
		{
		  author_ = holder_copyright_;
		}
		else
		{
		  if (holder_copyright_ == "")
		  {
		    holder_copyright_ = author_;
		  }
		  else
		  { // Copyright has been assigned to another, so list separately.
		    s_out << "<!-- " << author_ << " --> "<< std::endl;
		  }
		}
		if (holder_copyright_ != "")
		{ // Output copyright & date as both comment and meta data.
		  s_out << "<!-- SVG Plot Copyright " << holder_copyright_ << " " << date_copyright_ << " --> "<< std::endl;
		  s_out << "<meta name=\"copyright\" content=\"" << holder_copyright_ << "\" />" << std::endl;
		  //  Example:  \verbatim <meta name="copyright" content="Paul A. Bristow" /> \endverbatim
		  s_out << "<meta name=\"date\" content=\"" << date_copyright_ << "\" />" << std::endl;
		  // Example: \verbatim <meta name="Date" content="20071101"> \endverbatim
		}



		if (filename_ != "")
		{ // Example: <!-- File demo_1d_plot.svg -->
		  s_out << "<!-- File " << filename_ << " --> "<< std::endl;
		}


		if (is_license_ == true)
		{ // Add license information to the file.
		  // http://dublincore.org/documents/2000/07/16/usageguide/
		  // http://dublincore.org/documents/2000/07/16/usageguide/sectc.shtml#creator
		  s_out <<
		    "<metadata id = \"id0\">\n"
		      "<rdf:RDF>\n"
		         "<cc:Work rdf:about=\"" << filename_ << "\">\n" // Presumably .svg (or svgz?)
		           "<dc:format>image/svg+xml</dc:format>\n"
		           "<dc:type rdf:resource=\"http://purl.org/dc/dcmitype/StillImage\" />\n"
		           "<dc:title> " << (title_document_ != "" ? title_document_ : filename_) << "</dc:title>\n"
		           "<dc:creator> <cc:Agent> <dc:title>Boost.Plot</dc:title> </cc:Agent></dc:creator>\n"
		           "<dc:author><cc:Agent><dc:title>" << author_ << " </dc:title> </cc:Agent> </dc:author>\n"
		           "<dc:rights><cc:Agent><dc:title>" << holder_copyright_ << "</dc:title></cc:Agent></dc:rights>\n"
		           "<dc:date>" << date_copyright_ << "</dc:date>\n"
		           "<dc:identifier>" << filename_ << "</dc:identifier>\n" // URI for this svg document.
		           "<dc:source>" << "Boost.plot 0.5" << "</dc:source>\n"
		           "<dc:relation>" << "" << "</dc:relation>\n" // URI to a related document, perhaps user source program.
		           "<dc:publisher><cc:Agent><dc:title>" << holder_copyright_ << "</dc:title></cc:Agent></dc:publisher>\n"
		           "<dc:language>en_US</dc:language>\n" // Could be changed to suit, en-GB for example ;-)
		           "<dc:description>" << image_desc_ << "</dc:description>\n"
		           "<dc:contributor><cc:Agent><dc:title>" << author_ << "</dc:title></cc:Agent></dc:contributor>\n"
		           "<dc:subject><rdf:Bag><rdf:li>Boost svg plot keyword</rdf:li></rdf:Bag></dc:subject>\n"
		           // Could add keywords string here.
		           // License conditions URI: /by/ or /by_na/ ..
		           "<cc:license rdf:resource=\"http://creativecommons.org/licenses/\" />\n"
		           // Might instead select a specific license like http://creativecommons.org/licenses/by/3.0/  rather than a fully fexible combination as below.  Inkscape does this, for example.
		         "</cc:Work>\n"
		         "<cc:License rdf:about=\"http://creativecommons.org/licenses/\">\n"
		           "<cc:" << reproduction_ << " rdf:resource=\"http://web.resource.org/cc/Reproduction\"/>\n"
		           "<cc:" << distribution_ << " rdf:resource=\"http://web.resource.org/cc/Distribution\"/>\n"
		           "<cc:requires rdf:resource=\"http://web.resource.org/cc/Notice\"/>\n"
		           "<cc:" << attribution_ << " rdf:resource=\"http://web.resource.org/cc/Attribution\"/>\n"
		           "<cc:" << commercialuse_ << " rdf:resource=\"http://web.resource.org/cc/CommercialUse\"/>\n"
		           "<cc:" << derivative_works_ << " rdf:resource=\"http://web.resource.org/cc/DerivativeWorks\"/>\n"
		         "</cc:License>\n"
		      "</rdf:RDF>\n"
		     "</metadata>"
		   << std::endl;
		} // is_license
	}


}


#endif // BOOST_SVG_LICENSE_HPP





								_end_of_(namespace svg)
//} // namespace boost

#endif // BOOST_SVGRAPHIC_HPP
