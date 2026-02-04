#pragma once
#ifdef LPNG

#include "typses.h"
#include "voxelImage.h"
#include "IOUtils.h"

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <png.h>

#include <svg_color.hpp>

template <typename T>
int slice2RGBAPng(const voxelImageT<T>& imgRGB, char axis, std::string fnam, int iSlice, T bgnv, T endv, const voxelImageT<T>& imgAlf, T bgna, T enda)  {
	cout<<"Writing "<<fnam<<" \n";  // #WViz

	int width, height;
	auto nnn=imgRGB.size3();
	switch (axis) {
		case 'x':  width=nnn[2]; height=nnn[1]; ensure(iSlice<nnn[0]," i_slice ?<"+_s(nnn[0]),2); break;
		case 'y':  width=nnn[0]; height=nnn[2]; ensure(iSlice<nnn[1]," i_slice ?<"+_s(nnn[1]),2);  break;
		case 'z':  width=nnn[0]; height=nnn[1]; ensure(iSlice<nnn[2]," i_slice ?<"+_s(nnn[2]),2);  break;
	}

	FILE *fp = NULL;
	png_structp png_ptr = NULL;
	png_infop info_ptr = NULL;
	png_bytep row = NULL;

	// Open file for writing (binary mode)
	fp = fopen(fnam.c_str(), "wb");  ensure(fp, "Could not open file "+fnam+" for writing", -1);

	// Initialize write structure
	png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL); 	ensure (png_ptr,  "Could not allocate write struct",-1);

	// Initialize info structure
	info_ptr = png_create_info_struct(png_ptr); 	ensure (info_ptr,  "Could not allocate info struct\n",-1);

	// Setup Exception handling
	ensure(setjmp(png_jmpbuf(png_ptr))==0 , "Error during png creation\n",-1);

	png_init_io(png_ptr, fp);

	// Write header (8 bit colour depth)
	png_set_IHDR(png_ptr, info_ptr, width, height, 8, PNG_COLOR_TYPE_RGBA, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);


	png_write_info(png_ptr, info_ptr);

	// Allocate memory for one row (3 bytes per pixel - RGB)
	row = (png_bytep) malloc(4 * width * sizeof(png_byte));

	// Write image data
	T minvv=bgnv, maxvv=endv; 	if(minvv>maxvv) std::swap(minvv,maxvv);
	float ZdelV=1./(double(endv)-bgnv);
	T minaa=bgna, maxaa=enda; 	if(minaa>maxaa) std::swap(minaa,maxaa);
	float ZdelA=255./(double(enda)-bgna);
	int x, y;
	for (y=0 ; y<height ; y++) {
		for (x=0 ; x<width ; x++)
		{
			size_t ind=0;
			switch (axis) {
				case 'x': ind=imgRGB.index(iSlice,y,x); break;
				case 'y': ind=imgRGB.index(x,iSlice,y);  break;
				case 'z': ind=imgRGB.index(x,y,iSlice);  break;
			}
			float v = (float(min(max(minvv,imgRGB(ind)),maxvv))-bgnv)*ZdelV;
			//auto clr = svg::rgbtween(v, svg::svg_color(20,30,0), svg::svg_color(250,220,221)); // good for rock, goes through pink
			//auto clr = svg::rgbtween(v, svg::svg_color(0,0,100), svg::svg_color(255,110,100));
			auto clr = svg::rgbtween(v, svg::svg_color(10,11,150), svg::svg_color(210,190,170));
			row[x*4]=clr.r_; row[x*4+1]=clr.g_; row[x*4+2]=clr.b_; row[x*4+3]=(float(min(max(minaa,imgAlf(ind)),maxaa))-bgna)*ZdelA;
		}
		png_write_row(png_ptr, row);
	}

	// End write
	png_write_end(png_ptr, NULL);

	//finalise:
	if (fp)       fclose(fp);
	if (info_ptr) png_free_data(png_ptr, info_ptr, PNG_FREE_ALL, -1);
	if (png_ptr)  png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
	if (row)      free(row);

	return 0;
}


template <typename T>
int slice2RGBXPng(const voxelImageT<T>& imgRGB, char axis, std::string fnam, int iSlice, T bgnv, T endv, const voxelImageT<T>& imgAlf, T bgna, T enda, float wt)  {
	cout<<"Writing "<<fnam<<" \n";  // #WViz

	int width, height;
	auto nnn=imgRGB.size3();
	switch (axis) {
		case 'x':  width=nnn[2]; height=nnn[1]; ensure(iSlice<nnn[0]," i_slice ?<"+_s(nnn[0]),2); break;
		case 'y':  width=nnn[0]; height=nnn[2]; ensure(iSlice<nnn[1]," i_slice ?<"+_s(nnn[1]),2);  break;
		case 'z':  width=nnn[0]; height=nnn[1]; ensure(iSlice<nnn[2]," i_slice ?<"+_s(nnn[2]),2);  break;
	}

	FILE *fp = NULL;
	png_structp png_ptr = NULL;
	png_infop info_ptr = NULL;
	png_bytep row = NULL;

	// Open file for writing (binary mode)
	fp = fopen(fnam.c_str(), "wb");  ensure(fp, "Could not open file "+fnam+" for writing", -1);

	// Initialize write structure
	png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL); 	ensure (png_ptr,  "Could not allocate write struct",-1);

	// Initialize info structure
	info_ptr = png_create_info_struct(png_ptr); 	ensure (info_ptr,  "Could not allocate info struct\n",-1);

	// Setup Exception handling
	ensure(setjmp(png_jmpbuf(png_ptr))==0 , "Error during png creation\n",-1);

	png_init_io(png_ptr, fp);

	// Write header (8 bit colour depth)
	png_set_IHDR(png_ptr, info_ptr, width, height, 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);


	png_write_info(png_ptr, info_ptr);

	// Allocate memory for one row (3 bytes per pixel - RGB)
	row = (png_bytep) malloc(3 * width * sizeof(png_byte));
	float wv=1.0f-wt;
	// Write image data
	T minvv=bgnv, maxvv=endv; 	if(minvv>maxvv) std::swap(minvv,maxvv);
	float ZdelV=1./(double(endv)-bgnv);
	T minaa=bgna, maxaa=enda; 	if(minaa>maxaa) std::swap(minaa,maxaa);
	float ZdelA=127./(double(enda)-bgna);
	int x, y;
	for (y=0 ; y<height ; y++) {
		for (x=0 ; x<width ; x++)
		{
			size_t ind=0;
			switch (axis) {
				case 'x': ind=imgRGB.index(iSlice,y,x); break;
				case 'y': ind=imgRGB.index(x,iSlice,y);  break;
				case 'z': ind=imgRGB.index(x,y,iSlice);  break;
			}
			float vv = (float(min(max(minvv,imgRGB(ind)),maxvv))-bgnv)*ZdelV;
			auto clr = svg::rgbtween(vv, svg::svg_color(10,11,150), svg::svg_color(255,51,50));
			auto c2r = (float(min(max(minaa,imgAlf(ind)),maxaa))-bgna)*ZdelA;
			row[x*3]=clr.r_*wv+wt*c2r; row[x*3+1]=clr.g_*wv+wt*c2r; row[x*3+2]=clr.b_*wv+wt*c2r;
		}
		png_write_row(png_ptr, row);
	}

	// End write
	png_write_end(png_ptr, NULL);

	//finalise:
	if (fp)       fclose(fp);
	if (info_ptr) png_free_data(png_ptr, info_ptr, PNG_FREE_ALL, -1);
	if (png_ptr)  png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
	if (row)      free(row);

	return 0;
}


template <typename T>
int slice2RGBPng(const voxelImageT<T>& vImage, char axis, std::string fnam, int iSlice, T bgnv, T endv)  {

	cout<<"Writing "<<fnam<<" \n";  // #WViz

	int width, height;
	auto nnn=vImage.size3();
	switch (axis) {
		case 'x':  width=nnn[2]; height=nnn[1]; ensure(iSlice<nnn[0]," i_slice ?<"+_s(nnn[0]),2); break;
		case 'y':  width=nnn[0]; height=nnn[2]; ensure(iSlice<nnn[1]," i_slice ?<"+_s(nnn[1]),2);  break;
		case 'z':  width=nnn[0]; height=nnn[1]; ensure(iSlice<nnn[2]," i_slice ?<"+_s(nnn[2]),2);  break;
	}

	FILE *fp = NULL;
	png_structp png_ptr = NULL;
	png_infop info_ptr = NULL;
	png_bytep row = NULL;

	// Open file for writing (binary mode)
	fp = fopen(fnam.c_str(), "wb");  ensure(fp, "Could not open file "+fnam+" for writing", -1);

	// Initialize write structure
	png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL); 	ensure (png_ptr,  "Could not allocate write struct",-1);

	// Initialize info structure
	info_ptr = png_create_info_struct(png_ptr); 	ensure (info_ptr,  "Could not allocate info struct\n",-1);

	// Setup Exception handling
	ensure(setjmp(png_jmpbuf(png_ptr))==0 , "Error during png creation\n",-1);

	png_init_io(png_ptr, fp);

	// Write header (8 bit colour depth)
	png_set_IHDR(png_ptr, info_ptr, width, height, 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);


	png_write_info(png_ptr, info_ptr);

	// Allocate memory for one row (3 bytes per pixel - RGB)
	row = (png_bytep) malloc(3 * width * sizeof(png_byte));

	// Write image data
	T minvv=bgnv, maxvv=endv;  if(minvv>maxvv) std::swap(minvv,maxvv);
	float ZdelV=1./(double(endv)-bgnv);
	int x, y;
	for (y=0 ; y<height ; y++) {
		for (x=0 ; x<width ; x++)
		{
			size_t ind=0;
			switch (axis) {
				case 'x': ind=vImage.index(iSlice,y,x); break;
				case 'y': ind=vImage.index(x,iSlice,y);  break;
				case 'z': ind=vImage.index(x,y,iSlice);  break;
			}
			float v = float(min(max(minvv,vImage(ind)),maxvv)-bgnv)*ZdelV;
			//auto clr = svg::rgbtween(v, svg::svg_color(20,30,0), svg::svg_color(250,220,221)); // good for rock, goes through pink
			//auto clr = svg::rgbtween(v, svg::svg_color(0,0,100), svg::svg_color(255,110,100));
			auto clr = svg::rgbtween(v, svg::svg_color(140,40,150), svg::svg_color(250,220,210));
			row[x*3]=clr.r_; row[x*3+1]=clr.g_; row[x*3+2]=clr.b_;
		}
		png_write_row(png_ptr, row);
	}

	// End write
	png_write_end(png_ptr, NULL);

	//finalise:
	if (fp)       fclose(fp);
	if (info_ptr) png_free_data(png_ptr, info_ptr, PNG_FREE_ALL, -1);
	if (png_ptr)  png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
	if (row)      free(row);

	return 0;
}


template <typename T>
int slice2GrayPng(const voxelImageT<T>& vImage, char axis, std::string fnam, int iSlice, T bgnv, T endv){

	cout<<"Writing "<<fnam<<" \n";  // #WViz

	int width, height;
	auto nnn=vImage.size3();
	switch (axis) {
		case 'x':  width=nnn[2]; height=nnn[1]; ensure(iSlice<nnn[0]," i_slice ?<"+_s(nnn[0]),2); break;
		case 'y':  width=nnn[0]; height=nnn[2]; ensure(iSlice<nnn[1]," i_slice ?<"+_s(nnn[1]),2);  break;
		case 'z':  width=nnn[0]; height=nnn[1]; ensure(iSlice<nnn[2]," i_slice ?<"+_s(nnn[2]),2);  break;
	}

	FILE *fp = NULL;
	png_structp png_ptr = NULL;
	png_infop info_ptr = NULL;
	png_bytep row = NULL;

	// Open file for writing (binary mode)

	fp = fopen(fnam.c_str(), "wb");  ensure(fp, "Could not open file "+fnam+" for writing",-1);

	// Initialize write structure
	png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL); 	ensure (png_ptr, "Could not allocate write struct",-1);

	// Initialize info structure
	info_ptr = png_create_info_struct(png_ptr); 	ensure (info_ptr,  "Could not allocate info struct\n",-1);

	// Setup Exception handling
	ensure(setjmp(png_jmpbuf(png_ptr))==0 , "during png creation\n",-1);

	png_init_io(png_ptr, fp);

	// Write header (8 bit colour depth)
	png_set_IHDR(png_ptr, info_ptr, width, height, 8, PNG_COLOR_TYPE_GRAY, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);


	png_write_info(png_ptr, info_ptr);

	// Allocate memory for one row (3 bytes per pixel - RGB)
	row = (png_bytep) malloc(1 * width * sizeof(png_byte));

	// Write image data
	T minvv=bgnv, maxvv=endv; 	if(minvv>maxvv) std::swap(minvv,maxvv);
	float ZdelV=255./(double(endv)-bgnv);
	int x, y;
	for (y=0 ; y<height ; y++) {
		for (x=0 ; x<width ; x++)
		{
			size_t ind=0;
			switch (axis) { // move out of the loop
				case 'x': ind=vImage.index(iSlice,y,x); break;
				case 'y': ind=vImage.index(x,iSlice,y);  break;
				case 'z': ind=vImage.index(x,y,iSlice);  break;
			}
			row[x] = (float(min(max(minvv,vImage(ind)),maxvv))-bgnv)*ZdelV;
		}
		png_write_row(png_ptr, row);
	}

	// End write
	png_write_end(png_ptr, NULL);

	//finalise:
	if (fp)       fclose(fp);
	if (info_ptr) png_free_data(png_ptr, info_ptr, PNG_FREE_ALL, -1);
	if (png_ptr)  png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
	if (row)      free(row);

	return 0;
}


template <typename T>
static void fixData(const voxelImageT<T> &vImg, string& normalAxis,
                    std::string& fnam, int& iSlice, int& minv, int& maxv) {
	char axs = tolower(normalAxis[0]);
	normalAxis[0] = axs;

	int iDir=axs-'x'; iDir%=3;
	if (iSlice<-vImg.size3()[iDir]) iSlice=vImg.size3()[iDir]/2;

	if (maxv==-1000001)  { minv = accumulate(vImg,(std::min<T>));   maxv = accumulate(vImg,(std::max<T>));  }
	(cout<<" minmaxv: "<<minv<<" "<<maxv<<" ").flush();

	if (hasExt(fnam,".png"))  fnam = fnam.substr(0,fnam.size()-4);
	if (fnam[0]=='_') fnam=_s(axs)+_s(iSlice)+fnam+".png";
	else              fnam=fnam+"_"+_s(axs)+_s(iSlice)+".png";
	if (fnam.find('/')==std::string::npos) { mkdirs("fig"); fnam="fig/"+fnam; }
}


template<typename T>  bool sliceToPng(const voxelImageT<T> &vImg,
	string normalAxis, std::string fnam, int iSlice, int bgnv, int endv, std::string color="") {

	fixData(vImg, normalAxis, fnam, iSlice, bgnv, endv);

	if (color=="all") color="";
	if (color!="" && color!="RGB" && color!="gray") { cout<<"wrong color scheme, "<<color<<", expected 'gray' or 'RGB' or 'all'"<<endl;   color=""; }
	if (color=="" || color=="gray")
		slice2GrayPng(vImg, normalAxis[0], prepend("grey_",fnam), iSlice, T(bgnv), T(endv));  // grey_ is at begin for image slide-show order
	if (color=="" || color=="RGB")
		slice2RGBPng(vImg, normalAxis[0],                  fnam , iSlice, T(bgnv), T(endv));
	return 0;
}

template<typename T>  bool sliceMergeToPng(const voxelImageT<T> &vImg,
		string normalAxis, std::string fnam, int iSlice, int bgnv, int endv,
		const voxelImageT<T> &aImg, int mina, int maxa, float weight, std::string color="")
 {
	fixData(vImg, normalAxis, fnam, iSlice, bgnv, endv);
	if(maxa==-1000001)  { mina = accumulate(vImg,(std::min<T>));   maxa = accumulate(vImg,(std::max<T>));  }
	(cout<<"  minmaxa: "<<mina<<" "<<maxa<<"  w: "<<weight<<" ").flush();

	if (color=="all") color="";
	if (color!="" && color!="RGBA" && color!="RGBX") { cout<<"wrong color scheme, "<<color<<", expected 'RGBX' or 'RGBA' or 'all'"<<endl;   color=""; }
	if (color=="" || color=="RGBX")
		slice2RGBXPng(vImg, normalAxis[0], prepend("rgbx_",fnam), iSlice, T(bgnv), T(endv), aImg, T(bgnv), T(endv), weight);
	if (color=="" || color=="RGBA")
		slice2RGBAPng(vImg, normalAxis[0], prepend("rgba_",fnam), iSlice, T(bgnv), T(endv), aImg, T(bgnv), T(endv));
	return 0;
}


								namespace MCTProcessing _begins_


template<typename T>  bool sliceToPng( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("normalAxis(x/y/z) fileName iSlice  minv maxv  [al b e w]\n// if fileName starts with _, extra info are added as prefix otherwise as suffix to the filename");
	int minv(0); int maxv(-1000001); int iSlice(-1000000);
	string normalAxis("x");
	string fnam("aa.png"); string alphaImg;
	ins >> normalAxis;
	readOutFileNam(ins, fnam);
	ins >> iSlice >> minv >> maxv >> alphaImg;

	#ifdef _STOR_PUB
	if(len(alphaImg)) {
		int mina(0); int maxa(-1000001); float weight(0.9);
		ins  >> mina >> maxa >> weight;
		_dbgetOrReadImg(T,img2, alphaImg);
		sliceMergeToPng(vImg, normalAxis, fnam, iSlice, minv, maxv, img2, mina, maxa, weight, "");
	} //else
	#endif
	sliceToPng(vImg, normalAxis, fnam,  iSlice, minv, maxv, "");

	(cout<<".").flush();

	//slice2GrayPng(vImg,normalAxis[0],"grey_"+fnam+"_"+normalAxis+_s(iSlice)+".png", iSlice,T(minv),T(maxv));
	return true;
}


template<typename T>  bool sliceToPngBW( stringstream& ins, voxelImageT<T>& vImg) {
	if(ins.peek()=='?') { ins.str("normalAxis fileName iSlice  minv maxv"); return true;}
	int minv(0); int maxv(-1000001); int iSlice(-1000000);
	string normalAxis("x");
	string fnam("Img");
	ins >> normalAxis;
	readOutFileNam(ins, fnam);
	ins >> iSlice >> minv >> maxv;

	fixData(vImg, normalAxis, fnam, iSlice, minv, maxv);

	slice2GrayPng(vImg, normalAxis[0], prepend("grey_",fnam), iSlice,T(minv),T(maxv)); // grey_ is at begin for image slide-show order
	(cout<<".").flush();

	return true;
}


								_end_of_(namespace MCTProcessing)

#endif
