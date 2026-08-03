/*-------------------------------------------------------------------------*\

This file is part of libvoxel, a C++ template library for handelling 3D images.

Developed by:
 - Ali Q Raeini (2010-2022)

\*-------------------------------------------------------------------------*/

#include <memory>
#include <sstream>


#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "voxelPng_stbi.h"

#include "globals.h"  // ensure...
#include "voxelImage.h"
#include "voxelImageI.h"
#include "voxelEndian.h"


#include "InputFile.h"
#include "voxelImgUtils.h"

#ifdef _POCKETPY
namespace VxlPy { void execScript(const std::string& script, VoxelImagesBase* imgPtr, const std::string& nam); }
#endif

using namespace VoxLib;
using namespace std; //cin cout endl string stringstream  istream istringstream regex*

// use `reset maxNzGlobal 100` to limit the number of layers processed during fine-tuning of image processing params
int maxNzGlobal = 1000000, _maxNz = 1000000|12;


string plotAll_normalAxis="xyz";
int    plotAll_colrGrey=15;

template<typename T>
void VoxelImageT<T>::readFromHeader(const string& hdrNam, int procesKeys, int maxNz)  {
  //! read image from file header, format detected based on image extension

  maxNz = maxNz>0 ? std::min(maxNzGlobal, maxNz) : maxNzGlobal;

  if (hdrNam.empty() || hdrNam=="NO_READ")  return;

  //! read image from file header, format detected based on image extension
  auto& vImg=*this;

  string fnam;
  int3 nnn(0,0,0);

  string BinaryData="XXX", flipSigByt="False";
  bool X0read=false, dxread=false, autoUnit=true; //auto unit only applies to .mhd format
  double unit_=1.;
  int nSkipBytes(0);
  if (hasExt(hdrNam,".raw.gz") || hasExt(hdrNam,".raw") || hasExt(hdrNam,".dat") || hasExt(hdrNam,".txt"))  { // detect size and voxel size from image name.
      // vxlProcess(data,vImg,hdrNam); // removed from Python build for simplicity
      procesKeys=0;
      fileNameToImgInfo(hdrNam, nnn, vImg.dx_);
      // TODO support units nm and mm ?
      fnam = hdrNam;
  }
  else if (hasExt(hdrNam,".mhd") || hasExt(hdrNam,".py")) {
    getMhdHeaderData(hdrNam, fnam, nnn, vImg.dx_, vImg.X0_, nSkipBytes,
                     BinaryData, flipSigByt, unit_, X0read, dxread, autoUnit,
                     maxNz, sizeof(T));
  }
  #ifdef TIFLIB
  else if (hasExt(hdrNam,".tif"))  {  readTif(vImg, hdrNam, maxNz);  return;  }
  #endif
  else if (hasExt(hdrNam,".am"))  {
    fnam=hdrNam;
    procesKeys=0;
  }

  if (nnn.z > 1.1 * maxNz) { nnn.z = maxNz; }
  if(nnn.z) vImg.reset(nnn);
  int readingImage=0;
  if( !fnam.empty() && fnam!="NO_READ" && procesKeys!=2)  {
    if (hasExt(fnam,".tif")) {
      dbl3 dx=vImg.dx_, X0=vImg.X0_;
      readingImage = vImg.readBin(fnam, 0, maxNz);
      if(X0read) vImg.X0_=X0;
      if(dxread) vImg.dx_=dx;
    }
    else if ((hasExt(fnam,".raw") && BinaryData!="False") || BinaryData=="True")   {
      readingImage = vImg.readBin(fnam, nSkipBytes, maxNz);
    }
    else if (hasExt(fnam,".am"))    {
      int RLECompressed;
      dbl3 dx = vImg.dx_, X0 = vImg.X0_;
      getAmiraHeaderSize(fnam, nnn, vImg.dx_, vImg.X0_, nSkipBytes, RLECompressed);
      readingImage = vImg.readBin(fnam, nSkipBytes, maxNz); // VoxelField::readBin has no dX nor X0
      if(X0read) vImg.X0_=X0;
      if(dxread) vImg.dx_=dx;
    }
    else if (hasExt(fnam,".raw.gz")) {
      readingImage = vImg.readBin(fnam, 0, maxNz);
    }
    else {
      ensure(hasExt(fnam, ".dat") || hasExt(fnam, ".txt"), "assuming image is in ascii format");
      readAscii(fnam, nSkipBytes);
    }
  }
  ensure(readingImage==0, "cannot read image "+fnam,-1);

  if (procesKeys) {
    #ifdef _POCKETPY
    InputFile fil(hdrNam);
    VxlPy::execScript(fil.kwrd("VxlPy"), &vImg, hdrNam);
    #else
      cout<<"\n  procesKeys is not added to cpython builds"<<endl;
    #endif
  }

  if(autoUnit  && vImg.dx_[0]>1e-2)  {
    cout<<"\n\n  WARNING dx="<<vImg.dx_[0]<<"(>1e-2 -> assuming unit is um),\n  please set Unit manually if needed ****\n\n";
    unit_ = 1e-6;
  }
  if(abs(unit_-1.)>epsT(float) && (!autoUnit || vImg.dx_[0]>1e-2)) {
    vImg.dx_*=unit_;
    vImg.X0_*=unit_;
    cout<<"\n  unit= "<<unit_<<" => dx= "<<vImg.dx_<<", X0= "<<vImg.X0_<<endl;
  }
  cout<<"."<<endl;
}

template<class T>
std::unique_ptr<VoxelImagesBase> readToUnique(const string& hdrNam, int procesKeys, int maxNz) {
  // NB! this is also called from VoxelImage read constructor, please avoid syclic dependency
  VoxelImageT<T> vImg;
  vImg.readFromHeader(hdrNam, procesKeys, maxNz);
  return make_unique<VoxelImageT<T>>(std::move(vImg));
}


std::unique_ptr<VoxelImagesBase> readImage(string hdrNam,  int procesKeys, int maxNz)  {
  //! read or create image
  using namespace std;
  (cout<<"VoxelImage \""<<hdrNam<<"\": ").flush();

  if (hasExt(hdrNam, ".png")) { // grey-scale atm
    VoxelImage VImage;
    sliceFromPng(VImage, "z", hdrNam, 0, 0,255);
    return make_unique<VoxelImageT<unsigned char>>(std::move(VImage));
  }

  if (hasExt(hdrNam,".am")) {
    string vtype = getAmiraDataType(hdrNam);
    cout<<"reading '"<<vtype<<"'s from .am file"<<endl;

    IF_MACRO_NOTVX8_ONLY( // if consteval may instantiate the expensive template(?)
    if (vtype=="int")       return readToUnique<int>(hdrNam, 0, maxNz);
    if (vtype=="ushort")    return readToUnique<unsigned short>(hdrNam, 0, maxNz);
    ) // #ifndef VX8_ONLY
    IF_MACRO__ExtraVxlTypes(
    if (vtype=="short")     return readToUnique<short>(hdrNam, 0, maxNz);
    ) // #ifdef _ExtraVxlTypes
    if (vtype=="byte" ||
        vtype=="uchar")
                            return readToUnique<unsigned char>(hdrNam, 0, maxNz);

    alert("data type "+vtype+" not supported, when reading "+hdrNam, -1);
  }

  #ifdef TIFLIB
  if (hasExt(hdrNam,".tif"))  return readTifAnyT(hdrNam,maxNz);
  #endif

  string typ;
  std::ifstream fil(hdrNam); // header file
  if (!fil) {
    ensure(hdrNam.size()<4 || hdrNam[hdrNam.size()-4]!='.', "can not open header file '"+hdrNam+"', pwd: "+getpwd(), -1);
    typ = hdrNam; hdrNam="NO_READ";
  }
  else
     typ = fileNameToElemType(hdrNam);

  fil.close();

  if (typ=="MET_UCHAR")        return readToUnique<unsigned char >(hdrNam, procesKeys, maxNz);
  #ifndef VX8_ONLY
  if (typ=="MET_USHORT")       return readToUnique<unsigned short>(hdrNam, procesKeys, maxNz);
  if (typ=="MET_INT")          return readToUnique<int>           (hdrNam, procesKeys, maxNz);
  if (typ=="MET_FLOAT")        return readToUnique<float>         (hdrNam, procesKeys, maxNz);
  #endif
#ifdef _ExtraVxlTypes
  if (typ=="MET_CHAR")         return readToUnique<char>          (hdrNam, procesKeys, maxNz);
  if (typ=="MET_SHORT")        return readToUnique<short>         (hdrNam, procesKeys, maxNz);
  if (typ=="MET_UINT")         return readToUnique<unsigned int>  (hdrNam, procesKeys, maxNz);
  if (typ=="MET_DOUBLE")       return readToUnique<double>        (hdrNam, procesKeys, maxNz);
  if (typ=="MET_FLOAT_ARRAY")  return readToUnique<float3>        (hdrNam, procesKeys, maxNz);
  if (typ=="MET_DOUBLE_ARRAY") return readToUnique<dbl3>          (hdrNam, procesKeys, maxNz);
#endif //_ExtraVxlTypes
  return                              readToUnique<unsigned char> (hdrNam, procesKeys, maxNz);

}
