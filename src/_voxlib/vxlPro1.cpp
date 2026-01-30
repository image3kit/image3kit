/*-------------------------------------------------------------------------*\
You can redistribute this code and/or modify this code under the
terms of the GNU General Public License (GPL) as published by the
Free Software Foundation, either version 3 of the License, or (at
your option) any later version. see <http://www.gnu.org/licenses/>.

This file is part of voxelImage library, a C++ template library
developed by Ali Qaseminejad Raeini for handelling 3D raw images.

For further information please contact us by email:
Ali Q Raeini: a.q.raeini@gmail.com

\*-------------------------------------------------------------------------*/
#ifndef _VoxBasic
#include "voxelImageRegister.h" // order is important for typses.h
#endif //_VoxBasic

#include "voxelImage.h"
using namespace std; //cin cout endl string stringstream  istream istringstream regex*

#include "voxelImageProcess.h" // EXP:;

namespace MCTProcessing {
using namespace VoxLib;

//EXP_1:
template<typename T>
void addFuncs(VxlFuncs<T>& key_funs) requires(sizeof(T)>=3) {
		key_funs.insert({"end"   ,& end_here});
	}
template<typename T>
void addFuncs(VxlFuncs<T>& key_funs) requires(sizeof(T)<=2) {
	key_funs.insert({"averageWith_mBE"  , averageWith_mBE});
	key_funs.insert({"print_info"   , print_info});
	key_funs.insert({"mode26"       , mode26});
	key_funs.insert({"growingThreshold"  , growingThreshold});
	key_funs.insert({"labelImage"   , labelImage});
	key_funs.insert({"replaceOutSideValue"  , replaceOutSideValue});
	key_funs.insert({"smooth"       , smooth});
#ifndef _VoxBasic
	#ifndef SVG
	my bad, should have defined SVG C-macro
	#endif
	key_funs.insert({"svgHistogram" , svgHistogram});
	key_funs.insert({"svgZProfile"  , svgZProfile});
	key_funs.insert({"plotAll"  , plotAll});
	key_funs.insert({"flipEndian"   , flipEndian});
	key_funs.insert({"replaceRangeByImage"  , replaceRangeByImage});
	key_funs.insert({"replaceByImageRange"  , replaceByImageRange});
	key_funs.insert({"averageWith"  , averageWith});
	key_funs.insert({"readFromFloat", readFromFloat});
	key_funs.insert({"bilateralX"   , bilateralX});
	key_funs.insert({"bilateralGauss" , bilateralGauss});
	key_funs.insert({"meanWide"     , meanWide});
	key_funs.insert({"print_otsu"   , print_otsu});
	key_funs.insert({"segment"	    , segment});
	key_funs.insert({"segment2"     , segment2});
	key_funs.insert({"registerToImage" , registerToImage});
	key_funs.insert({"readTransFrom" , readTransFrom});
#endif // !_VoxBasic
	key_funs.insert({"dering"       , dering});
	key_funs.insert({"adjustBrightnessWith"  , adjustBrightnessWith});
	key_funs.insert({"adjustSliceBrightness" , adjustSliceBrightness});
	key_funs.insert({"cutOutside"   , vxlProCutOutside});
	key_funs.insert({"variance" , variance});
	key_funs.insert({"end"   , end_here});
}
//EXP_1;



#ifndef _noAddFuncs
 template void addFuncs(VxlFuncs<float>&);
 template void addFuncs(VxlFuncs<unsigned char>&);
#endif

} // namespace MCTProcessing
