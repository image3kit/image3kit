#define _noAddFuncs

#define STB_IMAGE_WRITE_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
 #include "vxlPro1.cpp"

								namespace MCTProcessing _begins_
#ifndef _VoxBasic8
 template void addFuncs(VxlFuncs<unsigned short>&);
 template void addFuncs(VxlFuncs<int>&);
#endif
								_end_of_(namespace MCTProcessing)
