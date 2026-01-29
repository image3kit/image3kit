/*-------------------------------------------------------------------------*\
You can redistribute this code and/or modify this code under the
terms of the GNU General Public License (GPL) as published by the
Free Software Foundation, either version 3 of the License, or (at
your option) any later version. see <http://www.gnu.org/licenses/>.

For further information please contact us by email:
Ali Q Raeini:    a.q.raeini@gmail.com
\*-------------------------------------------------------------------------*/

#include <fstream>
#include <iostream>
#include <vector>

#include <assert.h>


#include "voxelImage.h"
#include "KozenyCarman.h"

using namespace std;

int usage()  {
	std::cout<<"Ufraw2Uc  "<<std::endl;
		std::cout
		<<"Compute the gradient of psi field\n"
		<<"Usage examples, type:\n"
		<<" #cd PATH/TO/Ufx.*etc first"<< std::endl
		<<" voxelGrad  psi.raw  vxlImage.mhd z am  rocklbls.mhd  "<< std::endl;
	return 1;
}


int main(int argc, char** argv)  {

	if(argc<=3)		return usage();
	char dir( argc>3 ? argv[3][0] : 'x');
	std::string outFormat( argc>4 ? argv[4] : ".raw.gz");
	imgExt(outFormat);
	std::string headerName(argv[2]);
	if(headerName!="vxlImage.mhd") return usage();

	string psimg;
	if(argc>1) {
		psimg=string(argv[1]);
	}else psimg="psi.raw.gz";
	cout<<"//! psimg: "<<psimg<<", headerName: "<<headerName<<",  dir: "<<dir<<", outFormat: "<<outFormat <<"\n"<<endl;

	string RqfileName = (argc>5) ? string(argv[5]): "";


	voxelImage vimage("vxlImage.mhd", readOpt::procOnly);




	if(!vimage.nz()) {std::cout<<"Error: vxlImage.mhd not read"<<std::endl; return 1;}
	vimage.growBox(1);
	int3 n=vimage.size3();


	voxelImageT<float> fField(n[0]-2,n[1]-2,n[2]-2,0.);
	const float p2dx=0.5/vimage.dx()[0];
	fField.readBin(psimg);
	fField.growBox(1);
	if(dir=='x')  {
	 for (int k=0; k<n[2]; ++k)
	  for (int j=0; j<n[1]; ++j)  {
		if(vimage(2     ,j,k)==0) fField(0     ,j,k)=2*fField(1     ,j,k)-fField(2     ,j,k);
		if(vimage(n[0]-2,j,k)==0) fField(n[0]-1,j,k)=2*fField(n[0]-1,j,k)-fField(n[0]-2,j,k);
	  }
	}
	else if(dir=='y')  {
	 for (int k=0; k<n[2]; ++k)
	  for (int i=0; i<n[0]; ++i)  {
		if(vimage(i,2     ,k)==0) fField(i,0     ,k)=2*fField(i,1     ,k)-fField(i,2     ,k);
		if(vimage(i,n[1]-2,k)==0) fField(i,n[1]-1,k)=2*fField(i,n[1]-1,k)-fField(i,n[1]-2,k);
	  }
	}
	else if(dir=='z')  {
	 for (int i=0; i<n[0]; ++i)
	  for (int j=0; j<n[1]; ++j)  {
		if(vimage(i,j,2     )==0) fField(i,j,0     )=2*fField(i,j,1     )-fField(i,j,2     );
		if(vimage(i,j,n[2]-2)==0) fField(i,j,n[2]-1)=2*fField(i,j,n[2]-1)-fField(i,j,n[2]-2);
	  }
	}
	else cout<<"Error wrong flow direction, "<<dir<<endl;


	voxelField<float3> pGrad(n[0],n[1],n[2],{0.0f,0.0f,0.0f});
	{
	 float maxMagGrad=0.;
	 forAllkji_1_(vimage)  {
		size_t iii=vimage.index(i,j,k);
		if(vimage(iii)==0)  {
			size_t iip=vimage.I_i(1,iii), iim=vimage.I_i(-1,iii);
			if(vimage(iip)==0 && vimage(iim)==0)
			  pGrad(iii)[0]=(fField(iip)-fField(iim))*p2dx;
			else if(vimage(iip)==0)
			  pGrad(iii)[0]=(fField(iip)-fField(iii))*p2dx;
			else if(vimage(iim)==0)
			  pGrad(iii)[0]=(fField(iii)-fField(iim))*p2dx;

			iip=vimage.I_j(1,iii), iim=vimage.I_j(-1,iii);
			if(vimage(iip)==0 && vimage(iim)==0)
			  pGrad(iii)[1]=(fField(iip)-fField(iim))*p2dx;
			else if(vimage(iip)==0)
			  pGrad(iii)[1]=(fField(iip)-fField(iii))*p2dx;
			else if(vimage(iim)==0)
			  pGrad(iii)[1]=(fField(iii)-fField(iim))*p2dx;

			iip=vimage.I_k(1,iii), iim=vimage.I_k(-1,iii);
			if(vimage(iip)==0 && vimage(iim)==0)
			  pGrad(iii)[2]=(fField(iip)-fField(iim))*p2dx;
			else if(vimage(iip)==0)
			  pGrad(iii)[2]=(fField(iip)-fField(iii))*p2dx;
			else if(vimage(iim)==0)
			  pGrad(iii)[2]=(fField(iii)-fField(iim))*p2dx;

			maxMagGrad=max(maxMagGrad, mag(pGrad(iii)));
		}
	 }
	 cout<< "\nmaxMagGrad:  "<<maxMagGrad << endl;
	}


	if(RqfileName.size())  {
		voxelImage rockLabels(RqfileName, readOpt::procOnly);
		rockLabels.growBox(1);

		cout<< "Reading " << RqfileName <<":"<<std::endl;
		InputFile input(RqfileName);

		std::vector<int> segValues;
		std::vector<poroRange> RTs; ///. vValue-index pairs segValues

		KozenyCarman(input,rockLabels,segValues,RTs);

		forAlliii_(rockLabels)  {
			if(segValues[rockLabels(iii)]<len(RTs))  {
				int ri=segValues[rockLabels(iii)];
				double Kef=RTs[ri].Kef.a
			          + ( ri-RTs[ri].a+1e-18) /(RTs[ri].b-RTs[ri].a+2e-18)
						 * (RTs[ri].Kef.b-RTs[ri].Kef.a);
				pGrad(iii)*=Kef;
			}
		 }
	}


	imgExt(outFormat);

	pGrad.writeHeader("grad_psi"+imgExt(),{1,1,1},n-int3{1,1,1}, vimage.dx(), vimage.dx()+vimage.X0());
	pGrad.writeBin("grad_psi"+imgExt(),1,n[0]-1, 1,n[1]-1, 1,n[2]-1);

	std::cout<< "end" << std::endl;


	return 0;
}
