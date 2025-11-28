#include <filesystem>
namespace fs = std::filesystem;


template<typename T>  bool voxelTimeMedian( std::stringstream & ins, voxelImageT<T>& vImage)
{
  string inDir;     ins>>inDir;
  string sfx;       ins>>sfx;
  string ouDir;     ins>>ouDir;
  int   kernelN;    ins>>kernelN;

  set<string> imgNams
  auto lx=sfx.size();
  cout<<"files:"<<endl;
  for(auto& p: fs::directory_iterator(inDir))
  {  auto lp=p.path().size();
    if(lp>lx&& p.path().compare(lp-lx,lx,sfx) == 0)
    {
      imgNams.insert(p.path());
      cout<<p.path()<<endl;
    }
    else cout<<p.path()<<", skipped"<<endl;
  }
  cout<<"running voxelTimeMedian on "<<imgNams.size()<<" images, kernelN: "<<kernelN<<endl;
  outputFolder = fs::create_directory(ouDir);
  ensure(outputFolder, "can not create new directory "+ouDir);

  ensure(kernelN==3);
  size_t mid = kernelN/2;
  vector<voxelImageT<T> > imgs(kernelN);
  imgs[0]=vImage;
  imgs[0].readBin(*imgNams.begin());
  for(int i=1; i<imgs.size(); ++i) imgs[i]=imgs[0];
  for(auto imgnam=imgNams.begin()+1; imgnam!=imgNams.end()-1; ++imgnam)
  {  int i=0
    for(; i<imgs.size()-1; ++i) imgs[i]=std::move(imgs[i+1]);
    imgs[i+1].readBin(inDir+"/"+imgnam);
    forAlliii_(vImage)
    {
      std::vector<T> vvs(imgs.size(),0);
      forAlli(imgs)
      {  vvs[i]=imgs[i](iii);
        //std::array<T,7> vvs={{
            //vImage.v_i(-nW,vp), vImage.v_i( nW,vp),
            //vImage.v_j(-nW,vp), vImage.v_j( nW,vp),
            //vImage.v_k(-nW,vp), vImage.v_k( nW,vp)
             //}};
        //voxls(i,j,k) = std::nth_element(vvs.begin(),vvs.begin()+3,vvs.end())[3];
      }
      vImage(iii) = std::nth_element(vvs.begin(),vvs.begin()+mid,vvs.end())[mid];
    }
    vImage.writeBin(ouDir+"/"+imgnam);
  }

  return false;
}
