/*-------------------------------------------------------------------------*\
For further information please contact Ali Q Raeini:    a.q.raeini@gmail.com
\*-------------------------------------------------------------------------*/


#include "typses.h"
#include "voxelImageI.h" // IWYU pragma: keep
#include "InputFile.h"


class poroRange : public int2  {
 public:
	poroRange(std::string nam, int lower,int upper)
	: var2<int>(lower,upper),name(nam), phi(1.,1.), Kq(1e32,1e32), Kef(1.,1.){};
	bool outside(int value){return this->a > value  ||  value > this->b; };
	std::string name;
	dbl2 phi, Kq, Kef;
};



void inline KozenyCarman(const InputFile& inp, const voxelImage& rokImg, std::vector<int>& segValues, std::vector<poroRange>& RTs)  {
	///. \param RTs: vValue-index pairs segValues
	using namespace std;
	using std::string;
	using std::endl;
	using std::vector;

	int nRTypes(0);
	istringstream ins;
																						cout<<" rock type  indices:"<<endl;

	RTs.push_back(poroRange("void", 0, 0));								cout<<"  "<<0<<": void voxels "<<endl;
	if(inp.giv("porousRocks", ins))  {
		ins >> nRTypes;    cout<<"nRTypes: "<<nRTypes<<"\n";
		for(int i = 1; i <=  nRTypes; ++i) {
			string nam;   ins>>nam;  											cout<<"  "<<i<<" name "<<nam<<endl;
			RTs.push_back(poroRange(nam,i,i));
		}
	}																					cout<< "  number of rock types: "<<int(RTs.size())<<endl;

	segValues.resize(256, RTs.size());

	for(int i=0; i<=nRTypes; ++i)  {
		auto& rt = RTs[i];
		if(inp.giv(rt.name+"_range", rt)) 									ensure(rt.a <= rt.b, "Bad keyword '"+rt.name+"_range',  first value is higher than second value");

		if(i==0)  rt.phi = {1., 1.};
		else if(inp.giv(rt.name+"_porosity", rt.phi)) 					ensure (rt.phi.a >= rt.phi.b, "Bad keyword '"+rt.name+"_porosity':  first value is lower than second value");
		else 																			alert("keyword \""+rt.name+"_porosity not found",-1);
	}


	vector<dbl2> SmRt;
	dbls volFracs(RTs.size(),0.);
	if(inp.giv("SwPc", ins))  {
		double sigma;   string fnam;   ins>>fnam>>sigma;;

		ifstream file(fnam);
		int nrows(0); file >> nrows;  SmRt.resize(nrows, dbl2(0,0));
		file>>SmRt; 								  cout<<SmRt<<endl;  	ensure(SmRt.size()>1,"couldn't read SwPc");
		for(dbl2& dd:SmRt) { dd.a = 1.-dd.a;  dd.b = 2.*sigma/dd.b;  }
																						ensure(SmRt[0].a >= SmRt.back().a, "Sw in SwPc should be monotonically decreasing ");

		for(int i=0; i<=nRTypes; ++i)  {
			const auto& rt = RTs[i];
			volFracs[i] = rokImg.volFraction(rt.a, rt.b)*0.5*(rt.phi.a+rt.phi.b);
		}
																						cout<<" volFracs: " << volFracs <<endl;
		volFracs/=sum(volFracs);
	}
	else 																				cout<<"\n\n  SwPc keyword not given " <<endl;

																						cout<<SmRt<<"\n\n\n"<<endl;
	double cdf=0.;
	for(int i=0; i<=nRTypes; ++i)  {
		auto& rt = RTs[i];
																						cout<<"\n\n*************************************  del: "<<volFracs[i]<<"+ cdf: "<<cdf<<endl;
		if(inp.giv(rt.name+"_permeability", rt.Kq))						ensure(rt.Kq.a >= rt.Kq.b,"Bad keyword '"+rt.name+"_permeability', a value is lower than b value");
		else if(i!=0) {
																						cout<<" using Carman-Kozeny equation for rock   "<<rt.name<<" K_q "<<endl;
																						ensure(SmRt.size(),"Error can not find "+rt.name+"_permeability  nor  SwPc  keywords",-1);
			int midT = (int(rt.a)+1+rt.b)/2;
			double f1= rokImg.volFraction(rt.a, midT)  *rt.phi.a;
			double f2= rokImg.volFraction(midT+1, rt.b)*rt.phi.b;   	cout<<"\t phi_i:      "<<f1<<",  "<<f2<<endl;
			f1 = volFracs[i]*(f1/(f1+f2));             					cout<<"\t phi_i/phi:  "<<f1<<"  "<<volFracs[i]-f1<<endl;
			double d1=2.*averageCDF(cdf, cdf+f1, SmRt);
			rt.Kq = {d1*d1*rt.phi.a/180., d1*d1*rt.phi.b/180.};
																						cout<<"\tRp  : "<<d1/2.<<" \t vfrac: "<<cdf<<" \t "<<cdf+f1<<"\t*******\n kq= 2Rp*2Rp*Phi/180 " <<endl;
		}

		//! formationFactor
		if     (inp.giv(rt.name+"_conductivity",    rt.Kef)) 			ensure(rt.Kef.a >= rt.Kef.b, "Bad keyword '"+rt.name+"_conductivity', first value is lower than second value");
		else if(inp.giv(rt.name+"_formationFactor", rt.Kef)) {
			rt.Kef = {1./(rt.Kef.a+1e-32), 1./(rt.Kef.b+1e-32)}; 		ensure(rt.Kef.a >= rt.Kef.b,"Bad keyword '"+rt.name+"_formationFactor', first value is lower than second value");
		}
		else if(i!=0)	{ 															ensure(SmRt.size(),"cannot find "+rt.name+"_formationFactor nor  SwPc  keywords",-1);
																						cout<<" using Carman-Kozeny equation for rock "<<rt.name<<", K_e: "<<inp.kwrd(rt.name+"_formationFactor")<<endl;
			rt.Kef = {sqr(rt.phi.a)/0.6, sqr(rt.phi.b)/0.6};
		}

		for(int j=rt.a; j<=rt.b; ++j)  segValues[j] = i;

		cdf+=volFracs[i];
																						cout<<"  "<< rt.name<<" voxels: ["<<rt<<"];   phi: ["<<rt.phi<<"];    Kef: ["<<rt.Kef<<"];   Kq: ["<<rt.Kq<<"];  cdf:"<<cdf<<endl;
	}
																						cout<<"\n Voxel value indices:";
																						for(size_t i = 0; i <=  7; ++i)  cout<<" "<<segValues[i];
																						cout<<" ... "<<endl;
}
