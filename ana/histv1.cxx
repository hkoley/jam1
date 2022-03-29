#include <iostream>
#include <cstdlib>
#include <sstream>
#include <cstring>
#include <string>
#include <iomanip>
#include <fstream>
#include <cmath>

#include "Book1.h"


using namespace std;

int optHist=1;
int ntot=10000;
int iev, nv,nbary,nmeson;
double impactParam;
int nEvent=0;

double *mass,*px,*py,*pz,*pe,*vx,*vy,*vz,*ve,*vm;
int *kf;

int nbook=8;
Book1 *dndy;
Book1 *v1, *v2;

int ny=50;
double ymin=-3.0, ymax=3.0;
double wy = (ymax-ymin)/ny;

  string title1[] = {"dndy p","dndy anti-p","dndy pi-","dndy pi0","dndy pip",
                     "dndy nucleon", " dndy pion","dndy anti-N"};
  string title2[] = {"v1 p","v1 anti-p","v1 pi-","v1 pi0","v1 pip",
                       "v1 N","v1 pion","v1 anti-N"};
  string title3[] = {"v2 p","v2 anti-p","v2 pi-","v2 pi0","v2 pip",
                      "v2 N","v2 pion","v2 anti-N"};
void book1()
{

  dndy = new Book1[nbook];
  v1   = new Book1[nbook];
  v2   = new Book1[nbook];

  for (int i=0;i<nbook;i++) {
    dndy[i].setPar(title1[i],ny,ymin,ymax);
    v1[i].setPar(title2[i],ny,ymin,ymax);
    v2[i].setPar(title3[i],ny,ymin,ymax);
  }
}

void book1_print()
{
  string fname1[] = {"dndy_p.dat","dndy_ap.dat",
                     "dndy_pi-.dat","dndy_pi0.dat","dndy_pi+.dat",
                     "dndy_N.dat","dndy_pion.dat","dndy_aN.dat"};
  string fname2[] = {"v1_p.dat","v1_ap.dat",
      "v1_pi-.dat","v1_pi0.dat","v1_pi+.dat",
      "v1_N.dat","v1_pion.dat","v1_aN.dat"};
  string fname3[] = {"v2_p.dat","v2_ap.dat",
      "v2_pi-.dat","v2_pi0.dat","v2_pi+.dat",
      "v2_N.dat","v2_pion.dat","v2_antiN.dat"};

  double w=1.0/nEvent;
  for (int i=0;i<nbook;i++) {
    if(optHist==1) {
      Book1* h1 = new Book1(title2[i],ny,ymin,ymax);
      h1->operation(&v1[i],"/",&dndy[i]);
      h1->print(fname2[i],0);
      Book1* h2 = new Book1(title3[i],ny,ymin,ymax);
      h2->operation(&v2[i],"/",&dndy[i]);
      h2->print(fname3[i],0);
      delete h1,h2;
    } else {
      v1[i].scale(w);
      v2[i].scale(w);
      v1[i].print(fname2[i],0);
      v2[i].print(fname3[i],0);
    }

    dndy[i].scale(w);
    dndy[i].print(fname1[i],0);
  }

}

void book1_end()
{
  delete [] dndy;
  delete [] v1;
  delete [] v2;
}

double getY(double px, double py, double pz,double m)
{
    const double TINY = 1e-20;
    double mt=sqrt(m*m + px*px + py*py);
    double e=sqrt(mt*mt+pz*pz);
    double temp = log( (e + abs(pz)) / max( TINY, mt ) );
    return (pz > 0) ? temp : -temp;
}


bool acceptance(int kf, double pt)
{
    //return true;
    if(abs(kf)== 2212) {
	if(pt<0.4) return false;
	if(pt>2.0) return false;
    } else if(abs(kf)==221 || kf==111) {
	if(pt<0.2) false;
	if(pt>1.6) false;
    }
    return true;
}

void analyze()
{
    nEvent++;

    double dn[nbook][100];
    double dn2[nbook][100];

    for(int i=0;i<nbook;i++)
    for(int j=0;j<100;j++) {
	dn[i][j]=0.0;
	dn2[i][j]=0.0;
    }

    for(int i=0;i<nv;i++) {
	int id=-1;
	if(kf[i] == 2212) id=0;
	if(kf[i] == -2212)id=1;
	if(kf[i] == -211) id=2;
	if(kf[i] == 111)  id=3;
	if(kf[i] == 211)  id=4;
	int id2=-1;
	if(kf[i] == 2212 || kf[i]==2112) id2=5;
	if(kf[i] == 111 || abs(kf[i])== 211) id2=6;
	if(kf[i] == -2212 || kf[i]==-2112) id2=7;
	if(id >=0) {
	    double y=getY(px[i],py[i],pz[i],mass[i]);
	    dndy[id].fill(y,1.0/wy);
	    dndy[id2].fill(y,0.5/wy);
	    dndy[id2].fill(-y,0.5/wy);
	    double pt=sqrt(px[i]*px[i]+py[i]*py[i]);
            if(!acceptance(kf[i],pt)) continue;

	    int iy=(int)floor((y-ymin)/wy);
	    int iy2=(int)floor((-y-ymin)/wy);

	    if(iy>=0 && iy<ny) {

		if(optHist==1) {

	      double cos1=px[i]/pt;
	      double cos2=(px[i]*px[i] - py[i]*py[i])/(pt*pt);
	      v1[id].fill(y,cos1/wy);
	      v2[id].fill(y,cos2/wy);

	      v1[id2].fill(y,cos1/wy/2);
	      v2[id2].fill(y,cos2/wy/2);
	      v1[id2].fill(-y,-cos1/wy/2);
	      v2[id2].fill(-y,cos2/wy/2);


		} else {

		dn[id][iy] += 1.0;
		if(id2==-1) {
		    cout << " id2? " << id2 << endl;
		    exit(1);
		}
		dn[id2][iy] += 1.0;
		dn[id2][iy2] += 1.0;
		}

	    }
	}
    }

    if(optHist==1) return;

    for(int i=0;i<nv;i++) {
	int id=-1;
	if(kf[i]==2212) id=0;
	if(kf[i]==-2212)id=1;
	if(kf[i]==-211) id=2;
	if(kf[i]==111)  id=3;
	if(kf[i]==211)  id=4;
	int id2=-1;
	if(kf[i] == 2212 || kf[i]==2112) id2=5;
	if(kf[i] == 111 || abs(kf[i])== 211) id2=6;
	if(kf[i] == -2212 || kf[i]==-2112) id2=7;
	if(id >= 0) {
	    double ptsq=px[i]*px[i]+py[i]*py[i];
	    double pt=sqrt(ptsq);
            if(!acceptance(kf[i],pt)) continue;
	    double y=getY(px[i],py[i],pz[i],mass[i]);
	    int iy=(int)floor((y-ymin)/wy);
	    int iy2=(int)floor((-y-ymin)/wy);
	    if(iy>=0 && iy<ny) {
	      double cos1=px[i]/pt;
	      double cos2=(px[i]*px[i] - py[i]*py[i])/ptsq;
	      v1[id].fill(y,cos1/dn[id][iy]);
	      v2[id].fill(y,cos2/dn[id][iy]);
	      v1[id2].fill(y,cos1/dn[id2][iy]);
	      v2[id2].fill(y,cos2/dn[id2][iy]);
	      v1[id2].fill(-y,-cos1/dn[id2][iy2]);
	      v2[id2].fill(-y,cos2/dn[id2][iy2]);
	    }
	}
    }

}

int getParam(string templine,int ip)
{
    string com;
    istringstream is(templine);


    is >> kf[ip] >> mass[ip] >> px[ip] >> py[ip] >> pz[ip]>>pe[ip];
	 //>> vx[ip] >> vy[ip] >> vz[ip] >> ve[ip] >> vm[ip];

    /*
    cout << ip 
	 << " kf= " << kf[ip]
	 << " m= " << mass[ip]
	 << " px= " << px[ip]
	 << endl;
	 */

 
}

void read(string fname)
{
    string com;
    ifstream in;
    in.open(fname.c_str(),ios::in);
    if(!in) {
        cerr << "Error unable to open file " << fname << endl;
        exit(1);
    }
    string templine;
    int lineposition=0;

    getline(in,templine);
    istringstream is(templine);
    int nev,frame;
    double ylab, beta,gamma;
    is >> com >>  nev >> ylab >> beta >> gamma >> frame;

    cout << com << " nev= " << nev << " ylab= " << ylab << endl;

//  while(!in.eof()) {
    while(getline(in,templine)) {
        //in.getline(buffer,sizeof(buffer),'\n');
        lineposition++;
        int icomment1=templine.find('#');
	int ip;
        if(icomment1>=0) {
          istringstream is(templine);
          is >> com >>  iev >> nv >> nbary >> nmeson >> impactParam;
	  /*
	   cout << " icom= " << icomment1
	       << " com= " << com
	       << " iev= " << iev 
	       << " nv= " << nv
	       << " nbary= " << nbary
	       << endl;
	   cin.get();
	   */
	   ip=0;
	   if(nEvent % 500 == 0) cout << " nevent = " << nEvent << endl;
         } else {
           getParam(templine,ip);
	   ip++;
           if(ip==nv)  analyze();
	 }
    }

    in.close();
}


int main(int argc, char* argv[]) {

    int nfile=4;
    string fname="phase.dat";

    for(int i=1; i<argc; i++) {
      if(!strcmp(argv[i],"-o")) optHist = atoi(argv[i+1]);
      if(!strcmp(argv[i],"-n")) nfile = atoi(argv[i+1]);
      if(!strcmp(argv[i],"-f")) fname = argv[i+1];
    }

    kf = new int [ntot];
    mass = new double [ntot];
    px = new double [ntot];
    py = new double [ntot];
    pz = new double [ntot];
    pe = new double [ntot];
    //vx = new double [ntot];
    //vy = new double [ntot];
    //vz = new double [ntot];
    //ve = new double [ntot];
    //vm = new double [ntot];

    book1();

    stringstream ss;
    if(nfile==1) {
	read(fname);
    } else {
    for(int i=1;i<=nfile;i++) {
	ss << i;
	read(ss.str()+"/"+fname);
	ss.str("");
	ss.clear(stringstream::goodbit);
    }
    }

    book1_print();
    book1_end();
    cout << " optHist = " << optHist << endl;
    cout << " total event = " << nEvent << endl;

    delete [] kf;
    delete [] mass;
    delete [] px,py,pz,pe;
    //delete [] vx,vy,vz,ve,vm;
    return 0;
}

