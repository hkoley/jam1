#include "Book1.h"
#include <algorithm>
#include <cstdio>

using namespace std;


#ifdef ALPHA
#define abs(a) (((a) > 0) ? (a) : -(a))
#endif

void Book1::setPar(string title,int nx,double xl,double xu)
{
	Title = title;
	numberOfBin=nx;
	xMin=xl;
	xMax=xu;
	counter = 0;
	xLow = 0.0;
	xUpper = 0.0;
	binSize = (xMax-xMin)/numberOfBin;
	if(xMax <= xMin ) {
	    cerr << "(Book1:) Are you sure? xmin = " << xMin
	       << " xMax=" << xMax << endl;
	    exit(1);
	}
	w2 = 0;
	bin = new double[numberOfBin];
	hit = new int[numberOfBin];
	for(int i=0; i<numberOfBin;i++) {
	    bin[i]=0.0;
	    hit[i]=0;
	}
	setW2();
}

void Book1::setW2()
{
    w2 = new double[numberOfBin];
    for(int i=0; i<numberOfBin;i++) {
	w2[i]=0.0;
    }
}

//...Find bin in i, including under/overflow, and fill.
void Book1::fill(int ix, double w)
{

    if(ix < 0) {
	xLow += w;
    }else if(ix >= numberOfBin) {
	xUpper += w;
    } else {
        bin[ix] += w;
        hit[ix]++;
	if(w2) w2[ix] += w*w;
    }
}

//...Purpose: to accumulate statistics and fill 1-dim. histogram.
void Book1::fill(double x, double w)
{
    ++counter;

//...Find bin in x, including under/overflow, and fill.
    int ix = (int) floor((x-xMin)/binSize);
    if(ix < 0) {
	xLow += w;
    }else if(ix >= numberOfBin) {
	xUpper += w;
    } else {
        bin[ix] += w;
        hit[ix]++;
	if(w2) w2[ix] += w*w;
    }

}

// Rescale histgram.
void Book1::scale(double w)
{
    for(int i=0;i<numberOfBin;i++) {
	bin[i] *= w;
	if(w2) w2[i] *= w*w;
    }
}

// Print out accumulated statistics.
// norm=1: Normalize histgram
void Book1::print(string outFile,int norm)
{

    ofstream ofs(outFile.c_str());

    // Normalization.
    double fev=1.0;

    ofs << "# " <<  Title.c_str() << endl;
    ofs << "# bin size " << numberOfBin
       	<< " size = " << binSize << endl;
    ofs << "#" << "   x   " << "         y(x)    ";
    if(w2) ofs << "        error1    " << "      error2" << endl;
    else
       ofs << "       error    " << endl;
    ofs << "#" << endl;

    if(norm) {
	fev=0.0;
	for(int ix=0; ix<numberOfBin; ix++) {
	    fev += bin[ix];
	}
	fev *= binSize;
	if(fabs(fev) > 1e-30) {
	    fev = 1.0/fev;
	}else {
	    fev = 1.0;
	}
    }

    double ys=0.0;
    char str[1024];
    for(int ix=0; ix<numberOfBin; ix++) {
        double xx= xMin + (ix+0.5)*binSize;
        double yx=bin[ix];
        ys += yx;
        double err = (hit[ix] > 0) ? fev*fabs(yx)/sqrt((double)hit[ix]) : 0.0;
	if(w2) {
	    sprintf(str,"   %g     %g     %g     %g",
		    xx,fev*yx,fev*sqrt(w2[ix]),err);
	}else{
	    sprintf(str,"   %g     %g     %g",xx,fev*yx,err);
	}
	ofs << str << '\n';

    }
        //if(mform.eq.0) write(iunit,1600) fac*ys
        //if(mform.eq.1) write(iunit,1800) fac*ys
    
    ofs.close();

}


void Book1::operation(Book1* hist1,std::string oper, Book1* hist2)
{

    //if(!strcmp(oper.c_str(),"+")) {
    if(oper == "+") {
        for(int i=0;i<numberOfBin;i++) {
	    bin[i] = hist1->getBin(i) + hist2->getBin(i);
	}
    //} else if(!strcmp(oper.c_str(),"-")) {
    } else if (oper == "-") {
        for(int i=0;i<numberOfBin;i++) {
	    bin[i] = hist1->getBin(i) - hist2->getBin(i);
	}
    //} else if(!strcmp(oper.c_str(),"*")) {
    } else if (oper == "*") {
        for(int i=0;i<numberOfBin;i++) {
	    bin[i] = hist1->getBin(i) * hist2->getBin(i);
	}
    //} else if(!strcmp(oper.c_str(),"/")) {
    } else if (oper == "/") {
        for(int i=0;i<numberOfBin;i++) {
	    double h2 = hist2->getBin(i);
	    if(h2 > 1e-10) {
		double h1 = hist1->getBin(i);
		bin[i] = h1 / h2;
		hit[i] = hist1->getHit(i);
		if(w2) {
		    double e1 = sqrt(hist1->getW2(i));
		    double e2 = sqrt(hist2->getW2(i));
		    w2[i] = (e1*e1*h2*h2+e2*e2*h1*h1)/(h2*h2*h2*h2);

		}
	    } else {
		bin[i]=0.0;
		hit[i] = 0;
	    }
	}
    }
}

