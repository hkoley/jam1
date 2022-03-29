#ifndef Book1_h
#define Book1_h

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <cmath>

using namespace std;

class Book1
{
private:
    string Title;
    int    numberOfBin;
    double xMin;
    double xMax;
    int    counter;
    double *bin;
    double *w2;
    int    *hit;
    double xLow;
    double xUpper;
    double binSize;

public:
    Book1() {};
    Book1(string title,int nx,double xl,double xu) {
	setPar(title,nx,xl,xu);
    }
    ~Book1() {
	delete [] bin;
	delete [] hit;
	if(w2) delete [] w2;
    }

    void   setPar(string title,int nx,double xl,double xu);
    double getBinWidth() {return binSize;}
    double getBin(int i) {return bin[i];}
    double getW2(int i) {return w2[i];}
    int    getHit(int i) {return hit[i];}
    double getError(int i) {
	return (hit[i] > 0) ? bin[i]/sqrt((double)hit[i]) : 0.0;}
    void   fill(double x,double w);
    void   fill(int ix,double w);
    void   scale(double w);
    void   print(string title, int norm);
    void   operation(Book1* hist1,string oper, Book1* hist2);
    void   setW2();
};




#endif // Book1_h

