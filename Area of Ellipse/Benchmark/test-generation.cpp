#include<bits/stdc++.h>
#include "testlib.h"
using namespace std;

int main(int argc, char* argv[])
{
	registerGen(argc, argv, 1);

	int ta = rnd.next(50000,100000);
	double a = ta/200000.0;
	double b = rnd.next(1,100000)/200000.0;

	for(int i=1;i<=3;++i)
	{
		double x = rnd.next(0,ta)/200000.0;
		double r = 1.0-(x*x)/(a*a);
		double y = b*sqrt(r);
		int sgn = rnd.next(0,3);
		if(sgn == 0)
			cout << "-" << x << " -" << y << endl;
		else if(sgn == 1)
			cout << "-" << x << " " << y << endl;
		else if(sgn == 2)
			cout << "" << x << " -" << y << endl;
		else
			cout << "" << x << " " << y << endl;
	}

	return 0;
}