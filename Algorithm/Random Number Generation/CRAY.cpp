#include <stdio.h>
#include <string.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
using namespace std;

// Lehmar's LCG used in CRAY cannot be run on our 64 bit machines
// bad attempt

typedef long long ll;

vector<ll> cray_randomgen(ll n,ll seed)
{
	ll m = (1ll<<48);
	ll a = 44485709377909;
	
	seed = seed%m;
	static vector<ll> rnd_seq;
	rnd_seq.push_back(seed);
	for(ll i=1;i<n;++i)
	{
		seed = a*seed%m;
		rnd_seq.push_back(seed);
	}

	return rnd_seq;
}

void trace(vector<ll> td)
{
    cout << (ll) td.size() << "\n";
    for(ll i=0;i<(ll) td.size();++i)
        cout << td[i] << " ";
    cout << "\n";
}

int main()
{
	trace(cray_randomgen(100,1));
}
