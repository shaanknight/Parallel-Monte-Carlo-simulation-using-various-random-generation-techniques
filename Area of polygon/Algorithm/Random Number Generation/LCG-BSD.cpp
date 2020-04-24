#include <stdio.h>
#include <string.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
using namespace std;

typedef long long ll;
typedef double db;

vector<db> randomgen_BSDlcg(ll n,vector<ll> seed)
{
	ll m = (1ll<<31);
	ll a = 1103515245;
	ll c = 12345;
	ll y = seed[0]%m;
	static vector<db> rnd_seq;
	rnd_seq.push_back(y/(m-1.0));
	for(ll i=1;i<n;++i)
	{
		y = (a*y+c)%m;
		rnd_seq.push_back(y/(m-1.0));
	}

	return rnd_seq;
}

vector<ll> multiply(vector<vector<ll> > pre,vector<ll> base,ll m,ll mod)
{
	ll i,j,k;
	ll N = base.size();
	vector<ll> tmp(N);
	vector<vector<ll> > red(N);
	for(i=0;i<N;++i)
		red[i].resize(N,0);
	while(m)
	{
		if(m&1ll)
		{
			for(i=0;i<N;++i)
				tmp[i] = 0;
			for(i=0;i<N;++i)
				for(j=0;j<N;++j)
					tmp[i] = (tmp[i] + 1ll*pre[i][j]*base[j]%mod)%mod; 
			for(i=0;i<N;++i)
				base[i] = tmp[i];
		}
		for(i=0;i<N;++i)
			for(j=0;j<N;++j)
				red[i][j] = 0;
		for(i=0;i<N;++i)
			for(j=0;j<N;++j)
				for(k=0;k<N;++k)
					red[i][j] = (red[i][j] + 1ll*pre[i][k]*pre[k][j]%mod)%mod;
		for(i=0;i<N;++i)
			for(j=0;j<N;++j)
				pre[i][j] = red[i][j];
		m >>= 1;
	}
	return base;
}

vector<ll> seedgen_BSD(ll rank,ll iters,vector<ll> seed)
{
	return multiply({{1103515245,12345},{0,1}},seed,rank*iters,1ll<<31);
}

void uniformity(vector<db> td)
{
	// shows the distribution of numbers in the interval of 0.1
	vector<ll> cnt(10,0);
	for(ll i=0;i<(ll) td.size();++i)
	{
		int l = td[i]*10.0;
		cnt[l]++;
	}
	cout << "Uniformity of the sequence of size " << (ll) td.size() << " goes as : " << "\n";
	for(ll i=0;i<10;++i)
		cout << cnt[i] << " | ";
	cout << "\n";
}

void trace(vector<db> td)
{
    cout << (ll) td.size() << "\n";
    for(ll i=0;i<(ll) td.size();++i)
        cout << td[i] << " ";
    cout << "\n";
}

void tracell(vector<ll> td)
{
    cout << (ll) td.size() << "\n";
    for(ll i=0;i<(ll) td.size();++i)
        cout << td[i] << " ";
    cout << "\n";
}

int main()
{
	trace(randomgen_BSDlcg(11,{1,1}));
	tracell(seedgen_BSD(1,10,{1,1}));
}
