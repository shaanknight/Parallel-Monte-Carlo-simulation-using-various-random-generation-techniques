#include <stdio.h>
#include <string.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
using namespace std;

typedef long long ll;
typedef double db;

vector<db> randomgen_LEMRG(ll n, vector<ll> seed)
{
    ll a12 = 1403580;
    ll a13 = 810728;
    ll a21 = 527612;
    ll a23 = 1370589;
    ll m1 = (1ll << 32) - 209;
    ll m2 = (1ll << 32) - 22853;
    static vector<ll> y1(n), y2(n);
    static vector<db> rnd_seq(n);
    for (int i = 0; i < 3; i++)
    {
        y1[i] = seed[i];
        y2[i] = seed[i+3];
        rnd_seq[i] = (db)( ( y1[i] + y2[i]) % m1 ) / (db)(m1-1.0);
    }
    for (ll i = 3; i < n; i++){
        y1[i] = ( (a12 * y1[i-2])%m1 - (a13 * y1[i-3])%m1 + m1) % m1;
        y2[i] = ( (a21 * y2[i-1])%m2 - (a23 * y2[i-3])%m2 + m2) % m2;
        rnd_seq[i] = (db)( ( y1[i] + y2[i]) % m1 ) / (db)(m1-1.0);
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
                    tmp[i] = (tmp[i] + pre[i][j]*base[j])%mod; 
            for(i=0;i<N;++i)
                base[i] = tmp[i];
        }
        for(i=0;i<N;++i)
            for(j=0;j<N;++j)
                red[i][j] = 0;
        for(i=0;i<N;++i)
            for(j=0;j<N;++j)
                for(k=0;k<N;++k)
                {
                    if(pre[k][j] < 0 && pre[i][k] > 0)
                        red[i][j] = (red[i][j] - (pre[i][k]*(-1*pre[k][j]))%mod + mod)%mod;
                    else if(pre[i][k] < 0 && pre[k][j] > 0)
                        red[i][j] = (red[i][j] - (pre[k][j]*(-1*pre[i][k]))%mod + mod)%mod;
                    else
                        red[i][j] = (red[i][j] + 1ll*pre[i][k]*pre[k][j]%mod)%mod;
                }
        for(i=0;i<N;++i)
            for(j=0;j<N;++j)
                pre[i][j] = red[i][j];
        m >>= 1;
    }
    return base;
}

template<typename T>
vector<T> s(vector<T> const &v, int m, int n) 
{
   auto first = v.begin() + m;
   auto last = v.begin() + n + 1;
   vector<T> vector(first, last);
   return vector;
}

vector<ll> seedgen_Lecuyer(ll rank,ll iters,vector<ll> seed)
{
    ll a12 = 1403580;
    ll a13 = -810728;
    ll a21 = 527612;
    ll a23 = -1370589;
    ll m1 = (1ll << 32) - 209;
    ll m2 = (1ll << 32) - 22853;
    vector<ll> seed1 = s(seed,0,2);
    vector<ll> seed2 = s(seed,3,5);
    vector<ll> genseed1 = multiply({{0,1,0},{0,0,1},{a13,a12,0}},seed1,rank*iters,m1);
    vector<ll> genseed2 = multiply({{0,1,0},{0,0,1},{a23,0,a21}},seed2,rank*iters,m2);
    seed.clear();
    seed.insert(seed.begin(),genseed1.begin(),genseed1.end());
    seed.insert(seed.end(),genseed2.begin(),genseed2.end());
    return seed;
}

void uniformity(vector<db> td){
    // shows the distribution of numbers in the interval of 0.1
    vector<ll> cnt(10,0);
    for(ll i=0;i<(ll) td.size();++i){
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

void traceall(vector<ll> td){
    cout << (ll) td.size() << "\n";
    ll m1 = (1ll << 32) - 209;
    for(ll i=0;i<3;++i)
        cout << ((td[i]+td[i+3]+3*m1)%((1ll << 32) - 209))/((1ll << 32) - 208.0) << " ";
    cout << "\n";
}

int main(){
    vector<ll> seed(6);
    // First three terms of seed should be less than m1, and the last 3 values less than m2.
    for (int i = 0; i < 6; i++)
        seed[i] = 12345;
    uniformity(randomgen_LEMRG(10000000,seed));
    // trace(randomgen_LEMRG(7,seed));
    // traceall(seedgen_Lecuyer(1,6,seed));
}