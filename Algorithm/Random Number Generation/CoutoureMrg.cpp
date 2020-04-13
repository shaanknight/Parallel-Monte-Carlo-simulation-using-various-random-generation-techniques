#include <stdio.h>
#include <string.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
using namespace std;

typedef long long ll;
typedef double db;

vector<db> randomgen_CoutureMRG(ll n, vector<ll> seed)
{
    ll a1 = 1071064;
    ll a7 = 2113664;
    ll m = (1ll << 31) - 19;
    static vector<ll> y(n);
    static vector<db> rnd_seq(n,0);
    for (int i = 0; i < 7; i++)
    {
        y[i] = seed[i];
        rnd_seq[i] = db(y[i]) / db(m-1.0);
    }
    for (ll i = 7; i < n; i++)
    {
        y[i] = ( a1 * y[i-1] + a7 * y[i-7] ) % m;
        rnd_seq[i] = db(y[i]) / db(m-1.0);
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

vector<ll> seedgen_Couture(ll rank,ll iters,vector<ll> seed)
{
    vector<vector<ll> > pre((ll)seed.size());
    for(ll i=0;i<(ll)seed.size();++i)
        pre[i].resize((ll)seed.size(),0);
    for(ll i=0;i+1<(ll)seed.size();++i)
        pre[i][i+1] = 1;
    ll l = seed.size();
    pre[l-1][0] = 2113664;
    pre[l-1][l-1] = 1071064;
    return multiply(pre,seed,rank*iters,(1ll << 31) - 19);
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

void trace(vector<db> td){
    cout << (ll) td.size() << "\n";
    for(ll i=0;i<(ll) td.size();++i)
        cout << td[i] << " ";
    cout << "\n";
}

void tracell(vector<ll> td)
{
    cout << (ll) td.size() << "\n";
    for(ll i=0;i<(ll) td.size();++i)
        cout << td[i]/((1ll << 31) - 20.0) << " ";
    cout << "\n";
}

int main(){
    // vector<ll> seed(7);
    // for (int i = 0; i < 7; i++)
    //     seed[i] = 12345;
    // uniformity(randomgen_CoutureMRG(10000000,seed));
    vector<ll> seedx(7),seedy(7);
    for (ll i = 0; i < 7; i++)
    {
        seedx[i] = 12345;
        seedy[i] = 54321;
    }
    // seedx = seedgen_Couture(rank,samples,seedx);
    // seedy = seedgen_Couture(rank,samples,seedy);
    trace(randomgen_CoutureMRG(10000001,seedx));
    tracell(seedgen_Couture(1,10000000,seedx));
}