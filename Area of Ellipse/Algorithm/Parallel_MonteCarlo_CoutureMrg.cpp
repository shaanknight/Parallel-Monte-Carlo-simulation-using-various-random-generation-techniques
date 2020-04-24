/* MPI Program Template */

#include <stdio.h>
#include <string.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <queue>
#include "mpi.h"
using namespace std;

typedef long long ll;
typedef double db;

constexpr double PI = 3.14159265358979323846264338;

vector<db> vec,ar,rec,tmp;
string input;

// ************************************************************************* //
/* random generation of samples in parallel using fast exponentiation method */

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

// ************************************************************************* //
/* analytical finding of area of polygon to measure the error of monte carlo's method of finding area */

// Struct to hold positions
struct Point {
    db x,y;
    Point(){}
    Point(db x, db y){
        this->x = x;
        this->y = y;
    }
};

class Ellipse {
public:
    Ellipse(){}
    Ellipse(Point p, Point q, Point r){
        // Assuming d is handled
        double d = p.x*p.x * q.y*q.y - q.x*q.x * p.y*p.y;
        // If points are collinear with origin
        if(d==0)
        {
            Point temp = p;
            p = r;
            r = temp;
        }
        A = (q.y*q.y - p.y*p.y)/d;
        B = (p.x*p.x - q.x*q.x)/d;
    }
    // Function to check whether point is in ellipse or not
    bool is_interior(Point p){
        return A*p.x*p.x + B*p.y*p.y <= 1;
    }
    // Analytical method to calculate area of the ellipse
    double area(){
        double area = PI * std::sqrt(1/A) * std::sqrt(1/B);
        return area;
    }
private:
    // Assuming ellipse to be of form Ax^2 + By^2 = 1
    double A,B;
};

// ************************************************************************* //
/* MonteCarlo simulation */

class MonteCarlo {
public:
    MonteCarlo(vector<Point> vertices){
        ellipse = Ellipse(vertices[0],vertices[1],vertices[2]);
    }
    ll sampling(ll n, vector<ll> seedx,vector<ll> seedy)
    {
        ll a1 = 1071064;
        ll a7 = 2113664;
        ll m = (1ll << 31) - 19;
        ll count = 0;
        queue<ll> x,y;
        ll lx,ly;
        Point sample;
        for (ll i = 0; i < 7; i++)
        {
            x.push(seedx[i]);
            y.push(seedy[i]);
            lx = seedx[i];
            ly = seedy[i];
            sample = Point(lx/(m-1.0),ly/(m-1.0));
            if(ellipse.is_interior(sample))
                ++count;
        }
        for (ll i = 7; i < n; i++)
        {
            lx = ( a1 * lx + a7 * x.front() ) % m;
            ly = ( a1 * ly + a7 * y.front() ) % m;
            x.pop();
            y.pop();
            x.push(lx);
            y.push(ly);
            sample = Point(lx/(m-1.0),ly/(m-1.0));
            if(ellipse.is_interior(sample))
                ++count;
        }
        return count;
    }
    ll run_experiment(unsigned samples,int rank) {
        vector<ll> seedx(7),seedy(7);
        for (ll i = 0; i < 7; i++)
        {
            seedx[i] = 12345;
            seedy[i] = 54321;
        }
        vector<ll> resx = seedgen_Couture(rank,samples,seedx);
        vector<ll> resy = seedgen_Couture(rank,samples,seedy);
        return sampling(samples,resx,resy);
    }
    void analyse(ll hits,ll samples){
        anatical_area = ellipse.area();
        this->hit_count += hits;
        this->total_samples += samples;
        db calulated_area = calulate_area();
        std::cout<<"Actual area of ellipse: "<< anatical_area <<'\n';
        std::cout<<"After "<<total_samples<<" samples estimated area is: "<<calulated_area<<'\n';
        std::cout<<"Percentage error: "<<percent_error(calulated_area)<<"%\n";
        return;
    }
private:
    Ellipse ellipse;
    double anatical_area;
    unsigned int hit_count = 0;
    unsigned int total_samples = 0;

    // Estimate error:
    inline double percent_error(double calulated_area) {
        double error_percent = 100 * (1 - calulated_area/anatical_area);
        return (error_percent > 0.0) ? error_percent : -error_percent;
    }
    inline double calulate_area() {
        return 4*(double)hit_count/(double)total_samples;
    }
};

int main( int argc, char **argv ) {
    int rank, numprocs;

    /* start up MPI */
    MPI_Init( &argc, &argv );

    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &numprocs );

    /* user input from file */
    int n;
    if(!(rank))
    {
        input = argv[1];
        ifstream fin (input);
        db t;
        while(fin >> t)
            vec.push_back(t);
        n = vec.size();
        for(ll i=1;i<numprocs;++i)
        {
            for(ll j=0;j<n;++j)
            {
                vec.push_back(vec[j]);
            }
        }
        fin.close();
    }
    
    /*synchronize all processes*/
    MPI_Barrier( MPI_COMM_WORLD );
    double tbeg = MPI_Wtime();

    /* write your code here */
    MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
    ar.resize(n);
    MPI_Scatter(vec.data(),n,MPI_DOUBLE,ar.data(),n,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Status stat;

    ll samples = 10000000;
    vector<Point> vertices;
    for(ll i=0;i<(ll) ar.size();i+=2)
        vertices.push_back(Point(ar[i],ar[i+1]));
    MonteCarlo monte_carlo(vertices);
    ll hits = monte_carlo.run_experiment(samples/numprocs,rank);
    MPI_Allreduce(MPI_IN_PLACE,&hits,1,MPI_LONG_LONG_INT,MPI_SUM,MPI_COMM_WORLD);

    if(!rank)
        monte_carlo.analyse(hits,samples);

    MPI_Barrier( MPI_COMM_WORLD );
    double elapsedTime = MPI_Wtime() - tbeg;
    double maxTime;
    MPI_Reduce( &elapsedTime, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    if ( rank == 0 ) {
        printf( "Total time (s): %f\n", maxTime );
    }

    /* shut down MPI */
    MPI_Finalize();
    return 0;
}