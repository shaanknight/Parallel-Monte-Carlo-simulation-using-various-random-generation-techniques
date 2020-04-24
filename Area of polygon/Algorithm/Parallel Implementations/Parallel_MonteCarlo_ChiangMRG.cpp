/* MPI Program Template */

#include <stdio.h>
#include <string.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <queue>
#include "mpi.h"
using namespace std;

typedef long long ll;
typedef double db;

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

vector<ll> seedgen_Chiang(ll rank,ll iters,vector<ll> seed)
{
    vector<vector<ll> > pre((ll)seed.size());
    for(ll i=0;i<(ll)seed.size();++i)
        pre[i].resize((ll)seed.size(),0);
    for(ll i=0;i+1<(ll)seed.size();++i)
        pre[i][i+1] = 1;
    ll l = seed.size();
    pre[l-1][0] = 598295599;
    pre[l-1][1] = 1546795921;
    pre[l-1][2] = 885300443;
    return multiply(pre,seed,rank*iters,(1ll << 31) - 1);
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

// ax + by = c
class Line {
public:
    Line(){}
    // Form a line given two points
    Line(const Point p, const Point q){
        a = p.y - q.y;
        b = q.x - p.x;
        c = -a * p.x -b * p.y;
    }
    // Tell if a point lies on the positive side of a line
    bool is_positive(Point p){
        return a * p.x + b * p.y + c >= 0;
    }
    // Ensure that a point p lies on the positive side of the line
    void flip(Point p){
        if(!is_positive(p)){
            a = -a;
            b = -b;
            c = -c;
        }
    }
    void print(){
        std::cout <<a<<"x + "<<b<<"y + "<<c<<" = 0\n";
    }
private:
    db a,b,c;
};

class Polygon {
public:
    Polygon(int n=4){
        this->n = n;
        vertices.push_back(Point(0,0));
        vertices.push_back(Point(0.5,0));
        vertices.push_back(Point(0.5,0.5));
        vertices.push_back(Point(0,0.5));

        // Calculate centroid
        db sum_x = 0;
        db sum_y = 0;
        for(int i=0;i<n;i++){
            sum_x += vertices[i].x;
            sum_y += vertices[i].y;
        }
        centroid = Point(sum_x/n, sum_y/n);

        // Set eqation of n sides
        for(int i=0;i<n-1;i++){
            sides.push_back(Line(vertices[i],vertices[i+1]));
            sides[i].flip(centroid);
        }
        sides.push_back(Line(vertices[n-1],vertices[0]));
        sides[n-1].flip(centroid);
    }
    Polygon(std::vector<Point> vertices){
        n = vertices.size();
        this->vertices = vertices;
        // Calculate centroid
        db sum_x = 0;
        db sum_y = 0;
        for(int i=0;i<n;i++){
            sum_x += vertices[i].x;
            sum_y += vertices[i].y;
        }
        centroid = Point(sum_x/n, sum_y/n);
        // Set eqation of n sides
        for(int i=0;i<n-1;i++){
            sides.push_back(Line(vertices[i],vertices[i+1]));
            sides[i].flip(centroid);
        }
        sides.push_back(Line(vertices[n-1],vertices[0]));
        // Control the side
        sides[n-1].flip(centroid);
    }
    // Function to check whether point is in polygon or not
    bool is_interior(Point p){
        // O(n) algorithm
        for(int i=0;i<n;i++)
            if(!sides[i].is_positive(p))
                return false;
        return true;
    }
    // Analytical method to calculate area of the polygon
    db area(){
        db area = 0;
        for(int i=0;i<n-1;i++)
            area += area_triangle(centroid, vertices[i], vertices[i+1]);
        area += area_triangle(centroid, vertices[n-1], vertices[0]);
        return area;
    }
private:
    unsigned int n;
    std::vector<Point> vertices;
    std::vector<Line> sides;
    Point centroid;
    // Function to calculate area of triangle
    inline db area_triangle(Point a, Point b, Point c){
        db area = 0.5 * ((b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y));
        return (area > 0.0) ? area : -area;
    }
};

// ************************************************************************* //
/* MonteCarlo simulation */

class MonteCarlo {
public:
    MonteCarlo(vector<Point> vertices) {
        shape = Polygon(vertices);
    }
    ll sampling(ll n, vector<ll> seedx,vector<ll> seedy)
    {
        ll a1 = 885300443;
        ll a2 = 1546795921;
        ll a3 = 598295599;
        ll m = (1ll << 31) - 1;
        ll count = 0;
        queue<ll> x,y;
        ll lx,ly;
        Point sample;
        for (ll i = 0; i < 3; i++)
        {
            x.push(seedx[i]);
            y.push(seedy[i]);
            lx = seedx[i];
            ly = seedy[i];
            sample = Point(lx/(m-1.0),ly/(m-1.0));
            if(shape.is_interior(sample))
                ++count;
        }
        // cout<<sample.x<<" "<<sample.y<<endl;
        for (ll i = 3; i < n; i++)
        {
            ll p3x = x.front();
            ll p3y = y.front();
            x.pop(), y.pop();
            ll p2x = x.front();
            ll p2y = y.front();
            x.pop(), y.pop();
            ll p1x = x.front();
            ll p1y = y.front();
            x.pop(), y.pop();
            ll p0x = ( a1*p1x + a2*p2x + a3*p3x ) % m;
            ll p0y = ( a1*p1y + a2*p2y + a3*p3y ) % m;
            x.push(p2x), x.push(p1x), x.push(p0x);
            y.push(p2y), y.push(p1y), y.push(p0y);
            sample = Point(p0x/(m-1.0),p0y/(m-1.0));
            // cout<<sample.x<<" "<<sample.y<<endl;
            if(shape.is_interior(sample))
                ++count;
        }
        return count;
    }
    ll run_experiment(unsigned samples,int rank) {
        vector<ll> seedx(3),seedy(3);
        seedx[0] = 2858447;
        seedx[1] = 12345;
        seedx[2] = 45991;
        seedy[0] = 1000007;
        seedy[1] = 54321;
        seedy[2] = 8295599;
        vector<ll> resx = seedgen_Chiang(rank,samples,seedx);
        vector<ll> resy = seedgen_Chiang(rank,samples,seedy);
        return sampling(samples,resx,resy);
    }
    void analyse(ll hits,ll samples){
        anatical_area = shape.area();
        this->hit_count += hits;
        this->total_samples += samples;
        db calulated_area = calulate_area();
        std::cout<<"Actual area of polygon: "<<shape.area()<<'\n';
        std::cout<<"After "<<total_samples<<" samples estimated area is: "<<calulated_area<<'\n';
        std::cout<<"Percentage error: "<<percent_error(calulated_area)<<"%\n";
        return;
    }
private:
    Polygon shape;
    db anatical_area;
    unsigned int hit_count = 0;
    unsigned int total_samples = 0;
    
    // Estimate error:
    inline db percent_error(db calulated_area) {
        db error_percent = 100 * (1 - calulated_area/anatical_area);
        return (error_percent > 0.0) ? error_percent : -error_percent;
    }
    inline db calulate_area() {
        return (db)hit_count/(db)total_samples;
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