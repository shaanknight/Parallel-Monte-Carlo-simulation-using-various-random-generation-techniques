#include<iostream>
#include <random>
#include <vector>
#include <stack>
using namespace std;

typedef long long ll;
typedef double db;
constexpr db PI = 3.14159265358979323846264338;

// Struct to hold positions
struct Point {
    db x,y;
    Point(){}
    Point(db x, db y){
        this->x = x;
        this->y = y;
    }
};

// ************************************************************************* //
/* analytical finding of area of polygon to measure the error of monte carlo's method of finding area */

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
        anatical_area = shape.area();
    }
    ll sampling(ll n,vector<ll> seedx,vector<ll> seedy)
    {
        ll m = (1ll<<31);
        ll a = 1103515245;
        ll c = 12345;
        ll x = seedx[0]%m;
        ll y = seedy[0]%m;
        ll count = 0;
        Point sample = Point(x/(m-1.0),y/(m-1.0));
        if (shape.is_interior(sample))
            ++count;
        for(ll i=1;i<n;++i)
        {
            x = (a*x+c)%m;
            y = (a*y+c)%m;
            sample = Point(x/(m-1.0),y/(m-1.0));
            if(shape.is_interior(sample))
                ++count;
        }
        return count;
    }

    db run_experiment(unsigned samples) {
        this->hit_count += sampling(samples,{12345},{54321});
        this->total_samples += samples;
        db calulated_area = calulate_area();
        std::cout<<"Actual area of polygon: "<<shape.area()<<'\n';
        std::cout<<"After "<<total_samples<<" samples estimated area is: "<<calulated_area<<'\n';
        std::cout<<"Percentage error: "<<percent_error(calulated_area)<<"%\n";
        return calulated_area;
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

int main()
{
    vector<Point> vertices;
    db px,py;
    while(cin >> px)
    {
        cin >> py;
        vertices.push_back(Point(px,py));
    }

    MonteCarlo monte_carlo(vertices);
    monte_carlo.run_experiment(10000000);

    return 0;
}
