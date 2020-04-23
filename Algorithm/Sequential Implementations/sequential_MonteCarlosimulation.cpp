#include<iostream>
#include <random>
#include <vector>
#include <stack>
#include <cmath>
using namespace std;

constexpr double PI = 3.14159265358979323846264338;
typedef unsigned long long ll;

// Struct to hold positions
struct Point {
    double x,y;
    Point(){}
    Point(double x, double y){
        this->x = x;
        this->y = y;
    }
};

class BSDlcg_random_generator {
public:
    BSDlcg_random_generator(){}
    BSDlcg_random_generator(ll seed){
        this->seed_x = seed%m;
        this->seed_y = (~seed)%m;
    }
    Point get_random_point(){
        double x = (double)seed_x/(double)mod;
        seed_x = (a*seed_x+c)%m;
        double y = (double)seed_y/(double)mod;
        seed_y = (a*seed_y+c)%m;
        return Point(x,y);
    }
private:
    ll m = (1ll<<31);
	ll a = 1103515245;
	ll c = 12345;
    ll mod = m-1;
    ll seed_x, seed_y;
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
    double a,b,c;
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
        double sum_x = 0;
        double sum_y = 0;
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
        double sum_x = 0;
        double sum_y = 0;
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
    double area(){
        double area = 0;
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
    inline double area_triangle(Point a, Point b, Point c){
        double area = 0.5 * ((b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y));
        return (area > 0.0) ? area : -area;
    }
};

// ************************************************************************* //

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

// A global point needed for sorting points with reference to the first point
Point base;

// A utility function to find next to top in a stack
Point nextToTop(stack<Point> &S){
	Point p = S.top();
	S.pop();
	Point res = S.top();
	S.push(p);
	return res;
}

double distSq(Point p1, Point p2){
	return (p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y);
}

int orientation(Point p, Point q, Point r){
	double val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
    if(abs(val) < 0.0000000000005) return 0;
	return (val > 0)? 1: 2; // clock or counterclock wise
}

// A function used by library function qsort() to sort an array of
// points with respect to the first point
int compare(const void *vp1, const void *vp2){
    Point *p1 = (Point *)vp1;
    Point *p2 = (Point *)vp2;

    // Find orientation
    int o = orientation(base, *p1, *p2);
    if (o == 0)
    	return (distSq(base, *p2) >= distSq(base, *p1))? -1 : 1;

    return (o == 2)? -1: 1;
}

// Prints convex hull of a set of n points.
std::vector<Point> convexHull(Point points[], int n){
    // Find the bottommost point
    double ymin = points[0].y;
    int min = 0;
    for (int i = 1; i < n; i++)
    {
    	double y = points[i].y;

    	// Pick the bottom-most or chose the left
    	// most point in case of tie
    	if ((y < ymin) || (ymin == y &&
    		points[i].x < points[min].x))
    		ymin = points[i].y, min = i;
    }

    // Place the bottom-most point at first position
    base = points[min];
	points[min] = points[0];
	points[0] = base;
    qsort(&points[1], n-1, sizeof(Point), compare);

    int m = 1; // Initialize size of modified array
    for (int i=1; i<n; i++)
    {
    	while (i < n-1 && orientation(base, points[i],points[i+1]) == 0)
    		i++;
    	points[m] = points[i];
    	m++; // Update size of modified array
    }

    std::vector<Point> vertices;
    if (m < 3) return vertices;

    stack<Point> S;
    S.push(points[0]);
    S.push(points[1]);
    S.push(points[2]);
    for (int i = 3; i < m; i++)
    {
    	while (orientation(nextToTop(S), S.top(), points[i]) != 2)
            S.pop();
    	S.push(points[i]);
    }
    while (!S.empty()){
    	Point p = S.top();
        vertices.push_back(p);
    	S.pop();
    }
    return vertices;
}

// ************************************************************************* //

class MonteCarlo {
public:
    MonteCarlo(unsigned seed){
        srand(seed);
        const int mod = 1000;
        // Point p = Point((double)(rand()%mod)/mod, (double)(rand()%mod)/mod);
        // Point q = Point((double)(rand()%mod)/mod, (double)(rand()%mod)/mod);
        // Point r = Point((double)(rand()%mod)/mod, (double)(rand()%mod)/mod);

        Point p = Point((double)1, (double)0);
        Point q = Point((double)0, (double)1);
        Point r = Point((double)0, (double)-1);

        ellipse = Ellipse(p,q,r);
        anatical_area = ellipse.area();
    }
    MonteCarlo(int n, unsigned seed) {
        // static std::default_random_engine generator;
        // static std::uniform_real_distribution<double> dist(0, 1);
        srand(seed);
        const int mod = 100;
        Point points[1024];
        for(int i=0;i<n;i++){
            points[i] = Point((double)(rand()%mod)/mod,(double)(rand()%mod)/mod);
        }
        std::vector<Point> vertices = convexHull(points,n);
        // std::cout<<"Print polygon used\n";
        // for(auto i:vertices){
        //     std::cout<<'('<<i.x<<','<<i.y<<")\n";
        // }

        shape = Polygon(vertices);
        anatical_area = shape.area();
    }
    double monte_carlo_pi(unsigned samples) {
        unsigned count = 0;
        for (unsigned i = 0; i < samples; ++i) {
            Point sample = random_generator.get_random_point();
            if (in_circle(sample))
                ++count;
        }

        return 4.0 * count / samples;
    }
    double run_ellipse_experiment(unsigned samples) {
        static std::default_random_engine generator;
        static std::uniform_real_distribution<double> dist(0, 1);
        unsigned hit_count = 0;

        for (int i = 0; i < samples; i++) {
            Point sample = random_generator.get_random_point();
            if(ellipse.is_interior(sample))
                hit_count++;
        }

        this->hit_count += hit_count;
        this->total_samples += samples;
        double calulated_area = 4.0 * calulate_area();
        std::cout<<"Actual area of polygon: "<<ellipse.area()<<'\n';
        std::cout<<"After "<<total_samples<<" samples estimated area is: "<<calulated_area<<'\n';
        std::cout<<"Percentage error: "<<percent_error(calulated_area)<<"%\n";
        return calulated_area;
    }
    double run_experiment(unsigned samples) {
        static std::default_random_engine generator;
        static std::uniform_real_distribution<double> dist(0, 1);
        unsigned hit_count = 0;

        for (int i = 0; i < samples; i++) {
            Point sample = random_generator.get_random_point();
            if(shape.is_interior(sample))
                hit_count++;
        }

        this->hit_count += hit_count;
        this->total_samples += samples;
        double calulated_area = calulate_area();
        std::cout<<"Actual area of polygon: "<<shape.area()<<'\n';
        std::cout<<"After "<<total_samples<<" samples estimated area is: "<<calulated_area<<'\n';
        std::cout<<"Percentage error: "<<percent_error(calulated_area)<<"%\n";
        return calulated_area;
    }
private:
    Polygon shape;
    Ellipse ellipse;
    double anatical_area;
    unsigned int hit_count = 0;
    unsigned int total_samples = 0;
    BSDlcg_random_generator random_generator;

    // Function to check whether point is in circle or no;
    inline bool in_circle(Point p) {
        return p.x * p.x + p.y * p.y < 1 * 1;
    }
    // Estimate error:
    inline double percent_error(double calulated_area) {
        double error_percent = 100 * (1 - calulated_area/anatical_area);
        return (error_percent > 0.0) ? error_percent : -error_percent;
    }
    inline double calulate_area() {
        return (double)hit_count/(double)total_samples;
    }
};

int main()
{
    // int a,b,c,d;
    unsigned samples, points, seed=0;

    // std::cout << "Enter a number of vertices to be used: ";
    // std::cin>> points;
    std::cout << "Enter a random number for seed: ";
    std::cin>> seed;
    //
    // MonteCarlo monte_carlo(points,seed);
    //
    // std::cout << "Enter samples to use: ";
    // std::cin >> samples;
    //
    // double pi_estimate = monte_carlo.monte_carlo_pi(samples);
    // std::cout << "Pi = " << pi_estimate << '\n';
    // std::cout << "Percent error is: " << 100 * std::abs(pi_estimate - PI) / PI << " %\n";
    //
    // monte_carlo.run_experiment(samples);

    MonteCarlo ellipse_monte_carlo(seed);
    std::cout << "Enter samples to use: ";
    std::cin >> samples;
    ellipse_monte_carlo.run_ellipse_experiment(samples);

    return 0;
}
